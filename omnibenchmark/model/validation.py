"""Validation logic for benchmark models."""

import os
from pathlib import Path
from typing import Dict, List, TYPE_CHECKING, Any, Optional
from urllib.parse import urlparse

if TYPE_CHECKING:
    from .benchmark import SoftwareEnvironment, SoftwareBackendEnum


class ValidationError(Exception):
    """Exception raised for validation errors."""

    def __init__(self, *errors: Any) -> None:
        self.errors = errors[0] if len(errors) == 1 else errors
        super().__init__(*errors)

    def __str__(self) -> str:
        if isinstance(self.errors, str):
            return self.errors
        else:
            return "\n".join(map(str, self.errors))


class BenchmarkParseError(Exception):
    """Exception raised during benchmark parsing with location information.

    This exception carries additional context about where in the YAML file
    the error occurred, allowing the CLI layer to format it appropriately.
    """

    def __init__(
        self,
        message: str,
        yaml_file: Optional[Path] = None,
        line_number: Optional[int] = None,
        stage_id: Optional[str] = None,
        module_id: Optional[str] = None,
        parameter_index: Optional[int] = None,
        values: Optional[List[Any]] = None,
        original_error: Optional[Exception] = None,
    ):
        self.message = message
        self.yaml_file = yaml_file
        self.line_number = line_number
        self.stage_id = stage_id
        self.module_id = module_id
        self.parameter_index = parameter_index
        self.values = values
        self.original_error = original_error
        super().__init__(message)

    def __str__(self) -> str:
        """Format a compact, helpful error message with context."""
        parts = []

        # Add file and line number if available
        if self.yaml_file:
            location = str(self.yaml_file)
            if self.line_number:
                location += f":{self.line_number}"
            parts.append(f"Error in {location}")

        # Add stage/module context if available
        if self.stage_id or self.module_id:
            context_parts = []
            if self.stage_id:
                context_parts.append(f"stage '{self.stage_id}'")
            if self.module_id:
                context_parts.append(f"module '{self.module_id}'")
            if self.parameter_index is not None:
                context_parts.append(f"parameter #{self.parameter_index}")
            parts.append(f"  in {', '.join(context_parts)}")

        # Add the actual error message
        parts.append(f"  {self.message}")

        # Add parameter values if available (helps with debugging)
        if self.values:
            parts.append(f"  values: {self.values}")

        return "\n".join(parts)


def _find_duplicate_ids(items: List[str]) -> List[str]:
    """Find duplicate IDs in a list for validation purposes."""
    from collections import Counter

    counts = Counter(items)
    return [item for item, count in counts.items() if count > 1]


class BenchmarkValidator:
    """Base validator class containing pure model validation for benchmarks.

    This is a mixin class that expects the implementing class to provide:
    - stages: List of stage objects
    - software_environments: List of software environment objects
    - metric_collectors: Optional list of metric collector objects
    - software_backend: Software backend configuration
    - get_modules(): Method returning dict of modules
    - get_outputs(): Method returning dict of outputs
    """

    def validate_model_structure(self) -> None:
        """Validate pure model structure without execution context."""
        errors: List[str] = []

        # 1. Validate unique IDs
        stage_ids = [stage.id for stage in self.stages]  # type: ignore
        duplicate_stage_ids = _find_duplicate_ids(stage_ids)  # type: ignore
        if duplicate_stage_ids:
            errors.append(
                f"Found duplicate stage ids: {', '.join(duplicate_stage_ids)}"
            )

        all_modules = self.get_modules()  # type: ignore
        module_ids = list(all_modules.keys())  # type: ignore
        duplicate_module_ids = _find_duplicate_ids(module_ids)  # type: ignore
        if duplicate_module_ids:
            errors.append(
                f"Found duplicate module ids: {', '.join(duplicate_module_ids)}"
            )

        # NB: output ids may intentionally repeat across stages — that is the
        # gather shared-output-id contract (design 010 §3.1): a `gather.from`
        # collects every producer of the id. `get_outputs()` returns a dict
        # keyed by id, so cross-stage duplicates are collapsed here and this
        # check only ever sees unique keys. Do not "fix" it into a real
        # uniqueness check — it would break gather.
        all_outputs = self.get_outputs()  # type: ignore
        output_ids = list(all_outputs.keys())  # type: ignore
        duplicate_output_ids = _find_duplicate_ids(output_ids)  # type: ignore
        if duplicate_output_ids:
            errors.append(
                f"Found duplicate output ids: {', '.join(duplicate_output_ids)}"
            )

        # 2. Validate output file paths are relative
        for output in all_outputs.values():  # type: ignore
            if output.path.strip() == "":  # type: ignore
                errors.append(f"Output path for file {output.id} is empty")  # type: ignore
            if os.path.isabs(output.path):  # type: ignore
                errors.append(
                    f"Output path for file {output.id} must be relative, not absolute: {output.path}"  # type: ignore
                )

        # 3. Validate stage inputs reference valid outputs
        for stage in self.stages:  # type: ignore
            if stage.inputs:  # type: ignore
                for input_collection in stage.inputs:  # type: ignore
                    for input_id in input_collection.entries:  # type: ignore
                        if input_id not in output_ids:
                            errors.append(
                                f"Input with id '{input_id}' in stage '{stage.id}' is not valid"  # type: ignore
                            )

        # 4. Validate that software environment references exist
        self._validate_software_environments(errors)

        # NOTE: the "diamond" input-join check is api-gated in detect_diamond_input_joins()
        # and called (only for api < 0.7.0) from validate_execution_context(). From
        # api 0.7.0 a diamond is resolved as a fan-in join (`_select_input_bundles`
        # in cli/run.py, design 010 §3.9); older benchmarks are rejected up front.

        # Raise error if any validation failed
        if errors:
            raise ValidationError(errors)

    def detect_diamond_input_joins(self) -> List[str]:
        """Return an error per stage that joins two divergent branches (a "diamond").

        Fan-in — a stage collecting inputs from two stages on divergent branches
        (neither upstream of the other) — is gated on api ≥ 0.7.0 (design 010
        §3.9, issue #289), resolved by `_select_input_bundles`. Benchmarks below
        0.7.0 keep the old semantics, where the resolver linearises each stage
        onto a single upstream lineage and one branch falls out, surfacing only as
        an opaque "Could not resolve input <id>" at run time. This detects that
        case so the gate (validate_execution_context, api < 0.7.0) can reject it
        up front with an actionable message instead.
        """
        errors: List[str] = []

        output_producer: Dict[str, str] = {}
        for stage in self.stages:  # type: ignore
            for output in stage.outputs:  # type: ignore
                output_producer.setdefault(output.id, stage.id)  # type: ignore

        # Direct input-producing stages for each stage (deduplicated, order kept).
        direct_producers: Dict[str, List[str]] = {}
        for stage in self.stages:  # type: ignore
            producers: List[str] = []
            for input_collection in stage.inputs or []:  # type: ignore
                for input_id in input_collection.entries:  # type: ignore
                    producer = output_producer.get(input_id)
                    if (
                        producer is not None
                        and producer != stage.id  # type: ignore
                        and producer not in producers
                    ):
                        producers.append(producer)
            direct_producers[stage.id] = producers  # type: ignore

        def _is_upstream(ancestor: str, descendant: str) -> bool:
            """True if *ancestor* produces (transitively) an input of *descendant*."""
            seen: set = set()
            stack = list(direct_producers.get(descendant, []))
            while stack:
                current = stack.pop()
                if current == ancestor:
                    return True
                if current in seen:
                    continue
                seen.add(current)
                stack.extend(direct_producers.get(current, []))
            return False

        for stage in self.stages:  # type: ignore
            producers = direct_producers[stage.id]  # type: ignore
            reported = False
            for i in range(len(producers)):
                for j in range(i + 1, len(producers)):
                    left, right = producers[i], producers[j]
                    if not (_is_upstream(left, right) or _is_upstream(right, left)):
                        errors.append(
                            f"Stage '{stage.id}' collects inputs from stages "  # type: ignore
                            f"'{left}' and '{right}', which are on divergent branches "
                            f"(neither is upstream of the other). Joining branches in "
                            f"one stage requires api_version ≥ 0.7.0; bump the "
                            f"benchmark, or route the inputs through a shared upstream "
                            f"stage. See "
                            f"https://github.com/omnibenchmark/omnibenchmark/issues/289."
                        )
                        reported = True
                        break
                if reported:
                    break

        return errors

    def detect_unsatisfiable_excludes(self) -> List[str]:
        """Return an error per module that can never produce a node because its
        ``exclude`` rules leave no valid upstream lineage.

        ``exclude`` makes two module ids mutually incompatible within a lineage
        (see ``core._paths.is_lineage_excluded``). Combined with the data-flow
        dependencies, a module's excludes -- or those of the stages producing its
        inputs -- can eliminate every root ("initial"-stage) lineage that would
        supply its inputs. The module then expands to **zero nodes**, and every
        downstream consumer of its outputs fails at run time with the opaque
        "Could not resolve input <id>". Detect that here instead.

        Method: propagate, per module, the set of root-stage module ids from which
        a valid lineage to that module exists -- compatible under ``exclude``, with
        every declared input satisfiable by some compatible producer sharing that
        root. A module whose reachable-root set is empty can never run.

        The reachable set is an over-approximation (it does not track pairwise
        incompatibilities between two *different* producer modules on the same
        lineage), so an empty set is a **sound** signal: it always means zero
        nodes, never a false positive.
        """
        errors: List[str] = []
        stages = list(self.stages)  # type: ignore
        if not stages:
            return errors

        # exclude relation (symmetric in effect, matching is_lineage_excluded)
        excludes: Dict[str, set] = {}
        for stage in stages:
            for module in stage.modules:  # type: ignore
                declared = self.get_module_excludes(module.id)  # type: ignore
                if declared:
                    excludes[module.id] = set(declared)

        def incompatible(a: str, b: str) -> bool:
            return b in excludes.get(a, ()) or a in excludes.get(b, ())

        # output id -> producing stage id (first declarer)
        output_stage: Dict[str, str] = {}
        for stage in stages:
            for output in stage.outputs:  # type: ignore
                output_stage.setdefault(output.id, stage.id)  # type: ignore

        stage_by_id = {stage.id: stage for stage in stages}  # type: ignore

        def input_groups(stage) -> List[List[str]]:
            return [ic.entries for ic in (stage.inputs or [])]  # type: ignore

        def producer_stage_ids(stage) -> set:
            result: set = set()
            for group in input_groups(stage):
                for input_id in group:
                    producer = output_stage.get(input_id)
                    if producer is not None and producer != stage.id:  # type: ignore
                        result.add(producer)
            return result

        # Topologically order stages (producers before consumers). Kahn-style;
        # bail out on a cycle (reported elsewhere) to avoid spurious errors.
        deps = {stage.id: producer_stage_ids(stage) for stage in stages}  # type: ignore
        order: List[str] = []
        placed: set = set()
        while len(order) < len(stages):
            progressed = False
            for stage in stages:  # declaration order → deterministic
                if stage.id not in placed and deps[stage.id] <= placed:  # type: ignore
                    order.append(stage.id)  # type: ignore
                    placed.add(stage.id)  # type: ignore
                    progressed = True
            if not progressed:
                return errors

        # reachable roots per module id (over-approximation)
        reachable: Dict[str, set] = {}
        for stage_id in order:
            stage = stage_by_id[stage_id]
            groups = input_groups(stage)
            for module in stage.modules:  # type: ignore
                if not groups:
                    # initial stage: the module is itself a root
                    reachable[module.id] = {module.id}  # type: ignore
                    continue

                all_roots: set = (
                    set().union(*reachable.values()) if reachable else set()
                )
                module_roots: set = set()
                for root in all_roots:
                    if incompatible(module.id, root):  # type: ignore
                        continue
                    for group in groups:
                        if all(
                            output_stage.get(input_id) is not None
                            and any(
                                not incompatible(module.id, pm.id)  # type: ignore
                                and root in reachable.get(pm.id, ())
                                for pm in stage_by_id[output_stage[input_id]].modules  # type: ignore
                            )
                            for input_id in group
                        ):
                            module_roots.add(root)
                            break
                reachable[module.id] = module_roots  # type: ignore

                if module_roots:
                    continue

                culprit = self._first_unsatisfiable_input(
                    module.id,  # type: ignore
                    groups,
                    all_roots,
                    output_stage,
                    stage_by_id,
                    reachable,
                    incompatible,
                )

                # Suppress downstream cascades: if the blocking input's producer
                # stage cannot run *at all* (every module there is itself zero-node),
                # this failure is a consequence of that upstream module, which is
                # reported as its own primary cause. Only report primary causes —
                # where the input is genuinely produced somewhere but this module's
                # excludes keep it off every compatible lineage.
                if culprit is not None:
                    producers = stage_by_id[culprit[1]].modules  # type: ignore
                    if all(not reachable.get(pm.id) for pm in producers):  # type: ignore
                        continue

                detail = (
                    f" Input '{culprit[0]}' (produced by stage '{culprit[1]}') "
                    f"is only available on lineages excluded here."
                    if culprit
                    else ""
                )
                errors.append(
                    f"Module '{stage_id}/{module.id}' can never run: its "  # type: ignore
                    f"`exclude` rules (with those of its input-producing modules) "
                    f"leave no valid upstream lineage, so it expands to zero nodes "
                    f"and downstream consumers fail with 'Could not resolve input'."
                    f"{detail} Reconcile the `exclude` lists on this module and the "
                    f"stages producing its inputs."
                )
        return errors

    @staticmethod
    def _first_unsatisfiable_input(
        module_id, groups, all_roots, output_stage, stage_by_id, reachable, incompatible
    ):
        """Name the first declared input with no exclude-compatible producer lineage."""
        own_roots = {r for r in all_roots if not incompatible(module_id, r)}
        for group in groups:
            for input_id in group:
                producer_stage_id = output_stage.get(input_id)
                if producer_stage_id is None:
                    continue
                producers = stage_by_id[producer_stage_id].modules
                satisfiable_roots: set = set()
                for pm in producers:
                    if not incompatible(module_id, pm.id):
                        satisfiable_roots |= reachable.get(pm.id, set())
                if not (own_roots & satisfiable_roots):
                    return (input_id, producer_stage_id)
        return None

    def _validate_software_environments(self, errors: List[str]) -> None:
        """Validate software environment references without checking file paths."""
        env_ids = {env.id for env in self.software_environments}  # type: ignore

        # Check modules
        all_modules = self.get_modules()  # type: ignore
        for module in all_modules.values():  # type: ignore
            if module.software_environment not in env_ids:  # type: ignore
                errors.append(
                    f"Software environment with id '{module.software_environment}' is not declared. It should be listed as part of the stanza software_environments within the benchmarking YAML header."  # type: ignore
                )

        # Check metric collectors
        if self.metric_collectors:  # type: ignore
            for collector in self.metric_collectors:  # type: ignore
                if collector.software_environment not in env_ids:  # type: ignore
                    errors.append(
                        f"Software environment with id '{collector.software_environment}' for metric collector '{collector.id}' is not declared."  # type: ignore
                    )

                # Validate metric collector inputs
                all_outputs = self.get_outputs()  # type: ignore
                output_ids = list(all_outputs.keys())  # type: ignore
                for collector_input in collector.inputs:  # type: ignore
                    # Handle both string and IOFile inputs
                    if isinstance(collector_input, str):
                        input_id = collector_input
                    else:
                        input_id = collector_input.id  # type: ignore

                    if input_id not in output_ids:  # type: ignore
                        errors.append(
                            f"Input with id '{input_id}' for metric collector '{collector.id}' is not valid."  # type: ignore
                        )

    # Note: The old validate_structure method has been split into:
    # - validate_model_structure() for pure model validation (above)
    # - validate_execution_context() in the Benchmark class for path validation

    # Utility methods

    @staticmethod
    def is_initial(stage: Any) -> bool:
        """Check if a stage is an initial stage (has no inputs)."""
        return not stage.inputs or len(stage.inputs) == 0  # type: ignore

    @staticmethod
    def is_url(string: str) -> bool:
        """Check if the string is a valid URL using urlparse."""
        try:
            result = urlparse(string)
            return all(
                [result.scheme, result.netloc]
            )  # Valid if both scheme and netloc are present
        except ValueError:
            return False

    @staticmethod
    def is_absolute_path(string: str) -> bool:
        """Check if the string is an absolute path."""
        return Path(string).is_absolute()

    @staticmethod
    def get_environment_path(
        software_backend: "SoftwareBackendEnum",
        software: "SoftwareEnvironment",
        benchmark_dir: Path,
    ) -> Optional[str]:
        """Get the environment path based on software backend and environment configuration."""
        # Import here to avoid circular imports
        from .benchmark import SoftwareBackendEnum

        environment = None
        if (
            software_backend == SoftwareBackendEnum.apptainer
            or software_backend == SoftwareBackendEnum.docker
        ):
            environment = software.apptainer

        elif software_backend == SoftwareBackendEnum.conda:
            environment = software.conda

        elif software_backend == SoftwareBackendEnum.envmodules:
            environment = software.envmodule

        if not environment:
            return None

        if BenchmarkValidator.is_url(
            environment
        ) or BenchmarkValidator.is_absolute_path(environment):
            environment_path = environment
        elif software_backend == SoftwareBackendEnum.envmodules:
            environment_path = environment
        else:
            environment_path = os.path.join(benchmark_dir, environment)

        return environment_path
