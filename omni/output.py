"""Methods to automatically expand output file definitions"""

from pathlib import Path
from typing import List, Mapping, Optional, Union

import omni_schema.datamodel.omni_schema as model

from .utils import as_list


class IOCollection:
    """Collection of explicit Input-Output-Parameter mappings associated with a specific benchmark"""

    def __init__(
        self,
        benchmark: model.Benchmark,
        io_mappings: Union[List[Optional[IOMapping]], IOMapping] = [],
    ):
        self.id = benchmark.id
        self.name = benchmark.name
        self.version = benchmark.version
        self.io = {io.id: io for io in as_list(io_mappings)}

    def add_io_mapping(self, io_mapping=IOMapping):
        self.io[io_mapping.id] = io_mapping
        return self

    def get_outputs(self, io_id: str) -> List[ExplicitIO]:
        io_mapping = self.io.get(io_id)
        if not io_mapping:
            raise ValueError("IOMapping with id '{io_id}' does not exist.")
        ## FIX: Replace by output pattern
        return io_mapping.outputs


class IOMapping:
    """Explicit Input-Output-Parameter mapping that corresponds to an IOFile with implicit specification"""

    def __init__(
        self,
        io_file: model.IOFile,
        input: Optional[model.InputCollection],
        module: model.Module,
        stage: str,
        io_collection: IOCollection,
    ):
        self.id = io_file.id
        self.name = io_file.name
        self.template = io_file.path
        self.module = module
        self.stage = stage
        self.implicit_inputs = input.entries if input else []
        self.explicit_inputs = get_explicit_inputs(self.implicit_inputs, io_collection)


class ExplicitIO:
    """Explicit IO paths and their lineages"""

    def __init__(
        self,
        explicit_path: str,
        lineage: List[str],
        module: model.Module,
    ):
        self.explicit_path = explicit_path
        # Fix adjust to module name + params + all inputs? (do we need a group/hierarchy concept, e.g. module, group)
        self.lineage = lineage.append(module.id)


def get_explicit_inputs(
    implicit_inputs: List[Optional[str]], io_collection: IOCollection
) -> Mapping:
    missing_id = [
        input_id
        for input_id in implicit_inputs
        if input_id not in io_collection.io.keys()
    ]
    if len(missing_id) > 0:
        raise ValueError(
            f"Could not resolve {missing_id} in the provided IOCollection\n"
            f'Make sure to account for stage dependencies, e.g. "get_stage_dependencies()"'
        )

    unmatched = {
        input_id: io_collection.get_outputs(input_id) for input_id in implicit_inputs
    }
    matched = []
    ref_input = max(unmatched, key=lambda k: len(unmatched[k]))
    for output in unmatched.get(ref_input):
        IO_group = match_IO_groups(ref_input, unmatched)
        matched.append(IO_group)
    return matched


def match_IO_groups(
    ref: ExplicitIO, rel: Mapping[str : List[ExplicitIO]]
) -> Mapping[str, ExplicitIO]:
    pass


def get_stage_dependencies(benchmark: model.Benchmark) -> List:
    pass


def expand_benchmark(benchmark: model.Benchmark) -> IOCollection:
    pass
