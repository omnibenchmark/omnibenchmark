"""
Module resolution with entrypoint dereferencing.

This module implements the resolution phase that transforms abstract benchmark
definitions into concrete ResolvedNode entities ready for execution.

Resolution Process:
1. Clone/update repository to cache (~/.cache/omnibenchmark/git/)
2. Checkout specific commit to work directory (.snakemake/repos/)
3. Read entrypoint from omnibenchmark.yaml (or fallback to config.cfg)
4. Create ResolvedModule with concrete paths

The resolver acts as the bridge between abstract models (Benchmark, Module)
and concrete execution entities (ResolvedModule, ResolvedNode).
"""

import logging
import shutil
from pathlib import Path
from typing import Optional, Dict
import warnings

import configparser
import yaml

from omnibenchmark.git import clone_module_v2, get_or_update_cached_repo
from omnibenchmark.model import Module, SoftwareBackendEnum, SoftwareEnvironment
from omnibenchmark.model.resolved import ResolvedModule, ResolvedEnvironment


logger = logging.getLogger(__name__)


class ModuleResolver:
    """
    Resolves modules by cloning repositories and dereferencing entrypoints.

    This class coordinates the resolution process:
    - Uses two-tier git caching (cache + work layers)
    - Dereferences entrypoints from omnibenchmark.yaml or config.cfg
    - Creates ResolvedModule with concrete paths

    Usage:
        resolver = ModuleResolver(work_base_dir=Path(".snakemake/repos"))
        resolved = resolver.resolve(module, dirty=False)
    """

    def __init__(
        self,
        work_base_dir: Optional[Path] = None,
        cache_dir: Optional[Path] = None,
        output_dir: Optional[Path] = None,
        software_backend: Optional[SoftwareBackendEnum] = None,
        software_environments: Optional[Dict[str, SoftwareEnvironment]] = None,
        benchmark_dir: Optional[Path] = None,
    ):
        """
        Initialize the resolver.

        Args:
            work_base_dir: Base directory for module checkouts, relative to output dir
                          (default: .modules)
                          This will be used to store shallow clones for execution.
                          The output directory becomes self-contained and portable.
                          Example: out/.modules/abc123/
            cache_dir: Git cache directory for full clones
                      (default: ~/.cache/omnibenchmark/git)
            output_dir: Output directory for benchmark (for resolving environment paths)
            software_backend: Software backend type (conda, apptainer, etc.)
            software_environments: Dict of software environments {id: SoftwareEnvironment}
            benchmark_dir: Benchmark directory (where benchmark.yaml lives, for resolving env paths)
        """
        self.work_base_dir = work_base_dir or Path(".modules")
        self.cache_dir = cache_dir
        self.output_dir = output_dir or Path(".")
        self.software_backend = software_backend
        self.software_environments = software_environments or {}
        self.benchmark_dir = benchmark_dir or Path(".")

    def populate_cache(self, module: Module) -> Path:
        """
        Populate the git cache for a module (without checkout).

        This is useful for dry-run mode to fetch all repositories upfront.

        Args:
            module: Module to cache

        Returns:
            Path to cached repository
        """
        repo = module.repository
        logger.info(f"Populating cache for {module.id}: {repo.url}@{repo.commit}")

        cache_path = get_or_update_cached_repo(repo.url, self.cache_dir)
        logger.debug(f"Cached at: {cache_path}")

        return cache_path

    def resolve(
        self,
        module: Module,
        module_id: str,
        software_environment_id: str,
        dirty: bool = False,
    ) -> ResolvedModule:
        """
        Resolve a module to a ResolvedModule with concrete paths and entrypoint.

        Resolution steps:
        1. Clone repository to cache (or update if exists)
        2. Checkout specific commit to work directory
        3. Dereference entrypoint from omnibenchmark.yaml or config.cfg
        4. Create ResolvedModule

        Args:
            module: Module from benchmark definition
            module_id: Module identifier
            software_environment_id: Software environment reference
            dirty: If True, allow working copies without commits (development mode)
                  If False, enforce production mode with specific commits

        Returns:
            ResolvedModule with concrete paths and entrypoint

        Raises:
            RuntimeError: If resolution fails
        """
        # Resolve software environment first (prepare files/symlinks)
        logger.info(f"      Resolving environment '{software_environment_id}'...")
        resolved_env_path = self._resolve_environment(
            software_environment_id, module_id
        )
        if resolved_env_path:
            logger.info(
                f"        Environment: {resolved_env_path.backend_type.value} -> {resolved_env_path.reference}"
            )
        else:
            logger.info("        Environment: host (no environment file)")
        repo = module.repository

        # In production mode (dirty=False), we need a specific commit
        if not dirty and (not repo.commit or repo.commit.strip() == ""):
            raise RuntimeError(
                f"Module '{module_id}' has no commit specified. "
                "This is required in production mode (dirty=False). "
                "Use --dirty for development with working copies."
            )

        # Check if the ref is pinned (commit hash or tag) vs unpinned (branch)
        is_pinned = self._is_pinned_ref(repo.commit)

        if not dirty and not is_pinned:
            raise RuntimeError(
                f"Module '{module_id}' uses unpinned reference '{repo.commit}'. "
                "For reproducibility, use a commit hash (e.g., 'abc1234') or tag (e.g., 'v1.0.0'). "
                "Use --dirty to allow branch references during development."
            )

        logger.info(f"Resolving module '{module_id}': {repo.url}@{repo.commit}")

        # Clone and checkout to work directory
        # Work dir structure: work_base_dir/{repo_name}/{commit_hash}/
        # This allows multiple modules from the same repo at different commits
        from omnibenchmark.git.cache import parse_repo_url

        repo_path = parse_repo_url(repo.url)
        # Extract just the repo name (last component) for shorter paths
        # e.g., "github.com/user/repo" -> "repo"
        repo_name = Path(repo_path).name

        try:
            # First, clone/checkout to get the resolved commit
            # Use a temporary path based on the ref, then we'll know the actual commit
            temp_ref = repo.commit[:12] if repo.commit else "working"
            temp_work_dir = self.work_base_dir / repo_name / temp_ref

            actual_work_dir, resolved_commit = clone_module_v2(
                repository_url=repo.url,
                ref=repo.commit,
                work_dir=temp_work_dir,
                cache_dir=self.cache_dir,
            )

            # If the ref was a branch (dirty mode), the resolved_commit will be different
            # Move to the correct location based on resolved commit
            commit_prefix = resolved_commit[:7]
            final_work_dir = self.work_base_dir / repo_name / commit_prefix

            if actual_work_dir != final_work_dir and not final_work_dir.exists():
                # Move to the correct location
                import shutil

                if actual_work_dir.exists():
                    shutil.move(str(actual_work_dir), str(final_work_dir))
                    actual_work_dir = final_work_dir
            elif final_work_dir.exists() and actual_work_dir != final_work_dir:
                # Already exists at the correct location, remove temp
                import shutil

                if actual_work_dir.exists() and actual_work_dir != final_work_dir:
                    shutil.rmtree(actual_work_dir, ignore_errors=True)
                actual_work_dir = final_work_dir
        except Exception as e:
            logger.error(f"Failed to clone module '{module_id}': {e}")
            # Log full traceback for debugging
            import traceback

            logger.debug(f"Full traceback:\n{traceback.format_exc()}")
            raise RuntimeError(f"Failed to clone module '{module_id}': {e}")

        # Dereference entrypoint
        entrypoint_key = module.repository.entrypoint or "default"
        entrypoint = self._read_entrypoint(actual_work_dir, module_id, entrypoint_key)
        if entrypoint is None:
            raise RuntimeError(
                f"Failed to dereference entrypoint for module '{module_id}'. "
                f"Expected omnibenchmark.yaml or config.cfg in {actual_work_dir}"
            )

        # Entrypoint is relative to module_dir
        entrypoint_path = Path(entrypoint)

        # Make entrypoint executable
        entrypoint_full_path = actual_work_dir / entrypoint_path
        if entrypoint_full_path.exists():
            import stat

            entrypoint_full_path.chmod(
                entrypoint_full_path.stat().st_mode | stat.S_IEXEC
            )
            logger.debug(f"Made entrypoint executable: {entrypoint_full_path}")

        # Create ResolvedModule
        # module_dir should be relative to the Snakefile location (which is in out_dir)
        # Make path relative to work_base_dir to get just .modules/{commit}/
        if actual_work_dir.is_absolute():
            # Absolute path: make relative to cwd
            relative_module_dir = actual_work_dir.relative_to(Path.cwd())
        else:
            # Already relative: make relative to work_base_dir to strip the out_dir prefix
            # e.g., out/.repos/abc123 -> .repos/abc123
            relative_module_dir = actual_work_dir.relative_to(self.work_base_dir.parent)

        # Check entrypoint for shebang and infer interpreter
        entrypoint_full_path = actual_work_dir / entrypoint_path
        has_shebang, interpreter = self._check_shebang(entrypoint_full_path)

        resolved = ResolvedModule(
            repository_url=repo.url,
            commit=resolved_commit,
            module_dir=relative_module_dir,
            entrypoint=entrypoint_path,
            software_environment_id=software_environment_id,
            resolved_environment=resolved_env_path,
            has_shebang=has_shebang,
            interpreter=interpreter,
            dirty=dirty,
        )

        logger.debug(
            f"Resolved module '{module_id}': "
            f"dir={actual_work_dir}, entrypoint={entrypoint_path}"
        )

        return resolved

    def _is_pinned_ref(self, ref: str) -> bool:
        """
        Check if a git reference is pinned (commit hash or tag) vs unpinned (branch).

        Pinned references:
        - Commit hash: 7-40 character hex string (e.g., 'abc1234', 'abc1234567890...')
        - Semver tag: starts with 'v' followed by digits (e.g., 'v1.0.0', 'v2.3.4-beta')
        - Other tags: contains digits and dots like version numbers (e.g., '1.0.0', '2.3.4')

        Unpinned references (branches):
        - 'main', 'master', 'develop', 'feature/foo', etc.

        Args:
            ref: Git reference string

        Returns:
            True if the reference is pinned, False if it's a branch
        """
        import re

        if not ref:
            return False

        ref = ref.strip()

        # Check for commit hash (7-40 hex characters)
        if re.match(r"^[0-9a-fA-F]{7,40}$", ref):
            return True

        # Check for semver-like tags: v1.0.0, v2.3.4-beta, 1.0.0, etc.
        if re.match(r"^v?\d+\.\d+(\.\d+)?(-[\w.]+)?$", ref):
            return True

        # Everything else is considered a branch (unpinned)
        return False

    def _check_shebang(self, entrypoint_path: Path) -> tuple[bool, Optional[str]]:
        """
        Check if entrypoint has a valid shebang and infer interpreter.

        Args:
            entrypoint_path: Full path to entrypoint file

        Returns:
            Tuple of (has_shebang, interpreter)
            - has_shebang: True if file has valid shebang (#!)
            - interpreter: Inferred interpreter (python3, Rscript, bash) based on extension
        """
        has_shebang = False
        interpreter = None

        # Check for shebang
        shebang_line = None
        try:
            with open(entrypoint_path, "rb") as f:
                first_bytes = f.read(2)
                if first_bytes == b"#!":
                    # Read the full first line to parse the shebang
                    f.seek(0)
                    shebang_line = f.readline().decode("utf-8", errors="ignore").strip()
                    has_shebang = True
        except Exception as e:
            logger.warning(f"Could not read entrypoint file to check shebang: {e}")

        # Check for broken R shebangs (#!/usr/bin/env/R or #!/usr/bin/env R)
        # These should use Rscript instead of direct execution
        if has_shebang and shebang_line:
            # Check if shebang references R (not Rscript)
            if "/R" in shebang_line and "Rscript" not in shebang_line:
                logger.warning(
                    f"Detected shebang with 'R' instead of 'Rscript' in {entrypoint_path}. "
                    "Will use 'Rscript' interpreter instead of direct execution."
                )
                # Treat as if there's no shebang, force interpreter usage
                has_shebang = False

        # Infer interpreter from file extension
        entrypoint_str = str(entrypoint_path)
        if entrypoint_str.endswith(".py"):
            interpreter = "python3"
        elif entrypoint_str.endswith(".R"):
            interpreter = "Rscript"
        elif entrypoint_str.endswith(".sh"):
            interpreter = "bash"
        else:
            # Default fallback
            interpreter = "python3"

        return has_shebang, interpreter

    def _resolve_environment(
        self, environment_id: str, module_id: str
    ) -> Optional[ResolvedEnvironment]:
        """
        Resolve software environment by copying/symlinking files to output directory.

        This prepares the environment files in .out/envs/ so they can be referenced
        in the generated Snakefile.

        Strategy by backend:
        - conda: Copy yaml file to .envs/{env_id}.yaml
        - apptainer: Create symlink to local SIF in .envs/{env_id}.sif
                     (skip if URL-based image)
        - envmodules: Return module name as string
        - host: Return None

        Args:
            environment_id: Software environment identifier
            module_id: Module identifier (for logging)

        Returns:
            ResolvedEnvironment with backend type and reference,
            or None if no environment needed (host backend)
        """
        if not self.software_backend or not self.software_environments:
            logger.warning(
                f"No software backend/environments configured for module '{module_id}'. "
                "Skipping environment resolution."
            )
            return None

        # Get environment definition
        env = self.software_environments.get(environment_id)
        if not env:
            logger.error(
                f"Software environment '{environment_id}' not found for module '{module_id}'"
            )
            return None

        # Create envs directory (relative to output dir where Snakefile will be)
        envs_dir = self.output_dir / ".envs"
        envs_dir.mkdir(parents=True, exist_ok=True)

        # Handle based on backend type
        if self.software_backend == SoftwareBackendEnum.conda:
            return self._resolve_conda_environment(env, environment_id, envs_dir)
        elif self.software_backend == SoftwareBackendEnum.apptainer:
            return self._resolve_apptainer_environment(env, environment_id, envs_dir)
        elif self.software_backend == SoftwareBackendEnum.docker:
            return self._resolve_docker_environment(env, environment_id)
        elif self.software_backend == SoftwareBackendEnum.envmodules:
            return self._resolve_envmodules_environment(env, environment_id)
        elif self.software_backend == SoftwareBackendEnum.host:
            # Host backend doesn't need environment files
            return None
        else:
            logger.warning(
                f"Unsupported software backend '{self.software_backend}' for environment resolution"
            )
            return None

    def _resolve_conda_environment(
        self, env: SoftwareEnvironment, env_id: str, envs_dir: Path
    ) -> Optional[ResolvedEnvironment]:
        """
        Resolve conda environment by copying yaml file.

        Args:
            env: Software environment definition
            env_id: Environment identifier
            envs_dir: Directory for environment files (.envs/)

        Returns:
            ResolvedEnvironment with path to copied yaml file
        """
        if not env.conda:
            logger.error(f"Conda environment '{env_id}' has no conda field")
            return None

        # Source path (relative to benchmark_dir)
        source_path = self.benchmark_dir / env.conda

        if not source_path.exists():
            logger.error(
                f"Conda environment file not found: {source_path} (for env '{env_id}')"
            )
            return None

        # Destination path
        dest_filename = f"{env_id}.yaml"
        dest_path = envs_dir / dest_filename

        # Copy the file
        try:
            shutil.copy2(source_path, dest_path)
            logger.info(f"Copied conda environment: {source_path} -> {dest_path}")
        except Exception as e:
            logger.error(f"Failed to copy conda environment file: {e}")
            return None

        # Return ResolvedEnvironment with relative path
        relative_path = dest_path.relative_to(self.output_dir)
        return ResolvedEnvironment(
            backend_type=SoftwareBackendEnum.conda, reference=str(relative_path)
        )

    def _resolve_apptainer_environment(
        self, env: SoftwareEnvironment, env_id: str, envs_dir: Path
    ) -> Optional[ResolvedEnvironment]:
        """
        Resolve apptainer environment by creating symlink to local SIF.

        Strategy:
        - If apptainer field is a URL (oras://, docker://, http://, etc.): return URL as-is
        - If local SIF file: create symlink in .envs/ and return relative path

        Args:
            env: Software environment definition
            env_id: Environment identifier
            envs_dir: Directory for environment files (.envs/)

        Returns:
            ResolvedEnvironment with URL or path to SIF
        """
        if not env.apptainer:
            logger.error(f"Apptainer environment '{env_id}' has no apptainer field")
            return None

        # Check if URL-based image
        from urllib.parse import urlparse

        try:
            result = urlparse(env.apptainer)
            is_url = all([result.scheme, result.netloc])
        except ValueError:
            is_url = False

        if is_url:
            # URL-based image (oras://, docker://, http://, etc.)
            logger.info(
                f"Apptainer environment '{env_id}' uses remote image: {env.apptainer}"
            )
            return ResolvedEnvironment(
                backend_type=SoftwareBackendEnum.apptainer, reference=env.apptainer
            )
        else:
            # Local SIF file, create symlink
            source_path = self.benchmark_dir / env.apptainer
            if not source_path.exists():
                logger.error(
                    f"Apptainer SIF file not found: {source_path} (for env '{env_id}')"
                )
                return None

            # Destination path
            dest_filename = f"{env_id}.sif"
            dest_path = envs_dir / dest_filename

            # Create symlink (overwrite if exists)
            try:
                if dest_path.exists() or dest_path.is_symlink():
                    dest_path.unlink()
                dest_path.symlink_to(source_path.absolute())
                logger.info(
                    f"Created symlink for apptainer image: {source_path} -> {dest_path}"
                )
            except Exception as e:
                logger.error(f"Failed to create symlink for apptainer image: {e}")
                return None

            # Return path relative to output_dir
            relative_path = dest_path.relative_to(self.output_dir)
            return ResolvedEnvironment(
                backend_type=SoftwareBackendEnum.apptainer, reference=str(relative_path)
            )

    def _resolve_docker_environment(
        self, env: SoftwareEnvironment, env_id: str
    ) -> Optional[ResolvedEnvironment]:
        """
        Resolve docker environment.

        Docker environments are always URL-based (docker://, no local files).

        Args:
            env: Software environment definition
            env_id: Environment identifier

        Returns:
            ResolvedEnvironment with docker image URL
        """
        if not env.docker:
            logger.error(f"Docker environment '{env_id}' has no docker field")
            return None

        # Docker images are always URLs (docker://...)
        logger.info(f"Docker environment '{env_id}' uses image: {env.docker}")
        return ResolvedEnvironment(
            backend_type=SoftwareBackendEnum.docker, reference=env.docker
        )

    def _resolve_envmodules_environment(
        self, env: SoftwareEnvironment, env_id: str
    ) -> Optional[ResolvedEnvironment]:
        """
        Resolve envmodules environment.

        For envmodules, the module name is just a string that gets passed to 'module load'.

        TODO: Verify this approach is correct for envmodules backend
        - Consider if we need to copy easyconfig files to .out/envs/ for archival
        - For now, just returning the module name as the reference

        Args:
            env: Software environment definition
            env_id: Environment identifier

        Returns:
            ResolvedEnvironment with module name as reference
        """
        if not env.envmodule:
            logger.error(f"Envmodules environment '{env_id}' has no envmodule field")
            return None

        logger.debug(f"Environment modules backend for '{env_id}': {env.envmodule}")
        # TODO: Consider if easyconfig files should be copied for archival purposes
        # For now, envmodules are just strings (module names)
        return ResolvedEnvironment(
            backend_type=SoftwareBackendEnum.envmodules, reference=env.envmodule
        )

    def _read_entrypoint(
        self, module_dir: Path, module_id: str, entrypoint_key: str = "default"
    ) -> Optional[str]:
        """
        Read entrypoint from omnibenchmark.yaml or fall back to config.cfg.

        This implements the same logic as execution.py::_read_entrypoint
        but returns the entrypoint path instead of executing it.

        Args:
            module_dir: Path to module directory
            module_id: Module identifier (for error messages)
            entrypoint_key: Key to look up in the entrypoints dict (default: "default")

        Returns:
            Entrypoint path relative to module_dir, or None if not found
        """
        # Try new-style omnibenchmark.yaml first
        yaml_path = module_dir / "omnibenchmark.yaml"
        if yaml_path.exists():
            try:
                with open(yaml_path, "r") as f:
                    config = yaml.safe_load(f)

                if config and "entrypoints" in config:
                    entrypoints = config["entrypoints"]
                    if isinstance(entrypoints, dict) and entrypoint_key in entrypoints:
                        entrypoint = entrypoints[entrypoint_key]
                        logger.debug(
                            f"Found entrypoint '{entrypoint_key}' in omnibenchmark.yaml: {entrypoint}"
                        )
                        return entrypoint
                    elif isinstance(entrypoints, dict):
                        available = ", ".join(sorted(entrypoints.keys()))
                        logger.error(
                            f"Entrypoint '{entrypoint_key}' not found in module '{module_id}'. "
                            f"Available entrypoints: {available}"
                        )
                        return None
                    else:
                        logger.error(
                            f"Invalid omnibenchmark.yaml format in '{module_id}'. "
                            "Expected 'entrypoints' to be a dict."
                        )
                        return None
            except yaml.YAMLError as e:
                logger.error(
                    f"Failed to parse omnibenchmark.yaml in '{module_id}': {e}"
                )
                return None

        # Fall back to old-style config.cfg
        config_path = module_dir / "config.cfg"
        if config_path.exists():
            if entrypoint_key != "default":
                logger.error(
                    f"Module '{module_id}' uses deprecated config.cfg which only supports "
                    f"the 'default' entrypoint. Requested entrypoint: '{entrypoint_key}'. "
                    "Please migrate to omnibenchmark.yaml to use named entrypoints."
                )
                return None

            warnings.warn(
                f"Module '{module_id}' is using deprecated config.cfg. "
                "Please migrate to omnibenchmark.yaml. "
                "Support for config.cfg will be removed in a future version.",
                FutureWarning,
                stacklevel=3,
            )

            parser = configparser.ConfigParser()
            parser.read(config_path)

            if "DEFAULT" in parser and "SCRIPT" in parser["DEFAULT"]:
                entrypoint = parser["DEFAULT"]["SCRIPT"]
                logger.debug(f"Found entrypoint in config.cfg: {entrypoint}")
                return entrypoint
            else:
                logger.error(f"Invalid config.cfg format in '{module_id}'.")
                return None

        logger.error(
            f"No configuration file found in '{module_id}'. "
            "Expected omnibenchmark.yaml or config.cfg."
        )
        return None
