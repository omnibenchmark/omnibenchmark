"""CLI commands for creating new benchmarks and modules using templates"""

import sys
from pathlib import Path

import click
import copier
import git

from omnibenchmark import __version__
from omnibenchmark.cli.utils.logging import logger

from .debug import add_debug_option


def _check_target_directory(target_path: Path, no_input: bool) -> bool:
    """Check if target directory is suitable for creation. Returns True if OK to proceed."""
    if target_path.exists():
        # Check if directory has files other than .git
        non_git_files = [f for f in target_path.iterdir() if f.name != ".git"]
        if non_git_files:  # Directory exists and has non-git files
            logger.warning(f"Directory {target_path} already exists and is not empty.")
            if target_path.joinpath(".git").exists():
                logger.warning("Git repository already exists in this directory.")

            if no_input:
                logger.error(
                    "Cannot proceed with --no-input when target directory exists and is not empty."
                )
                return False
            else:
                if not click.confirm(
                    "Do you want to continue? This may overwrite existing files.",
                    abort=False,
                ):
                    logger.info("Aborted by user.")
                    return False
    else:
        # Ensure target directory exists
        target_path.mkdir(parents=True, exist_ok=True)

    return True


def _initialize_git_repo(
    target_path: Path, project_name: str, project_type: str
) -> None:
    """Initialize git repository and create initial commit."""
    git_dir = target_path / ".git"

    try:
        # Initialize git repository if it doesn't exist
        if not git_dir.exists():
            repo = git.Repo.init(target_path)
            logger.info("Initialized Git repository")
        else:
            repo = git.Repo(target_path)

        # Add all files and create initial commit
        repo.git.add(A=True)

        commit_msg = f"Initial commit: Create {project_name} {project_type}\n\nGenerated with OmniBenchmark v{__version__}"

        repo.index.commit(commit_msg)
        logger.info("Created initial commit")
    except git.GitCommandError as e:
        logger.warning(f"Git initialization failed: {e}")
    except Exception as e:
        logger.warning(f"Git initialization failed: {e}")


def _make_file_executable(file_path: Path) -> None:
    """Make a file executable."""
    if file_path.exists():
        file_path.chmod(0o755)


def _create_env_skeleton_files(
    target_path: Path, software_backend: str, benchmark_name: str
) -> None:
    """Create skeleton/recipe files and placeholders for non-conda environment backends."""
    envs_dir = target_path / "envs"
    envs_dir.mkdir(exist_ok=True)

    if software_backend == "apptainer":
        # Create Apptainer definition file (recipe to build SIF)
        def_file = envs_dir / "example.def"
        with open(def_file, "w") as f:
            f.write(f"""Bootstrap: docker
From: ubuntu:20.04

%labels
    Author {benchmark_name} Team
    Version v1.0

%help
    Container for {benchmark_name} benchmark

    This definition file can be used to build the Apptainer image:
        sudo apptainer build example.sif example.def

%post
    # Update system packages
    apt-get update && apt-get install -y \\
        python3 \\
        python3-pip \\
        r-base \\
        wget \\
        curl \\
        git

    # TODO: Install your benchmark-specific software here
    # Examples:
    # pip3 install numpy pandas scikit-learn
    # R -e "install.packages(c('tidyverse', 'ggplot2'))"

    # Clean up
    apt-get clean
    rm -rf /var/lib/apt/lists/*

%environment
    export LC_ALL=C

%runscript
    echo "Container for {benchmark_name} benchmark"
    exec "$@"
""")

        # Create conversion script for Docker images
        convert_script = envs_dir / "docker_to_sif.sh"
        with open(convert_script, "w") as f:
            f.write("""#!/bin/bash
# Convert Docker image to Apptainer SIF format
# Usage: ./docker_to_sif.sh <docker_image> [output_name]

DOCKER_IMAGE=${{1:-ubuntu:20.04}}
OUTPUT_NAME=${{2:-example.sif}}

echo "Converting Docker image $DOCKER_IMAGE to $OUTPUT_NAME"
apptainer build "$OUTPUT_NAME" "docker://$DOCKER_IMAGE"

# TODO: Replace with your actual Docker image
# Examples:
# apptainer build example.sif docker://continuumio/miniconda3
# apptainer build example.sif docker://rocker/r-ver:4.1.0
""")
        convert_script.chmod(0o755)

        # Create minimal placeholder SIF (just for validation)
        placeholder_sif = envs_dir / "example.sif"
        with open(placeholder_sif, "w") as f:
            f.write("""# PLACEHOLDER FILE - NOT A REAL SIF IMAGE
# This is a placeholder to prevent validation errors.
#
# To create the actual SIF file, use one of:
# 1. Build from definition: sudo apptainer build example.sif example.def
# 2. Convert from Docker: ./docker_to_sif.sh your-docker-image example.sif
# 3. Pull from registry: apptainer pull example.sif library://your-image
#
# Replace this file with your actual .sif image file.
""")

        logger.info(
            "Created Apptainer skeleton files: definition file, conversion script, and placeholder"
        )

    elif software_backend == "envmodules":
        # Create EasyBuild recipe skeleton
        eb_file = envs_dir / f"{benchmark_name}-1.0.0.eb"
        with open(eb_file, "w") as f:
            f.write(f"""easyblock = 'Bundle'

name = '{benchmark_name}'
version = '1.0.0'

homepage = 'https://github.com/your-org/{benchmark_name}'
description = "Software environment for {benchmark_name} benchmark"

toolchain = {{'name': 'foss', 'version': '2023a'}}

dependencies = [
    # TODO: Add your software dependencies here
    # Examples:
    # ('Python', '3.11.3'),
    # ('R', '4.3.0'),
    # ('numpy', '1.24.3', '-Python-%(pyver)s'),
    # ('pandas', '2.0.3', '-Python-%(pyver)s'),
]

# Allow system to figure out module hierarchy
moduleclass = 'tools'

# TODO: Add any additional configuration needed
""")

        logger.info("Created EasyBuild recipe file for Environment Modules")


def _convert_github_blob_to_raw(url: str) -> str:
    """Convert GitHub blob URL to raw URL.

    Examples:
        https://github.com/user/repo/blob/main/file.yml
        -> https://raw.githubusercontent.com/user/repo/main/file.yml

        https://github.com/user/repo/blob/branch/path/to/file.yml
        -> https://raw.githubusercontent.com/user/repo/branch/path/to/file.yml
    """
    import re

    # Pattern to match GitHub blob URLs
    pattern = r"https://github\.com/([^/]+)/([^/]+)/blob/(.+)"
    match = re.match(pattern, url)

    if match:
        user, repo, rest = match.groups()
        raw_url = f"https://raw.githubusercontent.com/{user}/{repo}/{rest}"
        logger.info(f"Converting GitHub blob URL to raw URL: {raw_url}")
        return raw_url

    return url


def _extract_stage_inputs(
    benchmark_path_or_url: str, stage_id: str, debug: bool = False
) -> list[str] | None:
    """Extract input IDs from a specific stage in the benchmark YAML.

    Args:
        benchmark_path_or_url: Path or URL to the benchmark YAML file
        stage_id: ID of the stage to extract inputs from
        debug: Whether to show debug information

    Returns:
        List of input IDs (e.g., ['data.raw', 'data.meta']) or None if error
    """
    import tempfile
    from urllib.parse import urlparse

    try:
        from omnibenchmark.model.benchmark import Benchmark

        # Convert GitHub blob URLs to raw URLs
        benchmark_path_or_url = _convert_github_blob_to_raw(benchmark_path_or_url)

        # Check if it's a URL
        parsed = urlparse(benchmark_path_or_url)
        is_url = parsed.scheme in ("http", "https")

        if is_url:
            # Fetch URL to temporary file
            import urllib.request

            logger.info(f"Fetching benchmark from URL: {benchmark_path_or_url}")

            with tempfile.NamedTemporaryFile(
                mode="w+", suffix=".yaml", delete=False
            ) as tmp_file:
                tmp_path = Path(tmp_file.name)

                try:
                    with urllib.request.urlopen(benchmark_path_or_url) as response:
                        content = response.read().decode("utf-8")
                        tmp_file.write(content)
                        tmp_file.flush()

                    # Load benchmark from temporary file
                    benchmark = Benchmark.from_yaml(tmp_path)
                finally:
                    # Clean up temporary file
                    if tmp_path.exists():
                        tmp_path.unlink()
        else:
            # Load from local path
            benchmark_path = Path(benchmark_path_or_url)
            if not benchmark_path.exists():
                click.echo(f"Benchmark file not found: {benchmark_path}", err=True)
                return None
            benchmark = Benchmark.from_yaml(benchmark_path)

        # Find the stage
        target_stage = None
        for stage in benchmark.stages:
            if stage.id == stage_id:
                target_stage = stage
                break

        if not target_stage:
            available_stages = [s.id for s in benchmark.stages]
            click.echo(
                f"Stage '{stage_id}' not found in benchmark. "
                f"Available stages: {', '.join(available_stages) if available_stages else 'none'}",
                err=True,
            )
            return None

        # Extract inputs
        if not target_stage.inputs:
            logger.info(
                f"Stage '{stage_id}' has no inputs (initial stage). "
                f"Module will only require --output_dir and --name."
            )
            return []

        # Flatten inputs from InputCollection objects
        input_ids = []
        for input_collection in target_stage.inputs:
            input_ids.extend(input_collection.entries)

        if debug:
            logger.debug(f"Extracted inputs for stage '{stage_id}': {input_ids}")

        return input_ids

    except FileNotFoundError:
        click.echo(f"Benchmark file not found: {benchmark_path}", err=True)
        return None
    except Exception as e:
        click.echo(f"Failed to extract stage inputs: {e}", err=True)
        if debug:
            raise
        return None


def _log_success_and_next_steps(target_path: Path, project_type: str) -> None:
    """Log success message and next steps."""
    logger.info(
        f"{project_type.title()} scaffolding created successfully with OmniBenchmark v{__version__}"
    )
    logger.info("Next steps:")

    # Check if we're already in the target directory
    current_dir = Path.cwd().resolve()
    target_dir = target_path.resolve()
    is_current_dir = current_dir == target_dir

    step_num = 1
    if not is_current_dir:
        logger.info(f"  {step_num}. cd {target_path}")
        step_num += 1

    if project_type == "benchmark":
        logger.info(f"  {step_num}. Review and customize the benchmark configuration")
        logger.info(f"  {step_num + 1}. Add your benchmark modules")
        logger.info(
            f"  {step_num + 2}. Run 'ob validate plan benchmark.yaml' to check your configuration"
        )
    else:  # module
        logger.info(
            f"  {step_num}. Implement your module logic in the entrypoint script"
        )
        logger.info(f"  {step_num + 1}. Update documentation and tests")
        logger.info(
            f"  {step_num + 2}. Run 'ob validate module' to check your configuration"
        )


@click.group(name="create")
@click.pass_context
def create(ctx):
    """Create new benchmarks or modules from templates."""
    ctx.ensure_object(dict)


@add_debug_option
@create.command(name="benchmark")
@click.argument(
    "path",
    required=False,
    type=click.Path(path_type=Path),
)
@click.option(
    "-d",
    "--destination",
    help="Destination directory for the new benchmark (deprecated, use positional argument)",
    type=click.Path(path_type=Path),
    default=None,
    hidden=True,
)
@click.option(
    "--no-input",
    help="Do not prompt for parameters and only use defaults",
    is_flag=True,
    default=False,
)
@click.option(
    "--non-interactive",
    help="Non-interactive mode (requires all mandatory parameters as flags)",
    is_flag=True,
    default=False,
)
@click.option(
    "--name",
    help="Name of the benchmark",
    type=str,
    default=None,
)
@click.option(
    "--author-name",
    help="Author name",
    type=str,
    default=None,
)
@click.option(
    "--author-email",
    help="Author email",
    type=str,
    default=None,
)
@click.option(
    "--license",
    help="License for the benchmark",
    type=click.Choice(["MIT", "Apache-2.0", "GPL-3.0", "BSD-3-Clause", "CC0-1.0"]),
    default=None,
)
@click.option(
    "--description",
    help="Description of the benchmark",
    type=str,
    default=None,
)
@click.pass_context
def create_benchmark(
    ctx,
    path,
    destination,
    no_input,
    non_interactive,
    name,
    author_name,
    author_email,
    license,
    description,
):
    """Create a new benchmark from a template.

    This command scaffolds a new benchmark project structure using copier templates.
    It creates all necessary configuration files, directory structure, and example
    modules to get you started with your benchmark.
    """
    ctx.ensure_object(dict)
    debug = ctx.obj.get("DEBUG", False)

    # Determine the target path - prioritize positional argument, then destination option
    target_destination = path or destination

    # Check for mandatory parameters when using parameter flags
    mandatory_flags = [name, author_name, author_email]
    has_mandatory_flags = any(flag is not None for flag in mandatory_flags)

    if has_mandatory_flags or non_interactive:
        # If any mandatory flag is provided or --non-interactive is used, all must be provided
        missing_flags = []
        if name is None:
            missing_flags.append("--name")
        if author_name is None:
            missing_flags.append("--author-name")
        if author_email is None:
            missing_flags.append("--author-email")

        if missing_flags:
            logger.error(
                f"When using --non-interactive or parameter flags, all mandatory parameters are required. Missing: {', '.join(missing_flags)}"
            )
            sys.exit(1)

        # Enable non-interactive mode when flags are provided
        non_interactive = True

    # Check if destination is required when using non-interactive mode
    if (no_input or non_interactive) and target_destination is None:
        flag_name = (
            "--no-input" if no_input and not non_interactive else "--non-interactive"
        )
        logger.error(f"Path is required when using {flag_name}")
        sys.exit(1)

    if debug:
        logger.debug(f"Path: {path}")
        logger.debug(f"Destination: {destination}")
        logger.debug(f"Target destination: {target_destination}")
        logger.debug(f"No input mode: {no_input}")
        logger.debug(f"Non-interactive mode: {non_interactive}")
        logger.debug(f"Parameter flags provided: {has_mandatory_flags}")

    try:
        # Get the template path
        template_path = Path(__file__).parent.parent / "templates" / "benchmark"

        if target_destination is None:
            logger.info(
                f"Starting copier questionnaire for new benchmark using OmniBenchmark v{__version__}..."
            )
            # Use current directory if no destination specified
            target_path = Path.cwd()
        else:
            logger.info(
                f"Creating new benchmark at {target_destination} using OmniBenchmark v{__version__}"
            )
            target_path = Path(target_destination)

        # Check target directory
        if not _check_target_directory(target_path, no_input or non_interactive):
            sys.exit(1)

        # Run copier
        copier_data = {"omnibenchmark_version": __version__}

        # If non-interactive, include provided parameters and enable defaults
        if non_interactive:
            # Use provided parameters for non-interactive mode
            copier_data.update(
                {
                    "benchmark_name": name,
                    "author_name": author_name,
                    "author_email": author_email,
                    "license": license or "MIT",
                    "description": description or "",
                }
            )

        # Determine whether to run copier with defaults (no questionnaire).
        # We want to use defaults when:
        #  - user explicitly requested non-interactive mode (flags),
        #  - user passed --no-input, or
        #  - stdin is not a TTY (e.g. running in CI or tests that pipe input).
        #
        # The last condition is necessary because the copier interactive
        # questionnaire requires a real TTY; in many automated/test environments
        # stdin may be present but is not a terminal. For those cases we fall
        # back to defaults to avoid copier raising an InteractiveSessionError.
        use_defaults = non_interactive or no_input

        if not use_defaults and not sys.stdin.isatty():
            # Documented fallback: tests/CI sometimes provide non-tty stdin, so
            # treat that the same as requesting defaults rather than failing.
            logger.warning(
                "Input is not a terminal (fd=0). Falling back to defaults for copier questionnaire."
            )
            use_defaults = True

        # Single, unified copier invocation to avoid duplicated calls.
        # Note: copier will show "Copying from template version None" for local templates
        # This is expected behavior and cannot be changed without using quiet=True
        copier.run_copy(
            src_path=str(template_path),
            dst_path=str(target_path),
            defaults=use_defaults,
            data=copier_data,
            quiet=False,
            pretend=False,
            unsafe=True,
        )

        # Read benchmark name and software backend for post-processing
        benchmark_name = name or "OmniBenchmark project"
        software_backend_used = "conda"  # default

        benchmark_yaml = target_path / "benchmark.yaml"
        if benchmark_yaml.exists():
            try:
                import yaml

                with open(benchmark_yaml, "r") as f:
                    config = yaml.safe_load(f)
                    if config:
                        if not name and "name" in config:
                            benchmark_name = config["name"]
                        if "software_backend" in config:
                            software_backend_used = config["software_backend"]
            except Exception:
                pass  # Fall back to defaults

        # Create skeleton/recipe files for non-conda backends
        if software_backend_used in ["apptainer", "envmodules"]:
            _create_env_skeleton_files(
                target_path, software_backend_used, benchmark_name
            )

        # Initialize git repository and create initial commit
        _initialize_git_repo(target_path, benchmark_name, "benchmark")

        # Log success and next steps
        _log_success_and_next_steps(target_path, "benchmark")

    except Exception as e:
        logger.error(f"Failed to create benchmark: {e}")
        if debug:
            raise
        sys.exit(1)


@add_debug_option
@create.command(name="module")
@click.argument(
    "path",
    required=True,
    type=click.Path(path_type=Path),
)
@click.option(
    "-d",
    "--destination",
    help="Destination directory for the new module (deprecated, use positional argument)",
    type=click.Path(path_type=Path),
    default=None,
    hidden=True,
)
@click.option(
    "-b",
    "--benchmark",
    help="Path or URL to benchmark YAML file (for stage-specific input parsing)",
    type=str,
    default=None,
)
@click.option(
    "--for-stage",
    help="Stage ID to generate module for (requires --benchmark)",
    type=str,
    default=None,
)
@click.option(
    "--no-input",
    help="Do not prompt for parameters and only use defaults",
    is_flag=True,
    default=False,
)
@click.option(
    "--non-interactive",
    help="Non-interactive mode (requires all mandatory parameters as flags)",
    is_flag=True,
    default=False,
)
@click.option(
    "--name",
    help="Name of the module",
    type=str,
    default=None,
)
@click.option(
    "--author-name",
    help="Author name",
    type=str,
    default=None,
)
@click.option(
    "--author-email",
    help="Author email",
    type=str,
    default=None,
)
@click.option(
    "--license",
    help="License for the module",
    type=click.Choice(
        ["MIT", "Apache-2.0", "GPL-3.0-or-later", "BSD-3-Clause", "CC0-1.0"]
    ),
    default=None,
)
@click.option(
    "--description",
    help="Description of the module",
    type=str,
    default=None,
)
@click.option(
    "--entrypoint",
    help="Main entrypoint script for the module",
    type=click.Choice(["run.R", "run.py", "run.sh"]),
    default=None,
)
@click.pass_context
def create_module(
    ctx,
    path,
    destination,
    benchmark,
    for_stage,
    no_input,
    non_interactive,
    name,
    author_name,
    author_email,
    license,
    description,
    entrypoint,
):
    """Create a new module from a template.

    This command scaffolds a new OmniBenchmark module project structure using
    copier templates. It creates the necessary configuration files (CITATION.cff,
    omnibenchmark.yaml), a sample entrypoint script, and documentation to get
    you started with your module development.

    When using --benchmark and --for-stage, the generated entrypoint will include
    CLI argument parsing for the specific inputs required by that stage.
    """
    ctx.ensure_object(dict)
    debug = ctx.obj.get("DEBUG", False)

    # Validate benchmark and for_stage options - they must be used together
    if for_stage and not benchmark:
        click.echo("--for-stage requires --benchmark to be specified", err=True)
        sys.exit(1)

    if benchmark and not for_stage:
        click.echo("--benchmark requires --for-stage to be specified", err=True)
        sys.exit(1)

    # Extract stage inputs if benchmark is provided
    stage_inputs = None
    stage_inputs_with_vars = None
    if benchmark and for_stage:
        stage_inputs = _extract_stage_inputs(benchmark, for_stage, debug)
        if stage_inputs is None:
            sys.exit(1)

        # Create a list of dicts with both original ID and sanitized variable name
        stage_inputs_with_vars = [
            {
                "id": input_id,
                "var_name": input_id.replace(".", "_"),
                "var_name_upper": input_id.replace(".", "_").upper(),
            }
            for input_id in stage_inputs
        ]

    # Path is mandatory, but destination option exists for backwards compatibility
    target_destination = path if path else destination

    # Check for mandatory parameters when using parameter flags
    mandatory_flags = [name, author_name, author_email]
    has_mandatory_flags = any(flag is not None for flag in mandatory_flags)

    if has_mandatory_flags or non_interactive:
        # If any mandatory flag is provided or --non-interactive is used, all must be provided
        missing_flags = []
        if name is None:
            missing_flags.append("--name")
        if author_name is None:
            missing_flags.append("--author-name")
        if author_email is None:
            missing_flags.append("--author-email")

        if missing_flags:
            logger.error(
                f"When using --non-interactive or parameter flags, all mandatory parameters are required. Missing: {', '.join(missing_flags)}"
            )
            sys.exit(1)

        # Enable non-interactive mode when flags are provided
        non_interactive = True

    if debug:
        logger.debug(f"Path: {path}")
        logger.debug(f"Destination: {destination}")
        logger.debug(f"Target destination: {target_destination}")
        logger.debug(f"No input mode: {no_input}")

    try:
        # Get the template path
        template_path = Path(__file__).parent.parent / "templates" / "module"

        logger.info(
            f"Creating new module at {target_destination} using OmniBenchmark v{__version__}"
        )
        target_path = Path(target_destination)

        # Check target directory
        if not _check_target_directory(target_path, no_input or non_interactive):
            sys.exit(1)

        # Run copier
        copier_data: dict[str, object] = {"omnibenchmark_version": __version__}

        # Add stage inputs if provided
        if stage_inputs_with_vars is not None:
            copier_data["stage_inputs"] = stage_inputs_with_vars

        if non_interactive:
            # Use provided parameters for non-interactive mode
            module_data = {
                "module_name": name,
                "module_title": name.replace("-", " ").replace("_", " ").title(),
                "author_name": author_name,
                "author_email": author_email,
                "license": license or "GPL-3.0-or-later",
                "description": description or "",
                "entrypoint": entrypoint or "run.sh",
            }

            copier_data.update(module_data)
            copier.run_copy(
                src_path=str(template_path),
                dst_path=str(target_path),
                defaults=True,
                data=copier_data,
                quiet=False,
                pretend=False,
                unsafe=True,
            )
        elif no_input:
            # Use defaults for non-interactive mode
            copier.run_copy(
                src_path=str(template_path),
                dst_path=str(target_path),
                defaults=True,
                data=copier_data,
                quiet=False,
                pretend=False,
                unsafe=True,
            )

            # For --no-input mode, make the default entrypoint executable
            default_entrypoint = target_path / "run.sh"
            _make_file_executable(default_entrypoint)

            # Make bash src files executable if they exist
            _make_file_executable(target_path / "src" / "main.sh")
        else:
            # Interactive mode
            copier.run_copy(
                src_path=str(template_path),
                dst_path=str(target_path),
                data=copier_data,
                quiet=False,
                pretend=False,
                unsafe=True,
            )

        # Make entrypoint and related files executable
        entrypoint_name = entrypoint or "run.sh"
        # In interactive mode, we need to check what entrypoint was actually created
        if not non_interactive:
            # Try to find the actual entrypoint file created
            possible_entrypoints = (
                list(target_path.glob("run.R"))
                + list(target_path.glob("run.py"))
                + list(target_path.glob("run.sh"))
            )
            if possible_entrypoints:
                entrypoint_name = possible_entrypoints[0].name

        # Make entrypoint executable
        entrypoint_file = target_path / entrypoint_name
        _make_file_executable(entrypoint_file)

        # Make bash files executable if they exist
        if entrypoint_name == "run.sh":
            _make_file_executable(target_path / "src" / "main.sh")

        # Read module name from the created file for commit message
        module_name_for_commit = name or "module"  # fallback
        if not name:
            omnibenchmark_yaml = target_path / "omnibenchmark.yaml"
            if omnibenchmark_yaml.exists():
                try:
                    import yaml

                    with open(omnibenchmark_yaml, "r") as f:
                        config = yaml.safe_load(f)
                        if config and "module" in config and "name" in config["module"]:
                            module_name_for_commit = config["module"]["name"]
                except Exception:
                    pass  # Fall back to default name

        # Initialize git repository and create initial commit
        _initialize_git_repo(target_path, module_name_for_commit, "module")

        # Log success and next steps
        _log_success_and_next_steps(target_path, "module")

    except Exception as e:
        logger.error(f"Failed to create module: {e}")
        if debug:
            raise
        sys.exit(1)
