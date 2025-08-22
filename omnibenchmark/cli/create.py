"""CLI commands for creating new benchmarks and modules using templates"""

import sys
from pathlib import Path
import subprocess

import click
import copier

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
            subprocess.run(
                ["git", "init"], cwd=target_path, check=True, capture_output=True
            )
            logger.info("Initialized Git repository")

        # Add all files and create initial commit
        subprocess.run(
            ["git", "add", "."], cwd=target_path, check=True, capture_output=True
        )

        commit_msg = f"Initial commit: Create {project_name} {project_type}\n\nGenerated with OmniBenchmark v{__version__}"

        subprocess.run(
            ["git", "commit", "-m", commit_msg],
            cwd=target_path,
            check=True,
            capture_output=True,
        )
        logger.info("Created initial commit")
    except subprocess.CalledProcessError as e:
        logger.warning(f"Git initialization failed: {e}")
    except FileNotFoundError:
        logger.warning("Git not found. Please install Git to initialize repository.")


def _make_file_executable(file_path: Path) -> None:
    """Make a file executable."""
    if file_path.exists():
        file_path.chmod(0o755)


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
    else:  # module
        logger.info(
            f"  {step_num}. Implement your module logic in the entrypoint script"
        )
        logger.info(f"  {step_num + 1}. Update documentation and tests")
    logger.info(
        f"  {step_num + 2}. Run 'ob validate' to check your configuration (coming soon!)"
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
            copier.run_copy(
                src_path=str(template_path),
                dst_path=str(target_path),
                defaults=True,
                data=copier_data,
                quiet=False,
                pretend=False,
                vcs_ref="HEAD",
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
                vcs_ref="HEAD",
                unsafe=True,
            )
        else:
            # Interactive mode
            copier.run_copy(
                src_path=str(template_path),
                dst_path=str(target_path),
                data=copier_data,
                quiet=False,
                pretend=False,
                vcs_ref="HEAD",
                unsafe=True,
            )

        # Read benchmark name for git commit
        benchmark_name = name or "OmniBenchmark project"
        if not name:
            benchmark_yaml = target_path / "benchmark.yaml"
            if benchmark_yaml.exists():
                try:
                    import yaml

                    with open(benchmark_yaml, "r") as f:
                        config = yaml.safe_load(f)
                        if config and "name" in config:
                            benchmark_name = config["name"]
                except Exception:
                    pass  # Fall back to default name

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
    required=False,
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
@click.option(
    "--environment",
    help="Environment specification type",
    type=click.Choice(["conda", "apptainer", "envmodules", "none"]),
    default=None,
)
@click.pass_context
def create_module(
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
    entrypoint,
    environment,
):
    """Create a new module from a template.

    This command scaffolds a new OmniBenchmark module project structure using
    copier templates. It creates the necessary configuration files (CITATION.cff,
    omnibenchmark.yaml), a sample entrypoint script, and documentation to get
    you started with your module development.
    """
    ctx.ensure_object(dict)
    debug = ctx.obj.get("DEBUG", False)

    # Determine the target path from positional argument or destination option
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

    try:
        # Get the template path
        template_path = Path(__file__).parent.parent / "templates" / "module"

        if target_destination is None:
            logger.info(
                f"Starting copier questionnaire for new module using OmniBenchmark v{__version__}..."
            )
            # Use current directory if no destination specified
            target_path = Path.cwd()
        else:
            logger.info(
                f"Creating new module at {target_destination} using OmniBenchmark v{__version__}"
            )
            target_path = Path(target_destination)

        # Check target directory
        if not _check_target_directory(target_path, no_input or non_interactive):
            sys.exit(1)

        # Run copier
        copier_data = {"omnibenchmark_version": __version__}

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
                "environment_type": environment or "conda",
            }

            # Handle environment type mapping and add environment-specific parameters
            if module_data["environment_type"] == "none":
                module_data["environment_type"] = (
                    "none (I will use the environments defined by the benchmarker)"
                )
                module_data["none_environment_note"] = (
                    "I will use the environments defined by the benchmarker"
                )
            elif module_data["environment_type"] == "conda":
                module_data["conda_environment"] = "env/conda.yaml"
            elif module_data["environment_type"] == "apptainer":
                module_data["apptainer_container"] = "env/container.sif"
            elif module_data["environment_type"] == "envmodules":
                module_data["envmodules_spec"] = "module-name/version"

            copier_data.update(module_data)
            copier.run_copy(
                src_path=str(template_path),
                dst_path=str(target_path),
                defaults=True,
                data=copier_data,
                quiet=False,
                pretend=False,
                vcs_ref="HEAD",
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
                vcs_ref="HEAD",
                unsafe=True,
            )

            # For --no-input mode, make the default entrypoint executable
            default_entrypoint = target_path / "run.sh"
            _make_file_executable(default_entrypoint)

            # Make bash src and test files executable if they exist
            _make_file_executable(target_path / "src" / "main.sh")
            _make_file_executable(target_path / "tests" / "test_main.sh")
        else:
            # Interactive mode
            copier.run_copy(
                src_path=str(template_path),
                dst_path=str(target_path),
                data=copier_data,
                quiet=False,
                pretend=False,
                vcs_ref="HEAD",
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
            _make_file_executable(target_path / "tests" / "test_main.sh")

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
