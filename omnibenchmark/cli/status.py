"""cli commands related to benchmark status"""

import sys

import click

from omnibenchmark.cli.utils.logging import logger


from jinja2 import Environment, FileSystemLoader

from omnibenchmark.benchmark.status.status import prepare_status, print_exec_path_dict


from .debug import add_debug_option


# @click.group(name="status")
# @click.pass_context
# def status(ctx):
#     """Show the status of benchmarks or benchmark modules."""
#     ctx.ensure_object(dict)


@add_debug_option
@click.command(
    name="status",
    context_settings=dict(
        ignore_unknown_options=True,
        allow_extra_args=True,
    ),
)
@click.option(
    "-b",
    "--benchmark",
    help="Path to benchmark yaml file or benchmark id.",
    required=True,
    envvar="OB_BENCHMARK",
    type=click.Path(exists=True),
)
@click.option(
    "--out-dir", type=str, default="out", help="Output folder name (local only)."
)
@click.option(
    "--missing_files",
    "show_missing_files",
    is_flag=True,
    help="Show missing files in the status report.",
    default=False,
)
@click.option(
    "--reason",
    "show_incomplete_reason",
    is_flag=True,
    help="Show reason for missing files in the status report.",
    default=False,
)
@click.option(
    "--logs",
    "show_logs",
    is_flag=True,
    help="Show logs for missing files in the status report.",
    default=False,
)
@click.option(
    "--json",
    "return_json",
    is_flag=True,
    help="Return the status report as JSON.",
    default=False,
)
@click.pass_context
def status(
    ctx,
    benchmark,
    out_dir,
    show_missing_files: bool = False,
    show_incomplete_reason: bool = False,
    show_logs: bool = False,
    return_json: bool = False,
):
    """Show the status of a benchmark."""
    ctx.ensure_object(dict)
    # extra_args = parse_extra_args(ctx.args)

    # b = validate_benchmark(benchmark, out_dir, echo=False)
    print(out_dir)
    status_dict, filedict, exec_path_dict = prepare_status(
        benchmark, out_dir, return_all=True
    )
    if return_json:
        import json

        print(json.dumps(status_dict, indent=4))
        # print(json.dumps(exec_path_dict, indent=4))
        # print(exec_path_dict)
        return
    result_file_str = "\n".join(
        [
            f"  {f} {f"({s})" if s=="missing" else ""}"
            for f, s in zip(
                status_dict["results"]["observed_output_files"]
                + status_dict["results"]["missing_output_files"],
                ["observed" for f in status_dict["results"]["observed_output_files"]]
                + ["missing" for f in status_dict["results"]["missing_output_files"]],
            )
        ]
    )
    stages = list(status_dict["stages"].keys())

    max_str_len_stage = max([len(st) for st in stages])
    max_file_len = len(str(status_dict["total"]["n"]))
    is_complete = status_dict["total"]["n_observed"] == status_dict["total"]["n"]

    env = Environment(loader=FileSystemLoader("omnibenchmark/benchmark/status"))
    template = env.get_template("cli_status.txt")
    result_str = template.render(
        name=status_dict["name"],
        version=status_dict["version"],
        benchmark_structure=" -> ".join(
            [
                f"{st:>{max_str_len_stage}} ({status_dict["stages"][st]["n_modules"]}, {status_dict["stages"][st]["n_nodes"]})"
                for st in stages
            ]
        ),
        files_total_n_observed=status_dict["total"]["n_observed"],
        files_total_n=status_dict["total"]["n"],
        max_str_len_stage=max_str_len_stage,
        max_file_len=max_file_len,
        stages=stages,
        stage_stats={st: status_dict["stages"][st] for st in stages},
        exec_paths_block=print_exec_path_dict(
            exec_path_dict,
            stages,
            threshold_n_missing=1,
            full=show_incomplete_reason,
            logs=show_logs,
        ),
        result_file_str=result_file_str,
        show_incomplete_reason=show_incomplete_reason and not is_complete,
        show_missing_files=show_missing_files and not is_complete,
    )
    logger.info(result_str)


def log_error_and_quit(logger, error):
    logger.error(error)
    sys.exit(1)


def log_result_and_quit(logger, success: bool, type: str):
    if success:
        logger.info(f"{type} run has finished successfully.")
        sys.exit(0)
    logger.error(f"{type} run has failed.")
    sys.exit(1)


def abort_if_user_does_not_confirm(msg: str, logger):
    _msg = f"Are you sure you want to {msg}?"
    if not click.confirm(_msg, abort=True):
        logger.debug("aborting")
        raise click.Abort()
