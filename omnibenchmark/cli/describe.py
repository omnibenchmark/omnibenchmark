"""cli commands related to benchmark infos and stats"""

import sys

import click
import json

import pandas as pd
from jinja2 import Environment, FileSystemLoader
from pathlib import Path

from omnibenchmark.benchmark.status.status import prepare_status, print_exec_path_dict
from omnibenchmark.benchmark import BenchmarkExecution
from omnibenchmark.dag import get_node_attributes

from omnibenchmark.cli.utils.logging import logger
from .debug import add_debug_option


@click.group(name="describe")
@click.pass_context
def describe(ctx):
    """Describe benchmarks and/or information about them."""
    ctx.ensure_object(dict)


@add_debug_option
@describe.command("snakemake")
@click.option(
    "--benchmark",
    "-b",
    required=True,
    type=click.Path(exists=True),
    help="Path to benchmark yaml file or benchmark id.",
    envvar="OB_BENCHMARK",
)
@click.pass_context
def snakemake_graph(ctx, benchmark: str):
    """Export a snakemake computational graph to dot format."""
    b = BenchmarkExecution(benchmark_yaml=Path(benchmark))
    if b is None:
        return
    dot = b.export_to_dot()
    click.echo(dot.to_string())


@add_debug_option
@describe.command("topology")
@click.option(
    "--benchmark",
    "-b",
    required=True,
    type=click.Path(exists=True),
    help="Path to benchmark yaml file or benchmark id.",
    envvar="OB_BENCHMARK",
)
@click.pass_context
def plot_topology(ctx, benchmark: str):
    """Export benchmark topology to mermaid diagram format."""
    b = BenchmarkExecution(benchmark_yaml=Path(benchmark))
    if b is None:
        return
    mermaid = b.export_to_mermaid()
    click.echo(mermaid)


@add_debug_option
@describe.command(
    name="status",
    context_settings=dict(
        ignore_unknown_options=False,
        allow_extra_args=False,
    ),
)
@click.argument(
    "benchmark",
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
@click.option(
    "--html",
    "return_html",
    is_flag=True,
    help="Return the status report as HTML.",
    default=False,
)
@click.option(
    "--html-file",
    type=str,
    default="status_report.html",
    help="Output HTML file name.",
)
@click.option(
    "--force",
    "overwrite_html_file",
    is_flag=True,
    help="If HTML file exists, overwrite it.",
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
    return_html: bool = False,
    html_file: str = "status_report.html",
    overwrite_html_file: bool = False,
):
    """Show the status of a benchmark.

    BENCHMARK: Path to benchmark YAML file.
    """
    ctx.ensure_object(dict)
    status_dict, filedict, exec_path_dict = prepare_status(
        benchmark, out_dir, return_all=True, cache_dir=Path(".snakemake") / "repos"
    )
    if return_json:
        logger.info(json.dumps(status_dict, indent=4))
        sys.exit(0)

    template_path = Path(__file__).parent.parent / "templates" / "status"
    env = Environment(loader=FileSystemLoader(template_path))
    stages = list(status_dict["stages"].keys())

    if not return_html:
        result_file_str = "\n".join(
            [
                f"  {f} {f"({s})" if s=="missing" else ""}"
                for f, s in zip(
                    status_dict["results"]["observed_output_files"]
                    + status_dict["results"]["missing_output_files"],
                    [
                        "observed"
                        for f in status_dict["results"]["observed_output_files"]
                    ]
                    + [
                        "missing"
                        for f in status_dict["results"]["missing_output_files"]
                    ],
                )
            ]
        )

        max_str_len_stage = max([len(st) for st in stages])
        max_file_len = len(str(status_dict["total"]["n"]))
        is_complete = status_dict["total"]["n_observed"] == status_dict["total"]["n"]

        template = env.get_template("cli_status.jinja")
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
        sys.exit(0)
    else:
        # Load benchmark to get the graph
        b = BenchmarkExecution(Path(benchmark))
        G = b.G

        # Extract graph data for vis.js
        node_stages = get_node_attributes(G, "stage", default="unknown")

        # Collect log files and file status from exec_path_dict and filedict
        logs_by_node = {}
        files_by_node = {}
        status_by_node = {}

        for exec_path_id in exec_path_dict.keys():
            for st in exec_path_dict[exec_path_id].stages:
                eps = exec_path_dict[exec_path_id].exec_path[st]
                node_key = (st, eps.node.get_id())
                if node_key not in logs_by_node:
                    logs = []
                    if eps.stdout_log is not None:
                        logs.append(str(Path(eps.stdout_log).absolute()))
                    if eps.stderr_log is not None:
                        logs.append(str(Path(eps.stderr_log).absolute()))
                    logs_by_node[node_key] = logs

        # Collect file status for each node
        for st in stages:
            for nd in filedict[st].keys():
                node_key = (st, nd.get_id())
                node_files = list(filedict[st][nd]["output_files"])
                files_by_node[node_key] = [str(Path(f).absolute()) for f in node_files]

                # Determine node status based on file validity
                if filedict[st][nd]["n_missing"] > 0:
                    status_by_node[node_key] = "missing"
                elif filedict[st][nd]["n_invalid"] > 0:
                    status_by_node[node_key] = "invalid"
                elif filedict[st][nd]["n_empty"] > 0:
                    status_by_node[node_key] = "empty"
                elif filedict[st][nd]["n_observed"] > 0:
                    status_by_node[node_key] = "complete"
                else:
                    status_by_node[node_key] = "unknown"

        # Define status colors (colorblind-friendly palette)
        status_color_dict = {
            "complete": "#999999",  # Grey
            "invalid": "#3200E6",  #
            "empty": "#F0E442",  # Yellow
            "missing": "#D50000",  #
            "unknown": "#222222",  # Dark Grey
        }

        # Create nodes data for vis.js
        graph_nodes = []
        for node in G.nodes:
            node_id = node.get_id()
            stage = node_stages.get(node, "unknown")
            node_key = (stage, node_id)

            # Get status and corresponding color
            status = status_by_node.get(node_key, "unknown")
            color = status_color_dict.get(status, "#cccccc")

            # Get files and logs for this node
            files = files_by_node.get(node_key, [])
            logs = logs_by_node.get(node_key, [])

            graph_nodes.append(
                {
                    "id": node_id,
                    "label": f"{node.module_id}\n{node.param_id if node.param_id != 'default' else ''}",
                    "color": color,
                    "stage": stage,
                    "module": node.module_id,
                    "status": status,
                    "files": files,
                    "logs": logs,
                    "title": f"Stage: {stage}\nModule: {node.module_id}\nParams: {node.param_id}\nStatus: {status}",
                }
            )

        # Create edges data for vis.js
        graph_edges = []
        for source, target in G.edges:
            graph_edges.append({"from": source.get_id(), "to": target.get_id()})

        # Create legend data for status
        graph_legend = [
            {"name": "Complete", "color": status_color_dict["complete"]},
            {"name": "Invalid", "color": status_color_dict["invalid"]},
            {"name": "Empty", "color": status_color_dict["empty"]},
            {"name": "Missing", "color": status_color_dict["missing"]},
            {"name": "Unknown", "color": status_color_dict["unknown"]},
        ]

        # Collect files for table (keeping existing code)
        logs_by_node_str = {}
        for node_key, logs in logs_by_node.items():
            logs_by_node_str[node_key] = ", ".join(logs) if logs else ""

        all_files = []
        for st in stages:
            for nd in filedict[st].keys():
                node_key = (st, nd.get_id())
                logs_str = logs_by_node_str.get(node_key, "")
                all_files += [
                    [
                        str(Path(f).absolute()),
                        f in filedict[st][nd]["observed_output_files"],
                        f in filedict[st][nd]["missing_output_files"],
                        f
                        in filedict[st][nd]["invalid_output_files_input_file_is_newer"],
                        f in filedict[st][nd]["invalid_output_files_repo_is_newer"],
                        nd.module_id,
                        nd.get_id(),
                        st,
                        logs_str,
                    ]
                    for j, f in enumerate(filedict[st][nd]["output_files"])
                ]
        alldf = pd.DataFrame(
            {
                "file": [a[0] for a in all_files],
                "observed": [a[1] for a in all_files],
                "missing": [a[2] for a in all_files],
                "dependent file newer": [a[3] for a in all_files],
                "dependent repo newer": [a[4] for a in all_files],
                "module": [a[5] for a in all_files],
                "node": [a[6] for a in all_files],
                "stage": [a[7] for a in all_files],
                "logs": [a[8] for a in all_files],
            }
        )

        sumdf = pd.DataFrame(
            [
                {
                    "Stage": st,
                    "Completed": status_dict["stages"][st]["n_observed"],
                    "Missing": status_dict["stages"][st]["n_missing"],
                    "Total": status_dict["stages"][st]["n"],
                    "Progress %": round(
                        status_dict["stages"][st]["n_observed"]
                        / status_dict["stages"][st]["n"]
                        * 100,
                        1,
                    ),
                }
                for st in stages
            ]
        )

        # Convert to HTML table
        all_table_html = alldf.to_html(
            index=False, classes="display", table_id="all_table", border=0
        )
        sum_table_html = sumdf.to_html(
            index=False, classes="display", table_id="status_table", border=0
        )

        # Add tfoot for column filters to the detailed table
        all_table_html = all_table_html.replace(
            "</table>",
            "<tfoot><tr>"
            + "".join([f"<th>{col}</th>" for col in alldf.columns])
            + "</tr></tfoot></table>",
        )

        template = env.get_template("cli_status.html.jinja")
        html = template.render(
            name=status_dict["name"],
            version=status_dict["version"],
            total_pct=round(
                status_dict["total"]["n_observed"] / status_dict["total"]["n"] * 100
            ),
            table1=sum_table_html,
            table2=all_table_html,
            graph_nodes=graph_nodes,
            graph_edges=graph_edges,
            graph_legend=graph_legend,
        )
        if Path(html_file).exists() and not overwrite_html_file:
            logger.error(
                f"HTML file {html_file} already exists. Use --force to overwrite."
            )
            sys.exit(1)
        # Save HTML to file
        with open(html_file, "w") as f:
            f.write(html)

        logger.info(f"HTML status report saved to {html_file}")
        sys.exit(0)
