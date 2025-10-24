"""cli commands related to benchmark/module execution and start"""

import os
import sys
from pathlib import Path

import click

from omnibenchmark.cli.utils.logging import logger
from omnibenchmark.cli.utils.validation import validate_benchmark
import pandas as pd
from jinja2 import Template


from datetime import datetime

from .debug import add_debug_option


def prepare_status(benchmark, out_dir, return_all=False):
    benchmark = "tests/data/Benchmark_001.yaml"
    out_dir = "out"
    b = validate_benchmark(benchmark, out_dir, echo=False)

    status_dict = {}
    status_dict["name"] = b.get_benchmark_name()
    status_dict["version"] = b.get_benchmark_version()
    aggr_results_files = list(b.get_metric_collector_output_paths())
    status_dict["results"] = {
        "observed_files": [f for f in aggr_results_files if Path(f).is_file()],
        "missing_files": [f for f in aggr_results_files if not Path(f).is_file()],
    }
    stages = list(b.get_stage_ids())

    from omnibenchmark.workflow.snakemake.scripts.utils import (
        generate_unique_repo_folder_name,
    )

    repositories_dir = Path(".snakemake") / "repos"
    modules_repo_timestamps = {}
    for st in stages:
        modules_repo_timestamps[st] = {}
        for node in b.get_nodes_by_stage_id(st):
            if node.module_id not in modules_repo_timestamps[st].keys():
                repo = node.get_repository()
                module_dir = repositories_dir / generate_unique_repo_folder_name(
                    repo["url"], repo["commit"]
                )
                if module_dir.is_dir():
                    timestamps = [
                        f.stat().st_mtime for f in module_dir.iterdir() if f.is_file()
                    ]
                    if len(timestamps) > 0:
                        timestamp = max(timestamps)
                    else:
                        timestamp = None
                else:
                    timestamp = None
                modules_repo_timestamps[st][node.module_id] = timestamp

    config = {
        "input": "",
        "output": "",
        "dataset": "[a-zA-Z0-9_.]*",
    }
    exec_paths = b.get_execution_paths()
    exec_path_dict = {}
    for exec_path_id in range(len(exec_paths)):
        exec_path_dict[exec_path_id] = {}
        exec_path = exec_paths[exec_path_id]

        # prepare init node
        node = exec_path[0]
        dataset = node.get_id().split("-")[1]

        tmp_config = config.copy()
        # tmp_config['input'] = os.path.commonpath(init_dirnames)
        tmp_config["dataset"] = dataset
        init_files = node.get_output_paths(tmp_config)
        init_files_exists = [Path(f).is_file() for f in init_files]

        init_timestamps = [
            Path(f).stat().st_mtime if init_files_exists[i] else None
            for i, f in enumerate(init_files)
        ]
        init_is_newer_file = init_files_exists

        init_timestamp_repo = modules_repo_timestamps[node.stage_id][node.module_id]
        if init_timestamp_repo is not None:
            init_is_newer_repo = [
                False if tst is None else tst > init_timestamp_repo
                for tst in init_timestamps
            ]
        else:
            init_is_newer_repo = [False for tst in init_timestamps]

        # dirnames of observed files for node
        init_dirnames = list(set([os.path.dirname(f) for f in init_files]))

        exec_path_dict[exec_path_id][node.stage_id] = {
            "node": node,
            "files": init_files,
            "exists": init_files_exists,
            "timestamps": init_timestamps,
            "dirnames": init_dirnames,
            "is_newer_file": init_is_newer_file,
            "is_newer_repo": init_is_newer_repo,
            "is_newer": [
                i and j for i, j in zip(init_is_newer_file, init_is_newer_repo)
            ],
        }

        # iterate over rest of nodes in exec path
        for iter in range(len(exec_path) - 1):
            node2 = exec_path[iter + 1]
            tmp_config = config.copy()
            tmp_config["input"] = os.path.commonpath(init_dirnames)
            tmp_config["dataset"] = dataset
            matched_files = node2.get_output_paths(tmp_config)
            matched_files_exists = [Path(f).is_file() for f in matched_files]

            matched_dirnames = list(set([os.path.dirname(f) for f in matched_files]))

            # check which dirnames match previous ones
            matched_timestamps = [
                Path(f).stat().st_mtime if matched_files_exists[i] else None
                for i, f in enumerate(matched_files)
            ]

            # compare timestamps with previous files from previous node
            # if any previous timestamp is None, consider all current timestamps as older
            max_init_timestamp = (
                datetime.now().timestamp() * 2
                if any([tst is None for tst in init_timestamps])
                else max([tst for tst in init_timestamps])
            )
            if all([iin for iin in init_is_newer_file if iin is not None]):
                matched_is_newer_file = [
                    None if tst is None else tst > max_init_timestamp
                    for tst in matched_timestamps
                ]
            else:
                matched_is_newer_file = [None for tst in matched_timestamps]

            matched_timestamp_repo = modules_repo_timestamps[node2.stage_id][
                node2.module_id
            ]
            if matched_timestamp_repo is not None:
                matched_is_newer_repo = [
                    False if tst is None else tst > matched_timestamp_repo and iin
                    for tst, iin in zip(matched_timestamps, init_is_newer_repo)
                ]
            else:
                matched_is_newer_repo = [False for tst in matched_timestamps]

            exec_path_dict[exec_path_id][node2.stage_id] = {
                "node": node2,
                "files": matched_files,
                "exists": matched_files_exists,
                "timestamps": matched_timestamps,
                "dirnames": matched_dirnames,
                "is_newer_file": matched_is_newer_file,
                "is_newer_repo": matched_is_newer_repo,
                "is_newer": [
                    i and j
                    for i, j in zip(matched_is_newer_file, matched_is_newer_repo)
                ],
            }
            init_dirnames = matched_dirnames
            init_timestamps = matched_timestamps
            init_is_newer_file = matched_is_newer_file
            init_timestamp_repo = matched_timestamp_repo

    filedict2 = {
        **{
            st: {
                nd: {
                    "observed_files": set(),
                    "missing_files": set(),
                    "invalid_files_file": set(),
                    "invalid_files_repo": set(),
                    "invalid_files": set(),
                }
                for nd in b.get_nodes_by_stage_id(st)
            }
            for st in stages
        }
    }
    for st in stages:
        for exec_path_id in exec_path_dict.keys():
            if st in exec_path_dict[exec_path_id].keys():
                de = exec_path_dict[exec_path_id][st]
                filedict2[st][de["node"]]["observed_files"] = filedict2[st][de["node"]][
                    "observed_files"
                ].union(
                    set(
                        [
                            f
                            for i, f in enumerate(de["files"])
                            if de["exists"][i] and de["is_newer"][i]
                        ]
                    )
                )
                filedict2[st][de["node"]]["missing_files"] = filedict2[st][de["node"]][
                    "missing_files"
                ].union(
                    set([f for i, f in enumerate(de["files"]) if not de["exists"][i]])
                )
                filedict2[st][de["node"]]["invalid_files"] = filedict2[st][de["node"]][
                    "invalid_files"
                ].union(
                    set(
                        [
                            f
                            for i, f in enumerate(de["files"])
                            if de["exists"][i] and (not de["is_newer"][i])
                        ]
                    )
                )
                filedict2[st][de["node"]]["invalid_files_file"] = filedict2[st][
                    de["node"]
                ]["invalid_files_file"].union(
                    set(
                        [
                            f
                            for i, f in enumerate(de["files"])
                            if not de["is_newer_file"][i]
                        ]
                    )
                )
                filedict2[st][de["node"]]["invalid_files_repo"] = filedict2[st][
                    de["node"]
                ]["invalid_files_repo"].union(
                    set(
                        [
                            f
                            for i, f in enumerate(de["files"])
                            if not de["is_newer_repo"][i]
                        ]
                    )
                )

    for st in stages:
        for nd in b.get_nodes_by_stage_id(st):
            filedict2[st][nd]["n_observed"] = len(filedict2[st][nd]["observed_files"])
            filedict2[st][nd]["n_missing"] = len(filedict2[st][nd]["missing_files"])
            filedict2[st][nd]["n_invalid"] = len(filedict2[st][nd]["invalid_files"])
            filedict2[st][nd]["n_invalid_file"] = len(
                filedict2[st][nd]["invalid_files_file"]
            )
            filedict2[st][nd]["n_invalid_repo"] = len(
                filedict2[st][nd]["invalid_files_repo"]
            )
            filedict2[st][nd]["n"] = (
                filedict2[st][nd]["n_observed"]
                + filedict2[st][nd]["n_missing"]
                + filedict2[st][nd]["n_invalid"]
            )

    status_dict["stages"] = {
        st: {
            "n": sum([filedict2[st][nd]["n"] for nd in filedict2[st].keys()]),
            "n_observed": sum(
                [filedict2[st][nd]["n_observed"] for nd in filedict2[st].keys()]
            ),
            "n_missing": sum(
                [filedict2[st][nd]["n_missing"] for nd in filedict2[st].keys()]
            ),
            "n_invalid": sum(
                [filedict2[st][nd]["n_invalid"] for nd in filedict2[st].keys()]
            ),
            "n_invalid_file": sum(
                [filedict2[st][nd]["n_invalid_file"] for nd in filedict2[st].keys()]
            ),
            "n_invalid_repo": sum(
                [filedict2[st][nd]["n_invalid_repo"] for nd in filedict2[st].keys()]
            ),
            "n_nodes": len(filedict2[st].keys()),
            "n_modules": len(set([nd.module_id for nd in filedict2[st].keys()])),
        }
        for st in stages
    }

    status_dict["total"] = {
        "n": sum([status_dict["stages"][st]["n"] for st in stages]),
        "n_observed": sum([status_dict["stages"][st]["n_observed"] for st in stages]),
        "n_missing": sum([status_dict["stages"][st]["n_missing"] for st in stages]),
        "n_nodes": sum([status_dict["stages"][st]["n_nodes"] for st in stages]),
        "n_modules": sum([status_dict["stages"][st]["n_modules"] for st in stages]),
    }
    if return_all:
        return status_dict, filedict2
    else:
        return status_dict


@click.group(name="status")
@click.pass_context
def status(ctx):
    """Show the status of benchmarks or benchmark modules."""
    ctx.ensure_object(dict)


@add_debug_option
@status.command(
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
@click.pass_context
def stat(
    ctx,
    benchmark,
    out_dir,
):
    """Run a benchmark as specified in the yaml."""
    ctx.ensure_object(dict)
    # extra_args = parse_extra_args(ctx.args)

    # b = validate_benchmark(benchmark, out_dir, echo=False)
    status_dict = prepare_status(benchmark, out_dir)
    result_file_str = "\n".join(
        [
            f"  {f} {f"({s})" if s=="missing" else ""}"
            for f, s in zip(
                status_dict["results"]["observed_files"]
                + status_dict["results"]["missing_files"],
                ["observed" for f in status_dict["results"]["observed_files"]]
                + ["missing" for f in status_dict["results"]["missing_files"]],
            )
        ]
    )
    stages = list(status_dict["stages"].keys())

    # ndstdf = pd.DataFrame(
    #     {
    #         "stage": [st for st in stages for nd in b.get_nodes_by_stage_id(st)],
    #         "module": [
    #             nd.get_id().split("-")[1]
    #             for st in stages
    #             for nd in b.get_nodes_by_stage_id(st)
    #         ],
    #     }
    # )
    # perfdf = pd.read_csv(f"{out_dir}/performances.tsv", sep="\t")
    # perfdf = perfdf.merge(ndstdf, on = "module")

    # perfdf_grouped = perfdf.groupby(['module', 'stage'])['s'].agg([
    #     'mean',
    #     'std',
    #     'min',
    #     'max',
    #     'median',
    #     ('q25', lambda x: x.quantile(0.25)),
    #     ('q75', lambda x: x.quantile(0.75))
    # ]).reset_index()

    max_str_len_stage = max([len(st) for st in stages])
    max_file_len = len(str(status_dict["total"]["n"]))

    return_str = f"""name: {status_dict['name']}
version: {status_dict['version']}

structure (#modules, #nodes):
{" -> ".join([f"{st:>{max_str_len_stage}} ({status_dict["stages"][st]["n_modules"]}, {status_dict["stages"][st]["n_nodes"]})" for st in stages])}

file completion: {int(status_dict['total']['n_observed'] / status_dict['total']['n'] * 100): 3d}% ({status_dict['total']['n_observed']}/{status_dict['total']['n']})
{"-" * (24 + max_file_len + 1 + max_file_len + 1)}
{"".join([f"    {st:>{max_str_len_stage}}: {int(status_dict["stages"][st]['n_observed'] / status_dict["stages"][st]['n'] * 100): 3d}% ({status_dict["stages"][st]['n_observed']:{max_file_len}d}/{status_dict["stages"][st]['n']:{max_file_len}d})\n" for st in stages])}
  incomplete reason:
    file missing:
{"".join([f"      {st:>{max_str_len_stage}}: {int(status_dict["stages"][st]['n_missing'] / status_dict["stages"][st]['n'] * 100): 3d}% ({status_dict["stages"][st]['n_missing']:{max_file_len}d}/{status_dict["stages"][st]['n']:{max_file_len}d})\n" for st in stages])}
    file invalid because of timestamps of dependent files:
{"".join([f"      {st:>{max_str_len_stage}}: {int(status_dict["stages"][st]['n_invalid_file'] / status_dict["stages"][st]['n'] * 100): 3d}% ({status_dict["stages"][st]['n_invalid_file']:{max_file_len}d}/{status_dict["stages"][st]['n']:{max_file_len}d})\n" for st in stages])}
    file invalid because of timestamps of dependent repo:
{"".join([f"      {st:>{max_str_len_stage}}: {int(status_dict["stages"][st]['n_invalid_repo'] / status_dict["stages"][st]['n'] * 100): 3d}% ({status_dict["stages"][st]['n_invalid_repo']:{max_file_len}d}/{status_dict["stages"][st]['n']:{max_file_len}d})\n" for st in stages])}


aggregated results: 
{result_file_str}

logs: ...
"""
    logger.info(return_str)


@add_debug_option
@status.command(
    name="report",
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
@click.pass_context
def report(
    ctx,
    benchmark,
    out_dir,
):
    ctx.ensure_object(dict)
    # extra_args = parse_extra_args(ctx.args)
    status_dict, filedict = prepare_status(benchmark, out_dir, return_all=True)
    stages = list(status_dict["stages"].keys())

    # filename, observed, missing, invalid_files_file, invalid_files_repo, module_id, node_id, stage
    all_files = []
    for st in stages:
        for nd in filedict[st].keys():
            all_files += [
                [
                    f,
                    True,
                    False,
                    f in filedict[st][nd]["invalid_files_file"],
                    f in filedict[st][nd]["invalid_files_repo"],
                    nd.module_id,
                    nd.get_id(),
                    st,
                ]
                for j, f in enumerate(filedict[st][nd]["observed_files"])
            ]
    alldf = pd.DataFrame(
        {
            "file": [a[0] for a in all_files],
            "observed": [a[1] for a in all_files],
            "missing": [a[2] for a in all_files],
            "dependent_file_newer": [a[3] for a in all_files],
            "dependent_repo_newer": [a[4] for a in all_files],
            "module": [a[5] for a in all_files],
            "node": [a[6] for a in all_files],
            "stage": [a[7] for a in all_files],
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

    template = Template("""
    <!DOCTYPE html>
    <html>
    <head>
        <title>{{ name }} - Status Report</title>
        <!-- DataTables CSS -->
        <link rel="stylesheet" href="https://cdn.datatables.net/1.13.6/css/jquery.dataTables.min.css">
        <link rel="stylesheet" href="https://cdn.datatables.net/buttons/2.4.1/css/buttons.dataTables.min.css">
    
        <!-- jQuery -->
        <script src="https://code.jquery.com/jquery-3.7.0.min.js"></script>
    
        <!-- DataTables JS -->
        <script src="https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"></script>
    
        <!-- Buttons extension -->
        <script src="https://cdn.datatables.net/buttons/2.4.1/js/dataTables.buttons.min.js"></script>
        <script src="https://cdnjs.cloudflare.com/ajax/libs/jszip/3.10.1/jszip.min.js"></script>
        <script src="https://cdnjs.cloudflare.com/ajax/libs/pdfmake/0.2.7/pdfmake.min.js"></script>
        <script src="https://cdnjs.cloudflare.com/ajax/libs/pdfmake/0.2.7/vfs_fonts.js"></script>
        <script src="https://cdn.datatables.net/buttons/2.4.1/js/buttons.html5.min.js"></script>
        <script src="https://cdn.datatables.net/buttons/2.4.1/js/buttons.print.min.js"></script>

    
    <style>
        body { 
            font-family: Arial, sans-serif; 
            margin: 20px; 
            background-color: #f5f5f5;
        }
        .container {
            max-width: 1400px;
            margin: 0 auto;
            background: white;
            padding: 30px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        .completion { 
            font-size: 24px; 
            font-weight: bold; 
            color: #2c5aa0; 
            margin: 20px 0;
        }
        h2 {
            color: #2c5aa0;
            margin-top: 40px;
            border-bottom: 2px solid #e0e0e0;
            padding-bottom: 10px;
        }
        .dataTables_wrapper {
            margin: 20px 0;
        }
        /* Column filter styling */
        tfoot input {
            width: 100%;
            padding: 3px;
            box-sizing: border-box;
        }
        tfoot select {
            width: 100%;
            padding: 3px;
        }
    </style>
    </head>
    <body>
    <div class="container">
        <h1>{{ name }} (v{{ version }})</h1>
        <p class="completion">Overall: {{ total_pct }}%</p>
    
        <h2>Stage Summary</h2>
        {{ table1|safe }}
    
        <h2>Detailed File Status</h2>
        {{ table2|safe }}
    </div>

    <script>
    $(document).ready(function() {
        // Stage summary table - simple with export buttons
        $('#status_table').DataTable({
            order: [[4, 'desc']],
            pageLength: 10,
            searching: false,
            paging: false,
            info: false,
            dom: 'Bfrtip',
            buttons: [
                'copy', 'csv', 'excel', 'pdf', 'print'
            ]
        });
        
        // Detailed file table - with column filtering and export
        $('#all_table thead tr').clone(true).addClass('filters').appendTo('#all_table thead');
        
        var table = $('#all_table').DataTable({
            orderCellsTop: true,
            fixedHeader: true,
            order: [[1, 'desc']],  // Sort by observed
            pageLength: 25,
            dom: 'Blfrtip',
            buttons: [
                'copy', 
                'csv', 
                'excel', 
                {
                    extend: 'pdf',
                    orientation: 'landscape',
                    pageSize: 'A4'
                },
                'print'
            ],
                columnDefs: [
                    {
                        // File column - wrap long paths
                        targets: 0,
                        render: function(data, type, row) {
                            if (type === 'display' && data.length > 50) {
                                return '<span title="' + data + '">' + data + '</span>';
                            }
                            return data;
                        },
                        createdCell: function(td, cellData, rowData, row, col) {
                            $(td).css({
                                'word-wrap': 'break-word',
                                'word-break': 'break-all',
                                'white-space': 'normal',
                                'max-width': '400px'
                            });
                        }
                    },
                    {
                        // Module and node columns - wrap if needed
                        targets: [5, 6],
                        createdCell: function(td, cellData, rowData, row, col) {
                            $(td).css({
                                'word-wrap': 'break-word',
                                'max-width': '200px'
                            });
                        }
                    }
                ],
            initComplete: function () {
                var api = this.api();
                
                // Add column filters in second header row
                api.columns().eq(0).each(function (colIdx) {
                    var cell = $('.filters th').eq(colIdx);
                    var title = $(cell).text();
                    
                    // Boolean columns (columns 1-4) - use select dropdown
                    if (colIdx === 1 || colIdx === 2 || colIdx === 3 || colIdx === 4) {
                        var select = $('<select><option value="">All</option></select>')
                            .appendTo($(cell).empty())
                            .on('change', function () {
                                var val = $.fn.dataTable.util.escapeRegex($(this).val());
                                api.column(colIdx).search(val ? '^' + val + '$' : '', true, false).draw();
                            });
                        
                        api.column(colIdx).data().unique().sort().each(function (d, j) {
                            select.append('<option value="' + d + '">' + d + '</option>');
                        });
                    } 
                    // Text columns - use text input
                    else {
                        $(cell).html('<input type="text" placeholder="Filter ' + title + '..." />');
                        
                        $('input', $('.filters th').eq(colIdx)).on('keyup change clear', function () {
                            if (api.column(colIdx).search() !== this.value) {
                                api.column(colIdx).search(this.value).draw();
                            }
                        });
                    }
                });
            }
        });
    });
</script>
    </body>
    </html>
    """)

    html = template.render(
        name=status_dict["name"],
        version=status_dict["version"],
        total_pct=round(
            status_dict["total"]["n_observed"] / status_dict["total"]["n"] * 100
        ),
        table1=sum_table_html,
        table2=all_table_html,
    )
    # Save HTML to file
    with open("status_report.html", "w") as f:
        f.write(html)


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
