"""cli commands related to benchmark/module execution and start"""

import sys

import click

from omnibenchmark.cli.utils.logging import logger


import pandas as pd
from jinja2 import Template
from jinja2 import Environment, FileSystemLoader

from omnibenchmark.benchmark.status.status import prepare_status, print_exec_path_dict


from .debug import add_debug_option


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
@click.pass_context
def statuss(
    ctx,
    benchmark,
    out_dir,
    report_html: bool = False,
    show_missing_files: bool = False,
    show_incomplete_reason: bool = False,
):
    """Run a benchmark as specified in the yaml."""
    ctx.ensure_object(dict)
    # extra_args = parse_extra_args(ctx.args)

    # b = validate_benchmark(benchmark, out_dir, echo=False)
    print(out_dir)
    status_dict, filedict, exec_path_dict = prepare_status(
        benchmark, out_dir, return_all=True
    )
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

    max_str_len_stage = max([len(st) for st in stages])
    max_file_len = len(str(status_dict["total"]["n"]))
    ascii_stages_header = "\n".join(
        [f"{'|' * (i + 1)}{st}" for i, st in enumerate(stages)]
    )
    ascii_stages_sep = "|" * len(stages)
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
        ascii_stages_header=ascii_stages_header,
        ascii_stages_sep=ascii_stages_sep,
        exec_paths_block=print_exec_path_dict(
            exec_path_dict, stages, threshold_n_missing=1, full=show_incomplete_reason
        ),
        result_file_str=result_file_str,
        show_incomplete_reason=show_incomplete_reason and not is_complete,
        show_missing_files=show_missing_files and not is_complete,
    )
    logger.info(result_str)


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
