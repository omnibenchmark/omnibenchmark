import json
from pathlib import Path
from typing import Dict, List, Any

import pandas as pd

from omnibenchmark.cli.utils.logging import logger


# Performance metrics configuration
PERFORMANCE_METRICS = {
    "s": {"class": "time", "flip": True, "transform": "[0,1]"},
    "max_rss": {"class": "memory", "flip": True, "transform": "[0,1]"},
    "max_vms": {"class": "memory", "flip": True, "transform": "[0,1]"},
    "max_uss": {"class": "memory", "flip": True, "transform": "[0,1]"},
    "max_pss": {"class": "memory", "flip": True, "transform": "[0,1]"},
    "io_in": {"class": "io", "flip": True, "transform": "[0,1]"},
    "io_out": {"class": "io", "flip": True, "transform": "[0,1]"},
    "mean_load": {"class": "load", "flip": True, "transform": "[0,1]"},
    "cpu_time": {"class": "time", "flip": True, "transform": "[0,1]"},
}

# Metric class colors
METRIC_CLASS_COLORS = {
    "time": "forestgreen",
    "memory": "blue",
    "io": "orange",
    "load": "purple",
}


def generate_bettr_dashboard(
    performance_file: Path, output_file: Path, debug: bool = False
) -> bool:
    """Generate a bettr dashboard JSON from performance data.

    Args:
        performance_file: Path to performances.tsv file
        output_file: Path where the dashboard JSON should be saved
        debug: Whether to enable debug mode

    Returns:
        True if successful, False otherwise
    """
    try:
        df = pd.read_csv(performance_file, sep="\t")

        if df.empty:
            logger.error("Performance file is empty")
            return False

        metric_cols = extract_metric_columns(df, id_col="module")

        if not metric_cols:
            logger.error("No metric columns found in performance data")
            return False

        bettr_data = convert_to_bettr_format(
            df=df,
            id_col="module",
            metric_cols=metric_cols,
        )

        with open(output_file, "w") as f:
            json.dump(bettr_data, f, indent=2)

        logger.info(
            f"Successfully generated bettr dashboard with {len(bettr_data['data'])} methods"
        )
        return True

    except Exception as e:
        logger.error(f"Error generating bettr dashboard: {e}")
        if debug:
            import traceback

            logger.debug(traceback.format_exc())
        return False


def extract_metric_columns(df: pd.DataFrame, id_col: str) -> List[str]:
    """Extract metric columns from the DataFrame.

    Args:
        df: DataFrame with performance data
        id_col: Name of the ID column

    Returns:
        List of metric column names
    """
    # Exclude non-metric columns
    exclude_cols = {id_col, "path", "params", "module", "dataset", "lineage"}

    metric_cols = []
    for col in df.columns:
        if col not in exclude_cols:
            # Check if column contains numeric data and is a known performance metric
            if pd.api.types.is_numeric_dtype(df[col]) and col in PERFORMANCE_METRICS:
                metric_cols.append(col)

    return metric_cols


def convert_to_bettr_format(
    df: pd.DataFrame,
    id_col: str,
    metric_cols: List[str],
) -> Dict[str, Any]:
    """Convert DataFrame to bettr JSON format without pivoting.

    Args:
        df: DataFrame with performance data
        id_col: Name of the ID column (module)
        metric_cols: List of metric column names

    Returns:
        Dictionary in bettr format
    """
    # First pass: assign measurement numbers per (module, dataset) pair
    module_dataset_counts: Dict[tuple, int] = {}
    row_assignments = []  # Store (module, dataset, measurement_num) for each row

    for idx, row in df.iterrows():
        module_name = str(row[id_col])
        dataset = str(row["dataset"])
        pair_key = (module_name, dataset)

        if pair_key not in module_dataset_counts:
            module_dataset_counts[pair_key] = 0
        module_dataset_counts[pair_key] += 1

        row_assignments.append(
            {
                "idx": idx,
                "module": module_name,
                "dataset": dataset,
                "measurement_num": module_dataset_counts[pair_key],
            }
        )

    # Second pass: group by (module, measurement_num) to collapse datasets
    collapse_key_to_rows: Dict[tuple, List[int]] = {}

    for assignment in row_assignments:
        collapse_key = (assignment["module"], assignment["measurement_num"])
        if collapse_key not in collapse_key_to_rows:
            collapse_key_to_rows[collapse_key] = []
        collapse_key_to_rows[collapse_key].append(assignment["idx"])

    # Third pass: create collapsed bettr data
    data: List[Dict[str, Any]] = []
    id_info: List[Dict[str, Any]] = []

    for collapse_key, row_indices in sorted(collapse_key_to_rows.items()):
        module_name, measurement_num = collapse_key

        # Create unique identifier for this collapsed row
        unique_module_name = f"{module_name}_{measurement_num}"
        record = {id_col: unique_module_name}

        # Collect datasets for this collapsed row
        datasets_in_row = set()

        # Aggregate metrics from all rows with this collapse key
        for row_idx in row_indices:
            row = df.iloc[row_idx]
            dataset = str(row["dataset"])
            datasets_in_row.add(dataset)

            # Add each metric as metric_dataset column
            for metric in metric_cols:
                metric_dataset_col = f"{metric}_{dataset}"
                value = row[metric]
                if pd.isna(value):
                    record[metric_dataset_col] = None
                elif isinstance(value, (int, float)):
                    record[metric_dataset_col] = float(value)
                else:
                    try:
                        record[metric_dataset_col] = float(value)
                    except (ValueError, TypeError):
                        record[metric_dataset_col] = None

        data.append(record)

        # Create id_info for this collapsed row
        id_info.append(
            {
                id_col: unique_module_name,
                "Measurement": f"#{measurement_num}",
            }
        )

    # Get all unique metric_dataset combinations from the data
    all_metric_dataset_cols = set()
    for record in data:
        for key in record.keys():
            if key != id_col:
                all_metric_dataset_cols.add(key)

    # Generate metricInfo for each metric_dataset column
    metric_info = []
    for metric_dataset_col in sorted(all_metric_dataset_cols):
        # Split on last underscore to get metric and dataset
        parts = metric_dataset_col.rsplit("_", 1)
        if len(parts) == 2:
            metric, dataset = parts
        else:
            metric = metric_dataset_col
            dataset = ""

        metric_config = PERFORMANCE_METRICS.get(
            metric, {"class": "unknown", "flip": False}
        )
        metric_info.append(
            {
                "Metric": metric_dataset_col,
                "Class": metric_config["class"],
                "Dataset": dataset,
            }
        )

    # Generate initialTransforms for metric_dataset columns that need flipping
    initial_transforms = {}
    for metric_dataset_col in all_metric_dataset_cols:
        parts = metric_dataset_col.rsplit("_", 1)
        metric = parts[0] if len(parts) == 2 else metric_dataset_col
        metric_config = PERFORMANCE_METRICS[metric]
        initial_transforms[metric_dataset_col] = {
            "flip": metric_config.get("flip", False),
            "transform": metric_config.get("transform", None),
        }

    # Generate metricColors using static configuration
    metric_colors = {"Class": METRIC_CLASS_COLORS.copy()}

    bettr_json = {
        "idCol": id_col,
        "data": data,
        "idInfo": id_info,
        "metricInfo": metric_info,
        "metricColors": metric_colors,
        "initialWeights": {},
        "initialTransforms": initial_transforms,
        "idColors": {},
    }

    return bettr_json


# def simulate_performance(orig: pd.DataFrame) -> pd.DataFrame:
#     import pandas as pd
#     import numpy as np
#
#     n = len(orig)
#
#     sim = pd.DataFrame({
#         's': np.random.uniform(1.0, 10.0, n),
#         'max_rss': np.random.uniform(1e6, 5e7, n),
#         'max_vms': np.random.uniform(5e7, 2e8, n),
#         'max_uss': np.random.uniform(5e5, 2e7, n),
#         'max_pss': np.random.uniform(5e5, 2e7, n),
#         'io_in': np.random.uniform(0, 1e6, n),
#         'io_out': np.random.uniform(0, 1e6, n),
#         'cpu_time': np.random.uniform(0.1, 5.0, n),
#         'mean_load': np.random.uniform(10.0, 90.0, n),
#         'module': orig['module'],
#         'dataset': orig['dataset'],
#         'lineage': orig['lineage'],
#         'path': orig['path'],
#         'params': orig['params']
#     })
#
#     return sim
#
#
# if __name__ == "__main__":
#     # performance_tsv = Path("/Users/dani/Documents/omni/omni-py/out/performances.tsv")
#     # performance_tsv = pd.read_csv(performance_tsv, sep="\t")
#     # simulated_tsv = simulate_performance(performance_tsv)
#     # simulated_tsv.to_csv(Path("/Users/dani/Documents/omni/omni-py/out/sim_performances2.tsv"), sep="\t", index=False)
#     performance_tsv = Path("/Users/dani/Documents/omni/omni-py/out/sim_performances2.tsv")
#     bettr_json = Path("/Users/dani/Documents/omni/omni-py/out/bettr_performances.json")
#     generate_bettr_dashboard(performance_tsv, bettr_json)
