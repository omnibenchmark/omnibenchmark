from typing import TYPE_CHECKING, Dict, List, Any

if TYPE_CHECKING:
    import pandas as pd

try:
    import pandas as pd

    EXPORT_AVAILABLE = True
except ImportError:
    EXPORT_AVAILABLE = False


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


def create_bettr_dashboard(
    df: "pd.DataFrame", id_col: str = "module"
) -> Dict[str, Any]:
    """Create bettr dashboard dict from performance DataFrame.

    Parses metrics into structured format, then serializes to bettr JSON.
    """
    if not EXPORT_AVAILABLE:
        raise ImportError(
            "Export functionality is not available. Install with: pip install omnibenchmark[export]"
        )

    if df.empty:
        raise ValueError("performance.tsv is empty")

    metric_cols = _extract_metric_columns(df, id_col=id_col)
    if not metric_cols:
        raise ValueError("No metric columns found in performance data")

    parsed_data = _parse_metrics(df, metric_cols, id_col=id_col)
    bettr_data = _serialize_to_bettr_format(parsed_data, id_col=id_col)

    return bettr_data


def _extract_metric_columns(df: "pd.DataFrame", id_col: str) -> List[str]:
    """Extract numeric columns that are defined in PERFORMANCE_METRICS.

    Excludes metadata columns (path, params, module, dataset, lineage).
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


def _parse_metrics(
    df: "pd.DataFrame",
    metric_cols: List[str],
    id_col: str,
) -> List[Dict]:
    """Parse metrics into structured records with explicit metric/dataset/value dicts.

    Groups rows by (module, measurement_num) and stores each metric as
    {"metric": str, "dataset": str, "value": float}.
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

    # Third pass: create structured data records
    records: List[Dict[str, Any]] = []

    for collapse_key, row_indices in sorted(collapse_key_to_rows.items()):
        module_name, measurement_num = collapse_key

        # Create unique identifier for this collapsed row
        module_id = f"{module_name}_{measurement_num}"
        metrics = []

        # Collect metrics from all rows with this collapse key
        for row_idx in row_indices:
            row = df.iloc[row_idx]
            dataset = str(row["dataset"])

            # Add each metric with explicit structure
            for metric in metric_cols:
                value = row[metric]
                if pd.isna(value):
                    parsed_value = None
                elif isinstance(value, (int, float)):
                    parsed_value = float(value)
                else:
                    try:
                        parsed_value = float(value)
                    except (ValueError, TypeError):
                        parsed_value = None

                metrics.append(
                    {
                        "metric": metric,
                        "dataset": dataset,
                        "value": parsed_value,
                    }
                )

        records.append(
            {
                "module_id": module_id,
                "module": module_name,
                "measurement_num": measurement_num,
                "metrics": metrics,
            }
        )

    return records


def _serialize_to_bettr_format(
    records: List[Dict],
    id_col: str,
) -> Dict[str, Any]:
    """Serialize structured metrics to bettr JSON format.

    Flattens metrics into metric_dataset columns and generates bettr metadata
    """

    # Flatten structured metrics into bettr data format
    data = []
    id_info = []
    seen_metric_dataset_combinations = set()

    for record in records:
        module_id = record["module_id"]
        measurement_num = record["measurement_num"]

        # Create flattened data row
        data_row = {id_col: module_id}

        # Flatten metrics into metric_dataset columns
        for metric_entry in record["metrics"]:
            metric = metric_entry["metric"]
            dataset = metric_entry["dataset"]
            value = metric_entry["value"]

            metric_dataset_col = f"{metric}_{dataset}"
            data_row[metric_dataset_col] = value
            seen_metric_dataset_combinations.add((metric, dataset))

        data.append(data_row)

        # Create id_info entry
        id_info.append(
            {
                id_col: module_id,
                "Measurement": f"#{measurement_num}",
            }
        )

    # Generate metricInfo from seen combinations
    metric_info = []
    for metric, dataset in sorted(seen_metric_dataset_combinations):
        metric_dataset_col = f"{metric}_{dataset}"
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

    # Generate initialTransforms
    initial_transforms = {}
    for metric, dataset in seen_metric_dataset_combinations:
        metric_dataset_col = f"{metric}_{dataset}"
        metric_config = PERFORMANCE_METRICS[metric]
        initial_transforms[metric_dataset_col] = {
            "flip": metric_config.get("flip", False),
            "transform": metric_config.get("transform", None),
        }

    # Generate metricColors using static configuration
    metric_colors = {"Class": METRIC_CLASS_COLORS.copy()}

    result = {
        "idCol": id_col,
        "data": data,
        "idInfo": id_info,
        "metricInfo": metric_info,
        "metricColors": metric_colors,
        "initialWeights": {},
        "initialTransforms": initial_transforms,
        "idColors": {},
    }

    return result
