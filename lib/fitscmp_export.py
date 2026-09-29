"""Export of numerical results: tables as CSV or Brotli-compressed Parquet, metadata as JSON."""
import csv
import json
import math
import os
from datetime import datetime

import numpy as np

from .fitscmp_tables import SECTION

EXPORT_FORMATS = ("csv", "parquet")


def to_jsonable(obj):
    """Recursively convert numpy types, tuples and non-finite floats to JSON-compatible values."""
    if isinstance(obj, dict):
        return {str(k): to_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return to_jsonable(obj.tolist())
    if isinstance(obj, (bool, np.bool_)):
        return bool(obj)
    if isinstance(obj, (int, np.integer)):
        return int(obj)
    if isinstance(obj, (float, np.floating)):
        return float(obj) if math.isfinite(obj) else None
    return obj


def write_json(path, obj):
    with open(path, "w") as f:
        json.dump(to_jsonable(obj), f, indent=2)
    return path


def write_table(base, columns, rows, fmt="csv"):
    """Write rows (lists aligned with `columns`) to base.csv or base.parquet; return the path."""
    if fmt not in EXPORT_FORMATS:
        raise ValueError(f"Unknown export format '{fmt}'")
    if fmt == "parquet":
        try:
            import pandas as pd  # optional dependencies (pandas + pyarrow)
        except ImportError as e:
            raise RuntimeError("Parquet export needs pandas and pyarrow (pip install pandas pyarrow)") from e
        path = base + ".parquet"
        pd.DataFrame(rows, columns=columns).to_parquet(path, compression="brotli", index=False)
        return path
    path = base + ".csv"
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(columns)
        w.writerows([["" if v is None or (isinstance(v, float) and not math.isfinite(v)) else v for v in r]
                     for r in rows])
    return path


def _split(v):
    """(numeric value or NaN, text) of a table cell."""
    if isinstance(v, (bool, np.bool_)):
        return math.nan, "yes" if v else "no"
    if isinstance(v, (int, float, np.integer, np.floating)):
        return float(v), ""
    return math.nan, "" if v is None else str(v)


def summary_rows(table):
    """Side-panel table as tidy rows: section, quantity, unit, column (image or B-A), value, text."""
    columns, rows, _ = table
    out, section = [], ""
    for row in rows:
        if row[0] == SECTION:
            section = row[1]
            continue
        label, values = row
        unit = label[label.rfind("[") + 1:-1] if label.endswith("]") else ""
        name = label[:label.rfind(" [")] if unit else label
        for col, v in zip(columns, values):
            value, text = _split(v)
            if math.isnan(value) and not text:
                continue  # empty cell
            if text.startswith("\u00d7"):  # flux ratio shown as xR
                value, text = float(text[1:]), "ratio B/A"
            out.append([section, name, unit, col, value, text])
    return out


SUMMARY_COLUMNS = ["section", "quantity", "unit", "column", "value", "text"]
PROFILE_COLUMNS = ["f_cycles_per_pix", "power_A", "power_B", "power_A_minus_B", "ratio_A_over_B", "frc"]


def export_result(base, table, meta, fc=None, fmt="csv"):
    """Result window export: summary table, radial profiles (comparison) and JSON metadata."""
    paths = [write_table(base + "_summary", SUMMARY_COLUMNS, summary_rows(table), fmt)]
    if fc is not None:
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = fc["prof_a"] / fc["prof_b"]
        prof = np.column_stack([fc["f"], fc["prof_a"], fc["prof_b"], fc["prof_d"], ratio, fc["frc"]])
        paths.append(write_table(base + "_profiles", PROFILE_COLUMNS, prof.tolist(), fmt))
    paths.append(write_json(base + "_meta.json", dict(meta, exported=datetime.now().isoformat(timespec="seconds"),
                                                     files=[os.path.basename(p) for p in paths])))
    return paths


def export_batch(base, columns, rows, meta, fmt="csv"):
    """Batch export: one row per image and JSON metadata."""
    paths = [write_table(base + "_batch", columns, rows, fmt)]
    paths.append(write_json(base + "_batch_meta.json", dict(meta, exported=datetime.now().isoformat(timespec="seconds"),
                                                           files=[os.path.basename(p) for p in paths])))
    return paths
