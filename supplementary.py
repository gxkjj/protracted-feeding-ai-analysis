"""Validate and reproduce the released S1 and S5 mouse-level endpoints."""

import csv
import math
from pathlib import Path


def statistics(data_dir, mean, sem, welch_t_test, paired_t_test):
    folder = data_dir / "supplementary"
    with (folder / "figure_map.csv").open(newline="", encoding="utf-8") as handle:
        specs = list(csv.DictReader(handle))
    expected = {(panel, col) for panel in ("S1A", "S1B", "S1C", "S1D", "S5A", "S5B", "S5C") for col in range(1, 5)}
    actual = [(s["figure"], int(s["column"])) for s in specs]
    if len(actual) != len(expected) or set(actual) != expected:
        raise ValueError("Supplementary panel map must cover S1 and S5 exactly once")
    files = [s["source_file"] for s in specs]
    actual_files = {p.relative_to(folder).as_posix() for p in folder.rglob("*.csv")}
    if len(set(files)) != 28 or actual_files != set(files) | {"figure_map.csv", "mouse_values.csv"}:
        raise ValueError("Unexpected supplementary source-data set")
    with (folder / "mouse_values.csv").open(newline="", encoding="utf-8") as handle:
        mouse_rows = list(csv.DictReader(handle))
    metric_keys = {"Mean IBI": "Mean IBI (s)", "Bout Number": "Bout Number", "Latency": "Latency (s)", "Total feeding time": "Total feeding time (s)"}
    results = []
    for spec in specs:
        relative = Path(spec["source_file"])
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError("Source paths must stay inside the data folder")
        names = [spec["group_1"], spec["group_2"]]
        paired = spec["figure"] == "S5C"
        test = "Paired t test" if paired else "Welch t test"
        if spec["test"] != test:
            raise ValueError("Unexpected supplementary statistical test")
        values = [[], []]
        with (folder / relative).open(newline="", encoding="utf-8") as handle:
            reader = csv.reader(handle)
            if next(reader) != names:
                raise ValueError(f"Group names do not match: {relative}")
            blank_seen = [False, False]
            for line, row in enumerate(reader, 2):
                if len(row) != 2 or not any(row) or (paired and not all(row)):
                    raise ValueError(f"Missing or malformed observation: {relative}:{line}")
                for i, raw in enumerate(row):
                    if not raw:
                        blank_seen[i] = True
                        continue
                    value = float(raw)
                    if blank_seen[i] or not math.isfinite(value) or value < 0:
                        raise ValueError(f"Invalid or re-aligned observation: {relative}:{line}")
                    if spec["metric"] == "Bout Number" and not value.is_integer():
                        raise ValueError("Bout Number must be an integer")
                    values[i].append(value)
        pair_ids = []
        for i, name in enumerate(names):
            if len(values[i]) != int(spec[f"n_{i + 1}"]):
                raise ValueError(f"Sample count mismatch: {relative}")
            records = sorted((r for r in mouse_rows if r["figure"] == spec["figure"] and r["group"] == name), key=lambda r: int(r["row_in_group"]))
            if [int(r["row_in_group"]) for r in records] != list(range(1, len(values[i]) + 1)):
                raise ValueError(f"Mouse row map mismatch: {relative}")
            if [float(r[metric_keys[spec["metric"]]]) for r in records] != values[i]:
                raise ValueError(f"Mouse values do not match panel table: {relative}")
            pair_ids.append([r["pair_id"] for r in records])
        if paired and (pair_ids[0] != pair_ids[1] or not all(pair_ids[0]) or len(set(pair_ids[0])) != len(pair_ids[0])):
            raise ValueError("DAT-Cre pair IDs must match uniquely in row order")
        statistic, df, p_value = (paired_t_test if paired else welch_t_test)(*values)
        results.append({
            "figure": spec["figure"], "column": spec["column"], "metric": spec["metric"], "unit": spec["unit"],
            "source_file": spec["source_file"], "test": test,
            "group_1": names[0], "n_1": len(values[0]), "mean_1": mean(values[0]), "sem_1": sem(values[0]),
            "group_2": names[1], "n_2": len(values[1]), "mean_2": mean(values[1]), "sem_2": sem(values[1]),
            "statistic": statistic, "df": df, "p_value": p_value,
        })
    return results
