#!/usr/bin/env python3
"""Reproduce the manuscript statistics from the released source data."""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Sequence


SOURCE_FILES = (
    "DAT-Cre.csv",
    "HFD female.csv",
    "HFD male female.csv",
    "HFD male.csv",
    "HFD rescue.csv",
    "Int-HFD male female manual.csv",
    "Restraint Mosapride.csv",
    "Restraint male female manual.csv",
    "isolation female.csv",
    "isolation male female manual.csv",
    "isolation male female.csv",
    "isolation male.csv",
    "isolation rescue.csv",
    "manual vs ai all exclude model.csv",
    "restraint female.csv",
    "restraint male female.csv",
    "restraint male.csv",
)

@dataclass(frozen=True)
class ComparisonSpec:
    figure: str
    source_file: str
    test: str
    group_1: str
    group_2: str
    reported_p: str


COMPARISONS = (
    ComparisonSpec("Figure 2A", "HFD male female.csv", "Welch t test", "Control", "Int-HFD", "0.0110"),
    ComparisonSpec("Figure 2B", "isolation male female.csv", "Welch t test", "Control", "Isolation", "0.0213"),
    ComparisonSpec("Figure 2C", "Int-HFD male female manual.csv", "Welch t test", "Control", "Int-HFD", "<0.0001"),
    ComparisonSpec("Figure 2D", "isolation male female manual.csv", "Welch t test", "Control", "Isolation", "0.0096"),
    ComparisonSpec("Figure 3B", "restraint male female.csv", "Welch t test", "Control", "Restraint", "0.1661"),
    ComparisonSpec("Figure 3C", "Restraint male female manual.csv", "Welch t test", "Control", "Restraint", "0.4492"),
    ComparisonSpec("Figure 3E", "Restraint Mosapride.csv", "Welch t test", "Restraint +Control", "Restraint +Mosapride", "0.0374"),
    ComparisonSpec("Figure 4E", "HFD rescue.csv", "Welch t test", "Int-HFD", "Int-HFD+DA", "0.0332"),
    ComparisonSpec("Figure 4F", "isolation rescue.csv", "Welch t test", "Isolation", "Isolation+DA", "0.0440"),
    ComparisonSpec("Figure 4H", "DAT-Cre.csv", "Paired t test", "DAT-Cre Control", "DAT-Cre CNO", "0.0360"),
    ComparisonSpec("Figure S1A (male)", "HFD male.csv", "Welch t test", "Control", "Int-HFD", "0.0720"),
    ComparisonSpec("Figure S1A (female)", "HFD female.csv", "Welch t test", "Control", "Int-HFD", "0.0380"),
    ComparisonSpec("Figure S1B (male)", "isolation male.csv", "Welch t test", "Control", "Isolation", "0.0470"),
    ComparisonSpec("Figure S1B (female)", "isolation female.csv", "Welch t test", "Control", "Isolation", "0.0410"),
    ComparisonSpec("Figure S1C (male)", "restraint male.csv", "Welch t test", "Control", "Restraint", "0.0684"),
    ComparisonSpec("Figure S1C (female)", "restraint female.csv", "Welch t test", "Control", "Restraint", "0.7205"),
)


def normalize_header(value: str) -> str:
    return " ".join(value.replace("\r", " ").replace("\n", " ").split())


def read_columns(path: Path) -> dict[str, tuple[float, ...]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = csv.reader(handle)
        try:
            headers = [normalize_header(value) for value in next(rows)]
        except StopIteration:
            raise ValueError(f"Empty CSV: {path.name}") from None
        columns: dict[str, list[float]] = {header: [] for header in headers}
        for row_number, row in enumerate(rows, start=2):
            padded = row + [""] * (len(headers) - len(row))
            for header, raw_value in zip(headers, padded):
                value = raw_value.strip()
                if not value:
                    continue
                try:
                    columns[header].append(float(value))
                except ValueError as error:
                    raise ValueError(
                        f"Non-numeric value in {path.name}, row {row_number}, column {header}: {value!r}"
                    ) from error
    return {name: tuple(values) for name, values in columns.items()}


def mean(values: Sequence[float]) -> float:
    return math.fsum(values) / len(values)


def sample_variance(values: Sequence[float]) -> float:
    center = mean(values)
    return math.fsum((value - center) ** 2 for value in values) / (len(values) - 1)


def sem(values: Sequence[float]) -> float:
    return math.sqrt(sample_variance(values) / len(values))


def _beta_continued_fraction(a: float, b: float, x: float) -> float:
    maximum_iterations = 300
    epsilon = 3.0e-14
    tiny = 1.0e-300
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < tiny:
        d = tiny
    d = 1.0 / d
    result = d
    for iteration in range(1, maximum_iterations + 1):
        twice = 2 * iteration
        aa = iteration * (b - iteration) * x / ((qam + twice) * (a + twice))
        d = 1.0 + aa * d
        c = 1.0 + aa / c
        if abs(d) < tiny:
            d = tiny
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        result *= d * c
        aa = -(a + iteration) * (qab + iteration) * x / ((a + twice) * (qap + twice))
        d = 1.0 + aa * d
        c = 1.0 + aa / c
        if abs(d) < tiny:
            d = tiny
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        result *= delta
        if abs(delta - 1.0) < epsilon:
            return result
    raise ArithmeticError("Incomplete beta calculation did not converge")


def _regularized_incomplete_beta(a: float, b: float, x: float) -> float:
    if x in (0.0, 1.0):
        return x
    factor = math.exp(
        math.lgamma(a + b)
        - math.lgamma(a)
        - math.lgamma(b)
        + a * math.log(x)
        + b * math.log1p(-x)
    )
    if x < (a + 1.0) / (a + b + 2.0):
        return factor * _beta_continued_fraction(a, b, x) / a
    return 1.0 - factor * _beta_continued_fraction(b, a, 1.0 - x) / b


def two_sided_t_p(t_statistic: float, degrees_of_freedom: float) -> float:
    x = degrees_of_freedom / (degrees_of_freedom + t_statistic**2)
    return _regularized_incomplete_beta(degrees_of_freedom / 2.0, 0.5, x)


def welch_t_test(first: Sequence[float], second: Sequence[float]) -> tuple[float, float, float]:
    variance_1 = sample_variance(first)
    variance_2 = sample_variance(second)
    term_1 = variance_1 / len(first)
    term_2 = variance_2 / len(second)
    t_statistic = (mean(first) - mean(second)) / math.sqrt(term_1 + term_2)
    degrees_of_freedom = (term_1 + term_2) ** 2 / (
        term_1**2 / (len(first) - 1) + term_2**2 / (len(second) - 1)
    )
    return t_statistic, degrees_of_freedom, two_sided_t_p(t_statistic, degrees_of_freedom)


def paired_t_test(first: Sequence[float], second: Sequence[float]) -> tuple[float, float, float]:
    if len(first) != len(second):
        raise ValueError("Paired groups must have equal sample sizes")
    differences = tuple(a - b for a, b in zip(first, second))
    t_statistic = mean(differences) / (math.sqrt(sample_variance(differences)) / math.sqrt(len(differences)))
    degrees_of_freedom = float(len(differences) - 1)
    return t_statistic, degrees_of_freedom, two_sided_t_p(t_statistic, degrees_of_freedom)


def linear_regression(x_values: Sequence[float], y_values: Sequence[float]) -> dict[str, float]:
    if len(x_values) != len(y_values) or len(x_values) < 3:
        raise ValueError("Regression requires equal groups with at least three observations")
    x_mean = mean(x_values)
    y_mean = mean(y_values)
    sum_xx = math.fsum((value - x_mean) ** 2 for value in x_values)
    sum_yy = math.fsum((value - y_mean) ** 2 for value in y_values)
    sum_xy = math.fsum((x - x_mean) * (y - y_mean) for x, y in zip(x_values, y_values))
    slope = sum_xy / sum_xx
    intercept = y_mean - slope * x_mean
    r_value = sum_xy / math.sqrt(sum_xx * sum_yy)
    r_squared = r_value**2
    degrees_of_freedom = len(x_values) - 2
    t_statistic = r_value * math.sqrt(degrees_of_freedom / (1.0 - r_squared))
    return {
        "n": float(len(x_values)),
        "slope": slope,
        "intercept": intercept,
        "r": r_value,
        "r_squared": r_squared,
        "p_value": two_sided_t_p(t_statistic, degrees_of_freedom),
    }


def matches_reported_p(p_value: float, reported: str) -> bool:
    if reported.startswith("<"):
        return p_value < float(reported[1:])
    decimals = len(reported.partition(".")[2])
    return abs(p_value - float(reported)) < 0.5 * 10 ** (-decimals)


def _comparison_rows(data_dir: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    test_functions: dict[str, Callable[[Sequence[float], Sequence[float]], tuple[float, float, float]]] = {
        "Welch t test": welch_t_test,
        "Paired t test": paired_t_test,
    }
    for spec in COMPARISONS:
        columns = read_columns(data_dir / spec.source_file)
        first = columns[spec.group_1]
        second = columns[spec.group_2]
        statistic, degrees_of_freedom, p_value = test_functions[spec.test](first, second)
        rows.append(
            {
                "figure": spec.figure,
                "source_file": spec.source_file,
                "test": spec.test,
                "group_1": spec.group_1,
                "n_1": len(first),
                "mean_1": mean(first),
                "sem_1": sem(first),
                "group_2": spec.group_2,
                "n_2": len(second),
                "mean_2": mean(second),
                "sem_2": sem(second),
                "statistic": statistic,
                "df": degrees_of_freedom,
                "p_value": p_value,
                "reported_p": spec.reported_p,
                "status": "PASS" if matches_reported_p(p_value, spec.reported_p) else "FAIL",
            }
        )
    return rows


def _write_statistics(path: Path, regression: dict[str, float], rows: list[dict[str, object]]) -> None:
    fields = (
        "figure", "source_file", "test", "group_1", "n_1", "mean_1", "sem_1",
        "group_2", "n_2", "mean_2", "sem_2", "statistic", "df", "p_value",
        "reported_p", "status",
    )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        regression_row: dict[str, object] = {
            "figure": "Figure 1D",
            "source_file": "manual vs ai all exclude model.csv",
            "test": "Linear regression",
            "group_1": "Manual",
            "n_1": int(regression["n"]),
            "group_2": "AI",
            "statistic": regression["r_squared"],
            "p_value": regression["p_value"],
            "reported_p": "R^2=0.5819; p<0.0001",
            "status": "PASS" if round(regression["r_squared"], 4) == 0.5819 and regression["p_value"] < 0.0001 else "FAIL",
        }
        writer.writerow(regression_row)
        writer.writerows(rows)


def _write_regression_svg(path: Path, x_values: Sequence[float], y_values: Sequence[float], regression: dict[str, float]) -> None:
    width, height = 760, 560
    left, right, top, bottom = 86, 34, 44, 76
    plot_width = width - left - right
    plot_height = height - top - bottom
    x_min, x_max = min(x_values), max(x_values)
    y_min, y_max = min(y_values), max(y_values)
    x_pad = (x_max - x_min) * 0.06
    y_pad = (y_max - y_min) * 0.06
    x_min, x_max = x_min - x_pad, x_max + x_pad
    y_min, y_max = y_min - y_pad, y_max + y_pad

    def sx(value: float) -> float:
        return left + (value - x_min) / (x_max - x_min) * plot_width

    def sy(value: float) -> float:
        return top + plot_height - (value - y_min) / (y_max - y_min) * plot_height

    points = "\n".join(
        f'<circle cx="{sx(x):.2f}" cy="{sy(y):.2f}" r="3.2" fill="#2563eb" fill-opacity="0.62"/>'
        for x, y in zip(x_values, y_values)
    )
    line_y_1 = regression["intercept"] + regression["slope"] * x_min
    line_y_2 = regression["intercept"] + regression["slope"] * x_max
    annotation = html.escape(f'n = {int(regression["n"])}; R² = {regression["r_squared"]:.4f}; p < 0.0001')
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<rect width="100%" height="100%" fill="white"/>
<text x="{width / 2}" y="28" text-anchor="middle" font-family="Arial, sans-serif" font-size="20" font-weight="700">Figure 1D reproduction</text>
<line x1="{left}" y1="{top + plot_height}" x2="{left + plot_width}" y2="{top + plot_height}" stroke="#111827" stroke-width="2"/>
<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_height}" stroke="#111827" stroke-width="2"/>
{points}
<line x1="{sx(x_min):.2f}" y1="{sy(line_y_1):.2f}" x2="{sx(x_max):.2f}" y2="{sy(line_y_2):.2f}" stroke="#dc2626" stroke-width="3"/>
<text x="{left + 16}" y="{top + 28}" font-family="Arial, sans-serif" font-size="16">{annotation}</text>
<text x="{left + plot_width / 2}" y="{height - 24}" text-anchor="middle" font-family="Arial, sans-serif" font-size="17">Manual mean feeding-bout duration (s)</text>
<text x="24" y="{top + plot_height / 2}" transform="rotate(-90 24 {top + plot_height / 2})" text-anchor="middle" font-family="Arial, sans-serif" font-size="17">AI mean feeding-bout duration (s)</text>
</svg>
'''
    path.write_text(svg, encoding="utf-8")


def reproduce(data_dir: Path, output_dir: Path) -> dict[str, object]:
    actual_files = {path.name for path in data_dir.glob("*.csv")}
    expected_files = set(SOURCE_FILES)
    if actual_files != expected_files:
        missing = sorted(expected_files - actual_files)
        extra = sorted(actual_files - expected_files)
        raise ValueError(f"Unexpected source-data set; missing={missing}, extra={extra}")

    regression_columns = read_columns(data_dir / "manual vs ai all exclude model.csv")
    x_values = regression_columns["X (Manual)"]
    y_values = regression_columns["Y (AI)"]
    regression = linear_regression(x_values, y_values)
    comparison_rows = _comparison_rows(data_dir)
    regression_pass = round(regression["r_squared"], 4) == 0.5819 and regression["p_value"] < 0.0001
    all_pass = regression_pass and all(row["status"] == "PASS" for row in comparison_rows)

    output_dir.mkdir(parents=True, exist_ok=True)
    _write_statistics(output_dir / "manuscript_statistics.csv", regression, comparison_rows)
    _write_regression_svg(output_dir / "Figure_1D_regression.svg", x_values, y_values, regression)
    with (output_dir / "source_checksums.sha256").open("w", encoding="utf-8") as handle:
        for filename in SOURCE_FILES:
            digest = hashlib.sha256((data_dir / filename).read_bytes()).hexdigest()
            handle.write(f"{digest}  {filename}\n")

    report_lines = [
        "# Reproduction report",
        "",
        f"Overall status: {'PASS' if all_pass else 'FAIL'}",
        "",
        f"Figure 1D: n={int(regression['n'])}, R^2={regression['r_squared']:.4f}, p={regression['p_value']:.6g} ({'PASS' if regression_pass else 'FAIL'})",
        "",
        "| Figure | Test | Computed p | Reported p | Status |",
        "|---|---|---:|---:|---|",
    ]
    for row in comparison_rows:
        report_lines.append(
            f"| {row['figure']} | {row['test']} | {row['p_value']:.6g} | {row['reported_p']} | {row['status']} |"
        )
    report_lines.extend(
        [
            "",
            "Reproduction covers the temporal analyses backed by the source files in this release.",
            "",
            "Source values for Figure 3F and Supplementary Figure S2 are held separately and are available from the corresponding author on reasonable request.",
            "",
        ]
    )
    (output_dir / "reproduction_report.md").write_text("\n".join(report_lines), encoding="utf-8")
    return {"status": "PASS" if all_pass else "FAIL", "regression": regression, "comparisons": comparison_rows}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--output-dir", type=Path, default=Path("results"))
    args = parser.parse_args()
    result = reproduce(args.data_dir, args.output_dir)
    regression = result["regression"]
    print(f"{result['status']}: Figure 1D n={int(regression['n'])}, R^2={regression['r_squared']:.4f}")
    print(f"Results written to {args.output_dir.resolve()}")
    return 0 if result["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
