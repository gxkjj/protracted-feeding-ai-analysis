# Protracted feeding: source data and reproduction

Minimal source-data and one-command reproduction release:

> Automated Temporal Analysis Reveals Stress-Induced Protracted Feeding Regulated by Mesolimbic Dopamine in Mice

## One-command reproduction

Python 3.9 or newer is the only requirement. No third-party packages are used.

```bash
./reproduce.sh
```

The command runs the tests, recalculates the manuscript statistics, and writes the following files to `results/`:

- `reproduction_report.md`: pass/fail comparison with the paper
- `manuscript_statistics.csv`: sample sizes, means, SEMs, test statistics, degrees of freedom, and p values
- `supplementary_statistics.csv`: the 28 comparisons in Supplementary Figures S3 and S5
- `Figure_1D_regression.svg`: manual-versus-AI regression check
- `source_checksums.sha256`: checksums for all released CSV files

The Figure 1D reproduction uses 170 mice and gives **R² = 0.5819, p < 0.0001**.

## Scope

Reproduction covers the temporal analyses backed by the source files in this release. Source values for Figure 3F and Supplementary Figure S2 are held separately and are available from the corresponding author on reasonable request. See [`data/README.md`](data/README.md) for the complete paper-to-file map.

Supplementary Figures S3 and S5 are included in [`data/supplementary/`](data/supplementary/README.md): IBI, Bout Number, Latency and Total feeding time. Bouts are merged across consecutive video files when gaps are shorter than 7 s. Total feeding time is the sum of merged-bout durations, including bridged gaps. Figure S4 food-intake data and reviewer-only additional analyses are not included in this release.

## Citation

GitHub: https://github.com/gxkjj/protracted-feeding-ai-analysis

Zenodo: https://doi.org/10.5281/zenodo.21775569
