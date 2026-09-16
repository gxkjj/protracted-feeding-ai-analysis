# Source-data map

The files in this directory were exported from GraphPad Prism. Blank cells represent unequal group sizes and are ignored by the reproduction script.

| Paper panel | Source file | Analysis |
|---|---|---|
| Figure 1D | `manual vs ai all exclude model.csv` | Linear regression, n = 170 |
| Figure 2A | `HFD male female.csv` | Welch's t test |
| Figure 2B | `isolation male female.csv` | Welch's t test |
| Figure 2C | `Int-HFD male female manual.csv` | Welch's t test |
| Figure 2D | `isolation male female manual.csv` | Welch's t test |
| Figure 3B | `restraint male female.csv` | Welch's t test |
| Figure 3C | `Restraint male female manual.csv` | Welch's t test |
| Figure 3E | `Restraint Mosapride.csv` | Welch's t test |
| Figure 4E | `HFD rescue.csv` | Welch's t test |
| Figure 4F | `isolation rescue.csv` | Welch's t test |
| Figure 4H | `DAT-Cre.csv` | Paired t test |
| Figure S2A | `HFD male.csv`, `HFD female.csv` | Sex-stratified Welch's t tests |
| Figure S2B | `isolation male.csv`, `isolation female.csv` | Sex-stratified Welch's t tests |
| Figure S2C | `restraint male.csv`, `restraint female.csv` | Sex-stratified Welch's t tests |

Reproduction covers the temporal analyses backed by the source files in this release. Manually recorded data not included in this repository are available from the corresponding author on reasonable request.

The source values for **Supplementary Figures S1 and S5** are in [`supplementary/`](supplementary/README.md). That directory contains the 28 panel CSVs for IBI, Bout Number, Latency and Total feeding time, plus a panel map and mouse-level table. DAT-Cre rows are paired. Total feeding time includes the short gaps bridged within merged bouts.
