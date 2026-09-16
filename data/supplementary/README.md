# Supplementary Figures S1 and S5: source data

The 28 two-column CSV files contain the individual values used in the four metric columns of Supplementary Figures S1 and S5. Values and English group labels follow the final figure source tables. Each value is one mouse; DAT-Cre has two conditions for each of 12 paired mice. Blank trailing cells in independent-group tables indicate unequal group sizes and are not zeros.

| Figure row | Cohort | n per group |
|---|---|---|
| S1A | HFD | 32, 31 |
| S1B | Isolation, males and females pooled | 32, 31 |
| S1C | Restraint | 23, 21 |
| S1D | Mosapride | 14, 15 |
| S5A | HFD rescue | 7, 7 |
| S5B | Isolation rescue | 6, 7 |
| S5C | DAT-Cre | 12 paired |

Columns within each figure row are **IBI, Bout Number, Latency, Total feeding time**, in that order. `figure_map.csv` specifies each panel file, metric, unit, group, sample size and test. `mouse_values.csv` links all four temporal endpoints to mouse IDs and preserves DAT-Cre pair IDs and row order. Mouse-condition records total 250, corresponding to 238 mice because the DAT-Cre mice contribute both conditions.

## Endpoint definitions

Sequential video files for each mouse were placed on a continuous session timeline before merging adjacent feeding bouts separated by gaps strictly shorter than 7 seconds, including gaps across video-file boundaries.

- **Mean IBI (s):** the mean interval between consecutive merged bouts, including intervals across video boundaries.
- **Bout Number:** the number of merged bouts in the session.
- **Latency (s):** time from session start to the onset of the first feeding bout.
- **Total feeding time (s):** the sum of merged-bout durations, **including short gaps bridged by the merging rule**. This endpoint measures time occupied by merged bouts, rather than feeding-positive frames alone.

Bout start and end times retain the values in the detailed event tables. Segment offsets use recorded frame counts and the exact video frame rate, **30000/1001 fps (approximately 29.97 fps)**, for all 877 segments. Times are in seconds; no conversion to 30 fps or additional rounding is applied. The classifier-training mouse is excluded from the HFD analysis.

## Statistics and reuse

Run `./reproduce.sh` from the repository root to calculate sample sizes, means, SEMs and two-sided P values in `results/supplementary_statistics.csv`. Independent groups use Welch's t test; DAT-Cre uses a paired t test, matching observations on the same row. The 28 P values are unadjusted secondary comparisons and are descriptive.

These are mouse-level figure source data. The reproduction script recalculates statistics from these endpoints; it does not rerun the classifier or reconstruct bouts from video. Manually recorded data not included in this repository are available from the corresponding author on reasonable request.
