# Spatiotemporal Signal Propagation Project

**Authors:** Aleksandra Burakowska, Mykyta Khrabust 

---

## Part A: Guided Analysis

### Task A1 — Multi-Mutation Spatiotemporal Comparison

#### Research question

Do different PI3K-AKT pathway mutations alter the strength of spatiotemporal ERK signal propagation?

#### Methods

In Task A1, we modified and ran compare_spatiotemporal_behavior.py to compare ERK spatiotemporal signal propagation across five cell lines: WT, AKT1_E17K, PIK3CA_E545K, PIK3CA_H1047R, and PTEN_del. We used ERKKTR_ratio as the signal column, with spatial_radius = 60, future_window_frames = 3, and jump_quantile = 0.9.

For each experiment-site block, the script computed the Relative Risk (RR) of a near-future ERK activity jump given current exposure to active neighbouring cells:

$$
RR = \frac{P(\text{future jump} \mid \text{neighbour exposed})}
{P(\text{future jump} \mid \text{not neighbour exposed})}
$$

The block-level RR values were aggregated by mutation. Mean RR and standard error were computed for each mutation:

$$
SE = \frac{s}{\sqrt{n}}
$$

where \(s\) is the standard deviation of block-level RR values and \(n\) is the number of analysed blocks. Each mutant was compared with WT using a two-sided Mann–Whitney U test. P-values were corrected for four comparisons using Bonferroni correction:

$$
p_{\text{corrected}} = \min(p \cdot 4, 1)
$$
#### Output files

The A1 analysis produced three main output files in `analysis_outputs_A1/comparison_mutation_ERKKTR_ratio/`.

`task_description.json` contains metadata describing the analysis setup. It records the main comparison question, the grouping variable, the analysed signal, parameter values, and the groups included in the analysis. In this task, the analysis compared mutations using `ERKKTR_ratio`, with `spatial_radius = 60` and `future_window_frames = 3`. The file confirms that 120 experiment-site blocks were analysed and that the included groups were WT, AKT1_E17K, PIK3CA_E545K, PIK3CA_H1047R, and PTEN_del.

`block_level_summary.csv` contains the detailed results for each individual experiment-site block. Each row corresponds to one analysed block and includes the experiment ID, site ID, mutation, signal column, jump threshold, number of cells/nodes, number of spatial and temporal edges, exposed and unexposed jump rates, risk difference, and block-level Relative Risk. This file was used for statistical testing, because Mann–Whitney U tests should compare the distributions of block-level RR values rather than only group means.

`group_level_summary.csv` contains the aggregated mutation-level results. It summarizes the block-level outputs by mutation and reports the number of analysed blocks, number of unique sites and experiments, total number of nodes, mean and median Relative Risk, mean risk difference, mean exposed and unexposed jump rates, and average numbers of spatial and temporal edges. This file was used to identify the overall propagation strength for each mutation and to generate the group-level comparison plot.
#### Results

The full comparison table was saved as `outputs/mutations_comparison_table.csv`.

All analysed cell lines showed mean Relative Risk values above 1. WT had a mean RR of 1.74. Among the mutants, PIK3CA_H1047R showed the highest mean RR of 3.20 and was significantly different from WT after Bonferroni correction $$(\(p = 1.23 \times 10^{-8}\))$$. PIK3CA_E545K had a mean RR of 1.71 and was not significantly different from WT $$(\(p = 0.73\))$$. AKT1_E17K and PTEN_del showed lower mean RR values than WT, with mean RR values of 1.58 and 1.56, respectively. Both were significantly different from WT after correction..

![Mean ERK spatiotemporal propagation across mutations](outputs/mutations_barplot.png)
**Figure 1.** Mean Relative Risk (RR) of ERK spatiotemporal propagation across PI3K-AKT pathway mutations. Bars show the mean block-level RR for each mutation, and error bars show the standard error. The dashed horizontal line marks \(RR = 1\), which corresponds to no neighbour-associated increase in the probability of a future ERK jump. Values above 1 indicate that cells exposed to active neighbouring cells were more likely to show a near-future ERK activity jump. Asterisks indicate mutations that were significantly different from WT after Bonferroni correction.
#### Interpretation

TODO: write 150–250 words after checking exact values from `mutations_comparison_table.csv`.

---

### Task A2 — Lagged Exposure Analysis Across Mutations

#### Research question

Do different mutations exhibit different spatiotemporal relay timescales?

#### Methods

TODO.

The lag-specific Relative Risk was defined as:

\[
RR(\tau) = \frac{P(\text{jump at } t+\tau \mid \text{neighbour exposed at } t)}{P(\text{jump at } t+\tau \mid \text{not exposed at } t)}
\]

The optimal lag was defined as:

\[
\tau^* = \arg\max_{\tau} RR(\tau)
\]

#### Results

TODO: insert lagged exposure plot and table.

#### Interpretation

TODO.

---

### Task A3 — Parameter Robustness Assessment

#### Research question

How sensitive is the Relative Risk metric to analysis parameter choices?

#### Methods

TODO.

#### Results

TODO: insert robustness plot.

#### Recommendation

TODO.

---

## Part B: Independent Research

### Research question and hypothesis

TODO.

### Methods

TODO.

### Results

TODO.

### Discussion

TODO.

---

## References

TODO.
