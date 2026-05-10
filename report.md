# Spatiotemporal Signal Propagation Project

**Authors:** Aleksandra Burakowska, Mykyta Khrabust 

---

## Part A: Guided Analysis

### Task A1 — Multi-Mutation Spatiotemporal Comparison

#### Research question

Do different PI3K-AKT pathway mutations alter the strength of spatiotemporal ERK signal propagation?
#### Methods

In Task A1, we modified and ran `compare_spatiotemporal_behavior.py` to compare ERK spatiotemporal signal propagation across five cell lines: WT, AKT1_E17K, PIK3CA_E545K, PIK3CA_H1047R, and PTEN_del. The analysis used `ERKKTR_ratio` as the signal column, with `spatial_radius = 60`, `future_window_frames = 3`, and `jump_quantile = 0.9`. Mutation was used as the grouping variable, and the final analysis included 120 experiment-site blocks.

For each experiment-site block, the script computed the Relative Risk (RR) of a near-future ERK activity jump given current exposure to active neighbouring cells:

$$
RR = \frac{P(\text{future jump} \mid \text{neighbour exposed})}
{P(\text{future jump} \mid \text{not neighbour exposed})}
$$

The script generated two main summary tables. `block_level_summary.csv` contains one row per experiment-site block, including the mutation, jump threshold, number of cells/nodes, spatial and temporal edge counts, exposed and unexposed jump rates, risk difference, and block-level RR. These block-level RR values were used for statistical testing. `group_level_summary.csv` aggregates the block-level results by mutation and reports summary statistics such as mean and median RR, mean risk difference, exposed and unexposed jump rates, and average edge counts. This file was used to summarize propagation strength across mutations and to generate the comparison plot. The analysis settings and included groups were documented in `task_description.json`.

For each mutation, we computed the mean RR and standard error from block-level RR values:

$$
SE = \frac{s}{\sqrt{n}}
$$

where \(s\) is the standard deviation of block-level RR values and \(n\) is the number of analysed blocks.

Each mutant was compared with WT using a two-sided Mann–Whitney U test. We used this test because the comparison involved independent groups of block-level RR values: WT blocks versus blocks from each mutant cell line. The Mann–Whitney U test is non-parametric, so it does not require RR values to follow a normal distribution. In our analysis, the test evaluated whether the distribution of block-level RR values for each mutant differed from the distribution observed in WT.

Since four mutant-vs-WT comparisons were performed, p-values were corrected using the Bonferroni method:

$$
p_{\text{corrected}} = \min(p \cdot 4, 1)
$$

Mutations with corrected \(p < 0.05\) were considered significantly different from WT.
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
