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

#### Results

The full comparison table was saved as `outputs/mutations_comparison_table.csv`.

All analysed cell lines showed mean Relative Risk values above 1. WT had a mean RR of 1.74. Among the mutants, PIK3CA_H1047R showed the highest mean RR of 3.20 and was significantly different from WT after Bonferroni correction (\(p = 1.23 \times 10^{-8}\)). PIK3CA_E545K had a mean RR of 1.71 and was not significantly different from WT (\(p = 0.73\)). AKT1_E17K and PTEN_del showed lower mean RR values than WT, with mean RR values of 1.58 and 1.56, respectively. Both were significantly different from WT after correction..

![Mean ERK spatiotemporal propagation across mutations](outputs/mutations_barplot.png)

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
