# Phylogenetic placement does not rescue AmphiBac function prediction

**Placing reads on a reference tree and predicting traits from it, PICRUSt2-style, performs no better than genus majority, and most variants perform worse.** We did this with EPA-ng placement and castor hidden-state prediction, plus continuous-score and study-aware models. With whole studies held out, no tree-based method beats the genus-majority baseline or majority voting over ≥99% hits (R4); no paired bootstrap CI excludes 0. In the leopard-frog data, no tree-based metric predicts Bd infection or load. The limit is the trait: once study effects are separated out, only about a quarter of the variance in inhibition score is phylogenetic.

## Reference package

- **Isolates and alignment.** 7,279 of the 7,873 tested isolates passed filters (≥300 bp, ≤1% ambiguity codes, aligns to 16S). They were aligned with Infernal cmalign to the Rfam bacterial SSU model RF00177, giving 1,533 match columns. Query fragments aligned to the same model land in the identical frame.
- **Tree.** FastTree (GTR+Γ), then raxml-ng `--evaluate` to optimise branch lengths and the model for EPA-ng.
- **Placement check.** Isolate fragments placed back into the tree landed on a tip with an identical fragment 74–76% of the time, and on a tip of the same genus 100% of the time.

## Benchmark: 30-fold leave-study-out

For each study, its isolates were pruned from the tree. Their V4 (254 bp) and V3–V4 (429 bp) fragments were then placed and predicted from the remaining isolates only.

| Predictor | V4, all isolates (n = 6,372): κ | V4, all isolates: AUC | V4, fixed panel (n = 1,397): κ | V4, fixed panel: AUC |
|---|---|---|---|---|
| Genus majority | 0.21 | 0.64 | 0.45 | 0.71 |
| R4: majority of ≥99% full-coverage hits | 0.13 (0.17 tuned) | 0.62 | 0.32 | 0.71 |
| R1: pipeline default top hit at 99% | 0.08 | 0.54 | — | — |
| M1: nearest reference tip (LWR vote) | 0.11 | 0.57 | 0.26 | 0.67 |
| M2: castor parsimony / empirical probabilities / subtree average | 0.11–0.12 | 0.57 | 0.25–0.26 | 0.67–0.68 |
| M3: Brownian / OU score models | degenerate (≈ prevalence) | ~0.5 | — | — |
| **M4: phylogeny + study Gaussian process, Platt-calibrated** | **0.22 (tuned)** | **0.63** | **0.31 (tuned)** | **0.69** |

κ is at a 0.5 threshold unless marked "tuned"; tuned thresholds were chosen on training folds only. "—" marks values not reported for that cell.

![Method comparison]({{artifact:b95974d0-c4fd-4253-8aee-a4be3b20e4ef}})

**What we learned:**

1. **Discrete tree methods are no better than the pipeline's top hit.** On AUC, the nearest-tip and castor hidden-state predictors (M1, M2) are significantly worse than genus majority (−0.06 to −0.07, CI excludes 0).
2. **The fitted trait models collapse because there is no conserved signal to fit.**
   - The castor Mk model fits transition rates of about 10⁵ per substitution per site.
   - The OU model fits α ≈ 10¹⁰.
   - Pure Brownian motion collapses to the mean.

   Each of these predicts the training prevalence or mean everywhere.
3. **Only the study-aware model (M4) ties the baselines.** It is the improvement that matters: separating study effects. It matches genus majority and R4 (AUC 0.63 vs 0.64 / 0.62 for V4) and gives a usable continuous prediction (Spearman 0.26 with the observed CFS score, 0.36 on the panel), but it does not exceed them.
4. **Variance in the CFS score splits roughly evenly between phylogeny and study** (REML, 5,740 scored isolates):

   | Component | Share of variance |
   |---|---|
   | Phylogeny | 0.26 |
   | Study | 0.25 |
   | Bd genotype | 0.08 |
   | Residual | 0.42 |

   A phylogeny-only model attributes 0.43 to phylogeny, because study differences are absorbed into clades. Adding study improves AIC by about 940.
5. **NSTI is not a useful reliability index here.** 75% of V4 queries sit at NSTI ≈ 10⁻⁴ (an identical sequence is in the tree) and are no more accurate than distant queries. Neither NSTI-based nor confidence-based abstention raised accuracy meaningfully.

## Real data: leopard frogs and bullfrogs

We placed the SERDP (150 bp), bullfrog (120 bp) and meta-analysis (90 bp) features into the full tree, predicted each feature with M1–M4, and built per-sample read-weighted metrics.

- **SERDP infection.** Analysed with a conditional logit within site × season (344 frogs; 34 informative strata). Odds ratios per SD were 0.89–1.27 for every tree metric, and every 95% CI includes 1. For example, M2 empirical probabilities gave 1.05 (0.73–1.50) and M4 gave 1.15 (0.77–1.72), against 1.02 (0.72–1.45) for the pipeline. Minimum BH q was 0.98 for infection and 0.72 for load among the 100 infected frogs. Results were the same after adjusting for coverage.
- **Shuffled-label null** (200 global and 200 within-study shuffles of tip labels). The observed SERDP effects sit inside the null: 33–100% of shuffles are at least as extreme.
- **Tree metrics are not just cultured coverage.** The pipeline metric correlated with its shuffled-label mean at r = 0.976. Tree-based per-frog metrics correlate with theirs at only r = −0.54 to 0.58, so the real labels do move the predictions. That label-driven variation is simply unrelated to Bd.
- **Most community features are far from any tested isolate.** Only 9.5% of SERDP features, 6.7% of bullfrog features and 3.9% of meta-analysis sOTUs have NSTI ≤ 0.05. The read share in well-placed features tracks the pipeline's cultured share (r = 0.92).
- **Bullfrog.** No metric survives BH correction. A nominal pond-level association of M4 with Bd (p = 0.02, q = 0.20) is positive in direction and disappears after adjusting for the well-placed share.

![SERDP forest plot]({{artifact:477b6ba3-3596-4a58-b804-43de3a1b56b0}})

## Conclusions

1. **Tree-based prediction fixes the mechanical problems of closed-reference matching.** It removes arbitrary tie-breaking and partial-coverage hits, and gives every read a probability. But it cannot create signal that the 16S tree does not carry: Bd inhibition is only weakly heritable along the tree.
2. **The most useful improvement was modelling study.** It turns a degenerate trait model into one that matches genus majority. It also shows that a quarter of the variance in inhibition is lab or study context, which any prediction from AmphiBac inherits.
3. **For practical use, report predictions as genus-level priors with uncertainty.** If a model is wanted, the study-aware phylogenetic model (M4) is the defensible choice. Validate against outcomes; in these data no method tracked infection.

## Caveats

- Variance components and EPA-ng model parameters were estimated once on all data and reused in every fold. No held-out labels were used in the predictions.
- Castor predictions at a placement use interpolation between the parent and child nodes, not exact re-attachment of the query as a new tip.
- An assay random effect was not identifiable, because all scored isolates come from 96-well CFS assays.
- The real-data predictions used the top placement only, not weighting across alternative placements.
- About 11% of SERDP reads belong to features with no sequence; these are excluded from the per-sample means.
- The meta-analysis was restricted to sOTUs present in at least 5 samples, and was used only for the null comparison.

## Files

**Reference and benchmark**
- [amphibac_placement_refpkg.tar.gz]({{artifact:e25b88fc-ad8c-4ed1-9cde-2b0c12cead2f}}) — alignment, tree, model, tip table, scripts
- [placement_benchmark_report.md]({{artifact:ff9e3dda-c0a3-490f-af17-5c03e9b345c8}}) — detailed benchmark write-up
- [amphibac_placement_benchmark.xlsx]({{artifact:e346d073-701c-4a08-b7a1-1647cd5c3631}})
- [per_isolate_predictions.csv]({{artifact:b342c435-1325-4a4f-b3d1-9642350b10b9}})
- [variance_components_full.csv]({{artifact:b05820a3-3926-4138-b1d1-f0abd541a338}})
- Calibration, NSTI and abstention figures: [V4]({{artifact:4013a01e-1a0d-4b48-ba76-3336c4c3f474}}), [V3–V4]({{artifact:788a20a8-775c-4b5e-95db-c36a0eaf3030}})
- [placement_benchmark_scripts.tar.gz]({{artifact:099028bb-74f9-458f-b493-282a1c03cac6}})

**Real data**
- [realdata_placement_report.md]({{artifact:6d7827cb-2745-42b0-8bf5-4b6e8416f8d2}}) — detailed real-data write-up
- [tree_metrics_Bd_tests.xlsx]({{artifact:ad10b9fc-a773-411c-befb-e00a329b1a5a}})
- [serdp_feature_predictions.csv]({{artifact:1690958f-9666-4e63-9434-728f01829477}})
- [serdp_tree_per_sample_metrics.csv]({{artifact:8984052b-7f58-4264-9824-d8fef34a16ce}})
- [null distributions figure]({{artifact:22103c1f-d5d8-49fa-a5f4-f715a86e335b}})
- [realdata_placement_scripts.tar.gz]({{artifact:3e747179-99c3-4460-962b-c6cf724d98b1}})
