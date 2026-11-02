# Tree-based (PICRUSt2-style) AmphiBac predictions applied to real communities

## What was run
- **Placement.** SERDP leopard frogs: 12,105 deblur V4 features, 150 bp. Bullfrog: 3,700 ASVs, 120 bp. Meta-analysis: 11,133 sOTUs at 90 bp, kept if present in >=5 of the 2,411 samples rarefied to 1,000 reads; these sOTUs hold 95.7% of rarefied reads. Each set was aligned with `cmalign --matchonly` to RF00177 and placed with EPA-ng 0.3.8 into the full 7,279-isolate reference tree (phase-0 refpkg). Runtimes at 2–3 threads: SERDP+bullfrog 61 s, meta 123 s.
- **NSTI** (top placement) is high: the median is 0.93 for SERDP+bullfrog and 1.29 for meta. Reads with NSTI <=0.05 make up 9.5% of SERDP features, 6.7% of bullfrog features and 3.9% of meta features. The SERDP well-placed read share correlates with the pipeline cultured share at r = 0.92.
- **M1.** Func of the nearest reference tip to the top placement. M1w weights the nearest tip of every placement with LWR >=0.01 (renormalised).
- **M2.** castor 1.8.7 on the reference tree, with every query grafted as a new tip at its top placement (`gappa examine graft`; exact attachment; query tips have NA states). Methods: emp_prob, subtree_average (identical to emp_prob for a 0/1 trait), mp, and Mk ER/ARD (ARD chosen by AIC: 9042 vs 9440).
  - **Mk is degenerate.** The fitted rates are about 1e5 per substitution/site, so every query prediction collapses to the stationary frequency (ER 0.5; ARD 0.365; SD 0 across features). Discordant labels among identical sequences on near-zero branches saturate the model. Mk was therefore excluded from the tests.
- **M3 and M4 (Gaussian process on the tree).** Covariances are built from exact path distances: BM uses a midpoint root, C = (h_x + h_y − D)/2; stationary OU uses exp(−αD). Variance components were fitted by REML once, on all 5,740 scored isolates. Queries are predicted at the attachment point plus the pendant edge as mean ± predictive SD, with P(inhibitory) = P(score > 80). M4 sets the study effect to 0 ("average study").
  - REML fits (σ² values):
    - BM: σ²_phylo = 1299, σ²_e = 2166.
    - BM+study: σ²_phylo = 936, σ²_study = 1447, σ²_e = 1838 (AIC 49365 vs 50278).
    - OU+study: α = 1.93 (AIC 49331, the best).
    - OU without study: α stayed at its starting value of 1.0, i.e. it is poorly identified.
  - The castor squared-change-parsimony score (unscored references pruned) is reported as M3-SCP.

## SERDP (372 frogs; conditional logit in species-site × season strata, n = 344; OLS load among the 100 infected)
- The previous pipeline result was OR 1.02 (0.72–1.45). No tree metric differs from it:
  - Unadjusted ORs ranged from 0.89 to 1.27, and every 95% CI includes 1.
  - The smallest BH q was 0.98 for infection and 0.72 for load.
  - Adjusting for the well-placed share or the cultured share does not change this.
- Examples: M2 emp_prob OR 1.05 (0.73–1.50); M4 BM+study P(score>80) 1.15 (0.77–1.72); M3 BM P 1.00 (0.65–1.53).
- Load among infected: M3 BM score β = −0.18 (−0.44 to 0.09). Restricting to NSTI <=0.1 gives β = −0.28 (p = 0.065, q = 0.72).
- The naive Mann–Whitney tests show no significant tree metric (q >= 0.37). The pipeline Inhibitory-99% metric gave p = 0.002, q = 0.10.

## Bullfrog (40 frogs, rarefied 1,883, seed 1; 10 ponds)
- No metric survives BH correction.
- Several tree score metrics have a nominal *positive* pond-mean β (higher predicted inhibition goes with more Bd): M4 OU+study P, β = 0.22/SD, p = 0.02, q = 0.20. This disappears after adjusting for the well-placed share (p >= 0.13). The pond-mean analysis has n = 10, so power is low.

## Label-permutation null (200 global + 200 within-study shuffles; fixed placements and GP variance components)
- **Design.** Func was shuffled among the 7,279 isolates. Scores were shuffled separately among the 5,740 scored isolates. Methods re-run: M1, M2 emp_prob (castor) and the M3/M4 BM GP (exact re-prediction via K⁻¹K_rq).
- **The observed SERDP effects sit inside the null.** The share of shuffles with an |log OR| at least as extreme was 33–100% (global) and 37–100% (within-study). For load β the range was 31–96%.
- **Unlike the pipeline, the per-frog tree metrics are not a coverage proxy.** The correlation of the observed per-frog metric with the shuffled-label mean was −0.46 to 0.26 (global) and −0.54 to 0.58 (within-study), against r = 0.976 for the pipeline.
  - The tree metrics are means over placed reads, so they vary between frogs because of the real labels: for M2 emp_prob, 2% of frogs lie above the null 95th percentile and 10% below the 5th.
  - That real-label signal is unrelated to Bd.
- **Bullfrog and meta.** The bullfrog r ranged from −0.31 to 0.87 and the meta r from −0.44 to 0.60. See [null_correlation_all_datasets.csv](null_correlation_all_datasets.csv).

## Interpretation
Placement removes the coverage artefact that dominated the closed-reference metrics: the observed metric no longer tracks the shuffled-label metric. It does not create a Bd signal. Two features of the data explain this:
1. Most community reads sit far from any tested isolate (median NSTI about 1).
2. The trait has weak phylogenetic signal: the Mk rates saturate, and the GP assigns large study (σ² = 1447) and residual (σ² = 1838) variance components.

Every tree-based infection OR is indistinguishable from 1 and from label-shuffled trees.

## Files
Per-feature tables: serdp_feature_predictions.csv, bullfrog_feature_predictions.csv, meta_feature_predictions.csv.
Per-sample tables: serdp_tree_per_sample_metrics.csv, bullfrog_tree_per_sample_metrics.csv.
Tests: tree_metrics_Bd_tests.xlsx (README, SERDP stratified/naive, bullfrog, summary, null).
Null outputs: serdp_tree_null_summary.csv, serdp_tree_null_permutations.csv, serdp_tree_null_per_sample.csv.
Scripts: realdata_placement_scripts.tar.gz (r00–r07).
