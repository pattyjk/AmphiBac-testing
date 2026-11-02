# AmphiBac: phylogenetic placement + trait prediction vs closed-reference rules (leave-study-out benchmark)

**Bottom line.** No tree-based method beats genus majority or R4 (majority of ≥99 % full-coverage hits), for either fragment or isolate set.
- The discrete hidden-state predictors (M1 nearest tip; M2 parsimony, empirical probabilities, subtree average, Mk) are on par with R1/R2. They are significantly worse than genus majority on AUC: the paired difference is about −0.04 to −0.07, with the 95 % CI excluding 0.
- The study-aware phylogenetic GP (M4), with a Platt mapping fitted on training folds, is the only placement method that ties the genus and R4 baselines (V4, all isolates: AUC 0.63 vs 0.64 / 0.62). Its κ at the training-tuned threshold is 0.22 vs 0.15 / 0.17. Every paired CI includes 0.
- NSTI is not a reliability index here. About 75 % of V4 queries sit at the EPA-ng pendant floor (NSTI ≈ 1e-4), and they are not predicted better than distant queries.

## Design
- **Folds.** 30 leave-study-out folds (study = `Reference` string as stored in AmphiBac 2025.1). Queries are all tested isolates with a fragment: V4_254 n = 6,372; V3V4_429 n = 6,148. This includes 196 / 175 isolates that are not in the reference tree.
- **Fixed panel.** The 1,397-isolate panel (`panel_all10`) is reported separately.
- **Placement.**
  - Fragments were aligned once with `cmalign --matchonly` against RF00177. Alignment is per-sequence against the CM, so it does not depend on the fold.
  - For each fold, the study's tips were pruned from the phase-0 tree and MSA (`placement_lib.prune_reference`). Fragments were then placed with EPA-ng 0.3.8 using the full-tree GTR+G model file.
  - Placements with LWR ≥ 0.01 were kept and renormalised; the top placement was always kept.
- **M1.** Nearest reference tip of each placement. P = LWR-weighted vote; the top-placement variant is also reported.
- **M2.** castor 1.8.7 HSP on the pruned tree:
  - `hsp_max_parsimony` (all_equal, weight_by_scenarios), `hsp_empirical_probabilities`, `hsp_subtree_averaging` (0/1 tip states) and `hsp_mk_model` (ER vs ARD by AIC; ARD won in all 30 folds).
  - **Prediction at the placement point** is a linear interpolation between the castor posteriors of the edge's child and parent nodes, by distal position. This is not an exact re-attachment of the query as a new tip.
  - For Mk, the attachment-point value is then propagated along the pendant edge with the fitted 2-state rate matrix. A variant without pendant propagation is also reported.
  - For binary states, empirical probabilities and subtree averaging are the same quantity (fraction of descendant tips inhibitory), so their results are identical.
- **M3/M4.** Exact Gaussian conditioning on the tree covariance (BM: C = shared root-to-MRCA path length; OU: s²·exp(−α·patristic distance)), evaluated at each placement point.
  - Covariance with tips below the placement edge's child = attachment-point depth; covariance with other tips = C(any tip under child, ·). Query variance = attachment depth + pendant.
  - Mean = LWR mixture over placements. Variance = top-placement conditional variance + between-placement variance.
  - Models:
    - BM (pure);
    - BM + residual;
    - OU (pure);
    - OU + residual;
    - M4s = BM-phylo + study + residual;
    - M4 = BM-phylo + study + Bd-genotype (spaces removed) + residual.
  - Assay method was dropped: all 5,740 scored isolates are 96-well CFS assays (the "97…296 well CFS assay" labels are spreadsheet autofill artefacts).
  - Variance components were fitted once by REML on all scored in-tree isolates using the full tree. Fold predictions condition on the fold's training tips with these components fixed (small leakage of variance parameters only).
  - The prediction for a new study sets study and genotype effects to 0 (average study) and puts their variances into the predictive variance. P(inh) = P(score > 80) under the Gaussian predictive.
  - Because that P is poorly ranked (see results), a Platt mapping logistic(a + b·mean) was also fitted per held-out study on the other studies' LOSO predictions ("Platt").
- **Baselines (recomputed on the same folds).** The held-out study's sequences were removed from each database.
  - R1: exact q2-vsearch call (`--id 0.99 --strand both --qmask none --minseqlength 1`, first hit) against FullDatabase.
  - R2: the same call against the Inhibitory subset.
  - R4: exhaustive (`--maxaccepts 0 --maxrejects 0`) hits to tested isolates with full query coverage (qilo = 1, qihi = ql), fraction Inhibitory, majority > 0.5.
  - Genus majority: the fraction inhibitory among other-study tested isolates of the query's full-length RDP genus (optimistic, as in round 3).
  - Unmatched queries get the fold's training prevalence as P and are called not inhibitory.
  - The recomputed V4 κ values match round 3: R1 panel 0.154, R4 panel 0.323, R4 all 0.134.
- **Metrics.**
  - Binary Inhibitory vs not: κ at P > 0.5, and κ at the threshold maximising balanced accuracy on the other studies' predictions.
  - AUC, Brier, and Spearman of predicted (mean score or P) vs observed CFS score.
  - 95 % study-cluster bootstrap CIs (B = 500), with paired bootstrap differences against genus majority and R4.

## Results

![Method comparison]({{artifact:art_39d682fd-623c-4a26-acdc-81b29a6c6469}})

Values are point estimates [95 % study-cluster CI]. κ_t = κ at the training-fold threshold. Full table: [summary_table.csv](/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/artifacts/proj_16c61e02b61f/da953a45-a220-41a4-9287-4c8ea79f3a82/vbdf5eda2_summary_table.csv).

| V4 (254 bp) | all: κ0.5 | all: κ_t | all: AUC | panel: κ0.5 | panel: κ_t | panel: AUC |
|---|---|---|---|---|---|---|
| Genus majority | 0.21 [0.09,0.43] | 0.15 [0.03,0.34] | 0.64 [0.58,0.71] | 0.45 [0.30,0.58] | 0.39 [0.23,0.52] | 0.71 [0.66,0.78] |
| R4 majority full-cov | 0.13 [0.04,0.28] | 0.17 [0.06,0.29] | 0.62 [0.56,0.69] | 0.32 [0.22,0.43] | 0.31 [0.19,0.49] | 0.71 [0.65,0.78] |
| R1 top hit 99 % | 0.08 [0.01,0.19] | 0.08 | 0.54 [0.51,0.59] | 0.15 [0.01,0.28] | 0.15 | 0.58 [0.51,0.64] |
| M1 nearest tip (LWR vote) | 0.11 [0.04,0.21] | 0.11 | 0.57 [0.53,0.63] | 0.26 [0.16,0.34] | 0.26 | 0.67 [0.63,0.72] |
| M2 max parsimony | 0.12 [0.04,0.22] | 0.10 | 0.57 [0.53,0.65] | 0.25 [0.16,0.34] | 0.23 | 0.68 [0.63,0.74] |
| M2 empirical prob. | 0.11 [0.04,0.23] | 0.13 | 0.57 [0.53,0.64] | 0.26 [0.17,0.33] | 0.25 | 0.67 [0.62,0.72] |
| M4 phylo+study+genotype (Platt) | 0.10 [0.03,0.22] | 0.22 [0.13,0.32] | 0.63 [0.55,0.70] | 0.20 [0.11,0.24] | 0.31 [0.20,0.39] | 0.69 [0.62,0.73] |

V3–V4 (429 bp) gives the same picture:
- **All isolates:** genus 0.21 / 0.13 / AUC 0.64; R4 0.11 / 0.20 / 0.60; M2 empirical 0.11 / 0.11 / 0.60; M4-Platt 0.10 / 0.23 / 0.63.
- **Panel:** genus 0.45 / 0.33 / 0.71; R4 0.27 / 0.38 / 0.70; M2 parsimony 0.21 / 0.27 / 0.67; M4-Platt 0.21 / 0.33 / 0.69.

**Paired differences vs genus majority (V4, all isolates):**
- AUC: M1 −0.07 [−0.11,−0.03]; M2 empirical −0.06 [−0.10,−0.02]; M2 parsimony −0.06 [−0.10,−0.02]; M4-Platt −0.01 [−0.07,+0.02].
- κ_t: M4-Platt +0.07 [−0.07,+0.12].

M4-Platt vs R4: AUC +0.01 [−0.06,+0.05], κ_t +0.05 [−0.05,+0.12]. No tree method has a paired CI above 0 against either baseline in any window × set cell ([paired_diffs.csv](/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/artifacts/proj_16c61e02b61f/b2258877-c3d9-441e-b595-7251ad658bc5/v8f61cf47_paired_diffs.csv)).

**Degenerate models (an informative failure):**
- **Mk rates.** The fitted Mk rates are about 0.9–1.2 × 10⁵ (non-inh→inh) and 1.8–2.1 × 10⁵ (inh→non-inh) per substitution/site, because identical sequences on ~1e-6 branches disagree.
- **Mk predictions.** Over EPA-ng's minimum pendant of 1e-4, any query is therefore at the stationary frequency. The predictions then equal the fold's training prevalence, which is anticorrelated with the held-out study's own prevalence, giving AUC 0.41 (V4 all).
- **Pure BM and pure OU** are degenerate in the same way. REML OU α ≈ 7 × 10¹⁰, so the OU is white noise; BM σ² explodes.
- **Gaussian P from M3/M4 ranks poorly** (AUC 0.50–0.53), even though the predicted means carry signal: M4 Spearman ρ with the CFS score is 0.26–0.27 for all isolates and 0.35–0.36 for the panel. The predictive SD (about 60 score units) varies between queries and scrambles the ranking of P(score > 80). The Platt mapping on the mean restores AUC 0.63 / 0.69.

### Variance components of the CFS score (REML, 5,740 scored isolates, full tree) — key result

Full table: [variance_components_full.csv](/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/artifacts/proj_16c61e02b61f/eb9d41c9-1ce3-4788-ac71-c09c250f5ddf/vb05820a3_variance_components_full.csv).

| model | AIC (REML) | phylogeny | study | Bd genotype | residual |
|---|---|---|---|---|---|
| M4: phylo + study + genotype + e | 59,881.7 | 0.26 (0.34) | 0.25 (0.25) | 0.08 (0.09) | 0.42 (0.32) |
| M4s: phylo + study + e | 59,912.4 | 0.29 (0.38) | 0.26 (0.27) | — | 0.45 (0.35) |
| OU + e | 60,775.0 | 0.22 | — | — | 0.78 |
| BM + e | 60,826.0 | 0.43 (0.56) | — | — | 0.57 |

- **How shares are computed.** Shares are expected sample-variance shares: phylogeny σ²_p·(mean diag C − mean C); study/genotype σ²·(1 − Σ p_k²). The values in parentheses use σ²_p × mean root-to-tip depth.
- **Fitted values.** M4 components: σ²_p = 883 per unit tree depth, σ²_study = 1,422, σ²_genotype = 503, σ²_e = 1,828 (score units²). The fitted OU + e half-life is 0.17 substitutions/site.
- **Interpretation.**
  - Study (lab/assay batch) explains about as much CFS-score variance as phylogeny, and adding it improves the REML AIC by about 900 over a phylogeny-only model.
  - Without the study term, the phylogeny absorbs it: the BM + e "heritability" of 0.43 roughly halves to 0.26 once study and genotype are modelled.
  - The model-implied total variance exceeds the observed score variance (2,561), so treat the shares as relative, model-based quantities.

### NSTI and abstention

![V4 calibration, NSTI strata, abstention]({{artifact:art_ffec1f35-7327-47e4-84b4-043cd05e0589}})
![V3-V4 calibration, NSTI strata, abstention]({{artifact:art_45810d57-24a4-4cf4-9caf-537dd13ec34a}})

**NSTI strata.**
- V4: 4,776 of 6,372 queries have NSTI ≈ 0 (< 5e-4). Their AUC is 0.63 for genus majority and 0.62 for M4-Platt, no better than at NSTI 0.01–0.05 or > 0.1 ([nsti_strata.csv](/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/artifacts/proj_16c61e02b61f/ae6c75f7-e420-45a8-ad92-7cc80caf27b2/v63e61b35_nsti_strata.csv)).
- Keeping only the lowest-NSTI 20–50 % of queries leaves AUC at about 0.52–0.66 for every method. Ties at the NSTI floor were broken at random.
- Confidence-based abstention does not help the calibrated models either: the most confident predictions are confident non-inhibitory calls.
- Raw accuracy by NSTI bin is confounded by prevalence: the far-NSTI bins are only 22–34 % inhibitory, vs 41–43 % at NSTI ≈ 0.

**Calibration.**
- Brier scores of M4-Platt and genus majority are 0.23, and R4's is 0.25 (V4 all).
- M1/M2 are over-dispersed (Brier 0.30–0.34) because tip-vote probabilities are near 0 or 1.
- R1/R2 are binary (Brier 0.43–0.51).

## Interpretation
Placement recovers the genus perfectly (phase 0), but within-genus Bd inhibition is not phylogenetically predictable across studies: identical fragments disagree between studies about as often as genus mates do. The tree adds nothing over the genus label for discrete prediction.

The one useful ingredient is the study random effect. Modelling the continuous score with study partialled out brings tree-based prediction up to the genus and R4 level. It does not exceed that level.

For a production rule, genus majority or R4 with calibrated probabilities is as good as anything tested here and is far simpler.

## Caveats
- Genus majority uses the full-length genus (optimistic upper bound).
- The EPA-ng model parameters and the GP variance components were estimated on the full tree (minor leakage of global parameters only; tip labels/scores of the held-out study are never used).
- The GP root is the basal trifurcation of each pruned tree (not re-estimated).
- castor HSP prediction uses parent/child interpolation, not exact re-attachment.
- Studies are `Reference` strings as stored. Several "unpublished"/"Not published" strings pool unrelated datasets, and three variants of "10.3354/dao02560, unpublished" are separate folds.
- The scores are raw: heavy left tail, min −390. Gaussian models were not robustified.

## Files
- [per_isolate_predictions.csv](/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/artifacts/proj_16c61e02b61f/4ee0aba9-fb54-4847-8ed2-5bbe88f594b2/vb342c435_per_isolate_predictions.csv): one row per isolate × fragment. It contains all methods' P, M3/M4 means and SDs, NSTI, the top-placement LWR, nearest tip, baseline hits, labels, score and panel flag.
- [metrics_tidy.csv](/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/artifacts/proj_16c61e02b61f/c61c58a1-6faf-4853-a914-f09d320ea046/v193f24de_metrics_tidy.csv), [summary_table.csv](/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/artifacts/proj_16c61e02b61f/da953a45-a220-41a4-9287-4c8ea79f3a82/vbdf5eda2_summary_table.csv), [paired_diffs.csv](/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/artifacts/proj_16c61e02b61f/b2258877-c3d9-441e-b595-7251ad658bc5/v8f61cf47_paired_diffs.csv), [nsti_strata.csv](/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/artifacts/proj_16c61e02b61f/ae6c75f7-e420-45a8-ad92-7cc80caf27b2/v63e61b35_nsti_strata.csv), [abstention.csv](/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/artifacts/proj_16c61e02b61f/d8f66312-3e2b-44f8-8745-322e4e17cb93/v561ba057_abstention.csv), [calibration.csv](/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/artifacts/proj_16c61e02b61f/62e1f11e-2d1e-428a-a3a1-aabeb01e0440/v03ab73bb_calibration.csv), [variance_components_full.csv](/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/artifacts/proj_16c61e02b61f/eb9d41c9-1ce3-4788-ac71-c09c250f5ddf/vb05820a3_variance_components_full.csv).
- Workbook with all of the above: [amphibac_placement_benchmark.xlsx](/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/artifacts/proj_16c61e02b61f/4b9ca9dd-1405-4d9b-98df-2d679ca87c30/ve346d073_amphibac_placement_benchmark.xlsx).
- [placement_benchmark_scripts.tar.gz](/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/artifacts/proj_16c61e02b61f/2edaba90-601a-4039-a579-6baf21d38027/v099028bb_placement_benchmark_scripts.tar.gz): scripts 01_loso_place.py, 02_baselines_vsearch.py, 03_varcomp.py, 04_predict_folds.py, 05_score.py, bench_lib.py, hsp_castor.R, plus the fold table and EPA-ng log. It relies on the phase-0 reference package and placement_lib.py.
- [castor_1.8.7_Rlib_osx-arm64.tar.gz](/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/artifacts/proj_16c61e02b61f/d82078e6-ef7a-44e3-946f-3a57cb14b03c/va208f97d_castor_1.8.7_Rlib_osx-arm64.tar.gz): castor built from CRAN source into a workspace library, because r-castor is not on conda-forge.
