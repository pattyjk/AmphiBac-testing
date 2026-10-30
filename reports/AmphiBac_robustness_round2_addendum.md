# AmphiBac robustness: round-2 addendum

This addendum follows the round-1 report, [AmphiBac_prediction_robustness_report.md]({{artifact:68283292-5522-4f52-929b-09571edb7201}}). It covers four tests:

1. A label-permutation null and a genus-only baseline.
2. Phylogenetic signal of Bd inhibition.
3. Stability across database releases, including predicting isolates added later.
4. Validation against Bd load in the bullfrog nitrogen-enrichment experiment.

## Summary

| Test | Result | What it means |
|---|---|---|
| **Permutation null** (1,000 label shuffles; hits fixed) | Observed per-sample Inhibitory_Proportion correlates with its shuffled-label mean at ρ = 0.92 (global shuffle), 0.94 (within genus) and 0.97 (within study). Only 1.9% of 2,411 meta-analysis samples exceed the null 95th percentile, against 5% expected by chance. 40% fall *below* the 5th percentile. For bullfrogs, 0/40 samples exceed it under any shuffle. | **The pipeline metric carries no information from the inhibition labels.** Random labels reproduce it, because it measures how much of the community has a cultured relative in the database (null mean vs Full_Proportion ρ = 0.93). |
| Null, majority-vote metric | 52% of meta samples beat the global null, but only 12% beat the within-genus null and 24% the within-study null. | Label consensus carries some signal, and most of it is genus-level or study-level. |
| Null, meta-analysis contrasts | Adult vs larvae and Madagascar vs USA come out significant in **100% of shuffles**. The observed Madagascar–USA effect (r = 0.31) is *below* the null median (0.50). *A. boreas* vs *L. catesbeianus* is significant in 13–18% of shuffles. | Under the pipeline metric, those differences come from cultured coverage and would appear with any labels. They say nothing about inhibitory function. |
| **Genus-only baseline** (SINTAX/RDP v18, 7,316/8,739 isolates assigned) | Cohen's κ for genus-majority prediction, leave-self-out / leave-study-out: **0.30 / 0.21** at 120 bp and **0.35 / 0.25** at full V4. The best 99% identity rule (majority vote) gives 0.32 / 0.12 and 0.32 / 0.11. Genus minus majority-vote κ, leave-study-out: +0.09 (95% CI 0.03–0.18) at 120 bp and +0.14 (0.05–0.26) at V4. | **Knowing only the genus predicts as well as 99% identity matching, and better across studies.** Fine-scale V4 matching adds no transferable information. |
| Within-genus heterogeneity | 46 of 50 genera with ≥20 tested isolates are 10–90% inhibitory. *Pseudomonas* is 60% inhibitory overall (n = 967, 26 studies), but its per-study share ranges from 0% to 100%. | Genus gives a prior, not a prediction. |
| **Function-decay curve** (5,378 V4 isolates, 14.5M pairs) | For identical-V4 pairs, κ beyond prevalence is 0.30 (0.20–0.42) within a study but 0.09 (−0.01 to 0.22) between studies. Between studies, P(partner inhibitory \| focal inhibitory) is a flat 0.52–0.57 from 0 to 10 V4 differences, then drops to the random level (~0.38) at ≥11. | V4 identity adds nothing within ~10 differences; only a shallow, genus-like lineage signal exists. Most identical-V4 agreement is a same-lab effect. |
| Full-length 16S (1,884 isolates, ≥1,000 bp shared) | Among between-study pairs with identical V4, P(partner inhibitory) is 0.78 when full-length identity is ≥99.5% and 0.48 when it is <99%. | **Near-full-length identity does carry signal that V4 cannot see.** This subset is provisional. |
| **Phylogenetic signal** | consenTRAIT τD = 0.0053 subst/site. That is deeper than global shuffles (p = 0.002) but not deeper than within-study shuffles (p = 0.71). Fritz–Purvis D = 0.58–0.76; typical microbial binary traits average about 0.07. Blomberg's K on the CFS score is ~3×10⁻⁵ and Pagel's λ is 0.54–0.71. | Inhibition is weakly clumped and shallow. The signal is far below that of typical conserved microbial traits, and the τD signal disappears once study is controlled. |
| **Release stability** | Full and Facilitating outputs are stable (ρ ≥ 0.999). Inhibitory changes little from 2023.2 to 2025.1 (ρ = 0.995; 3.6% of samples shift >0.05). **InhibitoryStrict from 2020.6 to 2025.1 is unstable**: ρ = 0.927, 26% of samples shift >0.05 and 22% change quartile. The cause is a change in how the subset is defined (550 isolates added, 88 removed); only 1 of 7,382 shared isolates changed its function label. Between 2020.6 and 2023.2, 11% of matched meta-analysis sOTUs changed their top-hit label (13% of matched reads), and only 3 of those changes involved new isolates. | Results published with different releases are not comparable for the Strict subset. Top-hit labels also change between releases through tie-breaking, not new data. |
| Release-dependent contrast | *A. boreas* vs *L. catesbeianus* is r = −0.20 (p = 0.0004) with 2020.6 Strict and n.s. (r = +0.05, p = 0.28–0.40) with later releases. Essentially one reference (AmphiBac_3521) drives the difference. | One reference can flip a published host comparison. |
| **Prospective prediction** (isolates added later, predicted from the older release) | From 2020.6 at 120 bp (n = 593): κ = 0.03 (top hit), 0.01 (Strict), 0.03 (majority). From 2023.2 (n = 239): κ between −0.14 and −0.03. Every rule, cohort and fragment length gives κ ≤ 0.15. Predicting the same isolates from the full 2025.1 database is no better (κ −0.04 to 0.15). | **True out-of-sample prediction is at chance.** Older releases are not missing the right references; the information is not in V4. |
| **Bd load, bullfrogs** (40 frogs, 10 ponds) | Bd load is completely separated by N treatment. No metric correlates with LogBd in any analysis: naive, within-treatment, pond-centred, mixed model, or pond means with treatment (15 metrics; minimum BH q = 0.29). All 15 naive correlations are *positive*: more "inhibitory" reads go with more Bd. Inhibitory and Full proportion correlate at r = 0.996 across frogs. With 10 ponds, 80% power requires a partial \|r\| ≥ 0.80. | No support for the metric as an infection predictor, though power is low. "No evidence" here is not "evidence of absence". |

## Bullfrog metadata issues

These come from `nut_enrich_16S_map.txt`; please check them against your original file.

- **`NitrateuM` is misaligned.** It equals the column sorted in ascending order. It takes 3–4 different values within every pond, and gives 4.8 µM to frogs in both treatments. Its correlation with LogBd (ρ = 0.995) is purely a sorting artefact.
- **Bd load and treatment are fully confounded.** Low N ranges from 25.9 to 1,066 ZE and High N from 20,914 to 506,726 ZE.
- **Six mapped frogs have no 16S data** (P10F1, P2F1, P2F2, P2F4, P2F5, P8F5), and two have fewer than 10 reads (P10F3, P4F1). The Low N ponds are therefore represented by only 1–4 frogs each; P2 has one.

## Figures

**Label-permutation null.**

![Permutation null]({{artifact:c97dafc4-3662-4743-acac-6b61fe7b2a71}})

Contrasts under shuffled labels: [fig_null_contrasts.png]({{artifact:3cd0d34f-288f-406a-8851-1d96581012a3}}).

**Genus-only baseline vs 99% identity rules, and label mixing within genera.**

![Genus baseline]({{artifact:b345b5be-def9-4292-912e-638073988986}})

**Function-decay curves**, within vs between study, for V4 and full-length 16S.

![Decay curves]({{artifact:3b809a47-3603-487e-a1a0-2995d18e7e88}})

Tree-based signal against nulls: [fig_r2_tree_signal_nulls.png]({{artifact:064eb4af-e007-40e3-812a-d496b9c29672}}).

Release stability: [fig_r2_release_stability.png]({{artifact:af2e749b-e7ea-459a-b3f2-23a26b0e7d84}}). Prospective prediction: [fig_r2_prospective.png]({{artifact:a7d24475-bcd3-47fe-8384-9dc30b59c168}}).

**Bd load vs AmphiBac metrics** (bullfrog).

![Bd validation]({{artifact:a3094ba4-400d-4264-9553-942866638a1b}})

## Updated recommendations

These add to the round-1 list.

1. **Do not interpret per-sample "inhibitory proportion" from subset mode as function.** It is indistinguishable from random labels. Any functional claim should be tested against a label-permutation null, ideally within genus and within study.
2. **Benchmark against a genus-only baseline.** If a 16S-matching method cannot beat genus-majority prediction leave-study-out, it adds nothing. None of the 99% V4 rules beat it.
3. **If identity matching is used, prefer near-full-length 16S** (≥99.5%) over V4. That is the only setting where agreement between studies rose clearly above the genus-level rate (provisional, n = 1,884).
4. **State the database release and subset.** InhibitoryStrict results from 2020.6 are not comparable with 2023.2/2025.1.
5. **Validate against outcomes.** In the one dataset with Bd loads, no metric tracked infection. An adequately powered test needs many more independent units (ponds or sites) than 10.

## Caveats and declared limits

- **Permutation null.**
  - Subset mode is approximated by "any ≥99% full-coverage hit is Inhibitory". This reproduces the pipeline for 276 of its 291 inhibitory features; the per-sample difference is at most 0.036 (meta) and 0.007 (bullfrog).
  - A single rarefaction draw was used.
  - 1,423 isolates unassigned at genus level were pooled as one group.
- **Genus baseline.** RDP v18 SINTAX was used at a 0.8 cutoff. Unassigned isolates were predicted non-inhibitory, and scoring is binary only.
- **Phylogenetic signal.**
  - Single ML FastTree trees were used, with no bootstrap trees.
  - consenTRAIT was re-implemented in R.
  - D and the unique-V4 τD use one random isolate per identical-V4 group.
  - Full-length results rest on 1,884 isolates of ≥1,200 bp and are provisional.
  - Pair-level proportions are descriptive, because pairs are not independent.
- **Release stability.**
  - The prospective tests have only 5 (after 2020.6) and 2 (after 2023.2) study clusters, so the cluster CIs are coarse.
  - About half of the added isolates lack a 515F site. This includes the whole 522-isolate 8F/1429R Brazil batch.
  - 2020.6 reference labels were taken from the FASTA header suffix.
- **Species contrasts.** The species contrast uses the two most-sampled *named* species; "Boophis sp" (n = 170) was skipped as a genus-level label. The contrasts ignore study and site confounding, so they illustrate method sensitivity, not biology.
- **Bd validation.**
  - The pond variance component was estimated at zero, so the mixed models reduce to OLS plus treatment.
  - Frog-level power ignores clustering, so it is optimistic.

## Files

- Null and genus: [addendum]({{artifact:aebfc937-2553-4929-b296-f926de1e067f}}), [summary xlsx]({{artifact:425f654c-4028-4320-839e-98cb6910588d}}), [genus_assignments.csv]({{artifact:5f91e1d1-4bab-4e34-903c-fda11446522b}}), [null_per_sample_summary.csv]({{artifact:a8c0fc3b-a34f-4fcc-bca7-86c21b69d51b}}), [null_contrasts.csv]({{artifact:269c14a4-41f7-4e77-ace6-db3f9306f8ea}}), [genus_baseline_metrics_long.csv]({{artifact:e39e1fd6-97a0-4b2c-a2a8-5bd67acb364a}}), [scripts]({{artifact:eeb11567-e054-4544-8bcf-41fe7a372658}}).
- Phylogenetic signal: [addendum]({{artifact:295bf8d0-4475-4e8f-80ff-5dcd2210d5a6}}), [summary xlsx]({{artifact:9ba920fb-3c5b-45ac-82ee-3582c75f5193}}), [decay_v4_bins.csv]({{artifact:2a4c0164-478e-421d-90f1-b22d238c7398}}), [V4 tree]({{artifact:26c8b5c2-197b-485a-9de1-de4e4446c7d8}}), [full-length tree]({{artifact:72d0c5ee-9894-4bb9-9e01-c9ed845a6602}}), [scripts]({{artifact:90f4dade-a5dc-4dca-860c-a908783232f5}}).
- Release stability: [addendum]({{artifact:494e6896-7c42-4b3a-98fd-e3eeae3915cc}}), [summary xlsx]({{artifact:2fcd0bbb-7230-4184-bb96-6f6470fe1b1c}}), [release_comparison_per_sample.csv]({{artifact:a9d0b997-0f9c-4ddf-ad74-4c50eba6a3c3}}), [prospective_metrics_long.csv]({{artifact:34424b76-4d8b-4a9d-aa04-70f65e96d970}}), [relabelling_summary.csv]({{artifact:382f3263-b2d4-422d-9c32-f1f9086ba472}}), [scripts]({{artifact:9bf9b060-cf92-4381-9fe8-84e588464d41}}).
- Bd validation: [note]({{artifact:62b2f0b2-a8a5-4bb7-a173-852f7caf3483}}), [tests xlsx]({{artifact:4cd17f3e-ee0e-486a-9346-e704ebf9c55d}}), [bf_metrics_with_Bd.csv]({{artifact:08691095-b96f-49ef-8a6d-e966332b8b13}}), [bd_power.csv]({{artifact:cb6f2da9-00a2-41aa-bc98-8099cc36d71a}}), [script]({{artifact:24d6f2d9-53cc-4dca-9b29-3e8b2ea48061}}).
