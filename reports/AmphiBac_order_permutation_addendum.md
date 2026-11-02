# Addendum: reference-order (tie-break) uncertainty in AmphiBac closed-reference calls

**Scope.** This addendum covers only the uncertainty that comes from the arbitrary order of sequences in the AmphiBac 2025.1 FASTA. It does not cover label/assay uncertainty, rarefaction, or the choice of identity threshold.
All searches used the exact q2-vsearch closed-reference call (vsearch 2.32.0, `--id 0.99 --strand both --qmask none`, default maxaccepts 1 / maxrejects 32, default dust masking of the database).

## 1. Emulator: how vsearch actually breaks ties
- Defaults are maxaccepts 1 and maxrejects 32. vsearch ranks candidate targets by the number of unique 8-mers they share with the query (dust-masked database index), then by **shorter target length**, then by file position. It keeps the top 33 and aligns them in that order. The first target with identity ≥99% is accepted; the search stops after 32 rejects. When hits are found on both strands, the one with the higher exact identity wins.
- Reshuffling the reference file therefore only re-orders targets that tie on k-mer count **and** target length. Most max-identity ties are never reachable (reachable set < max-identity tie set for 57% (meta), 74% (bullfrog) and 64% (SERDP) of matched features). Among order-dependent features the median number of reachable references is 2 (max 15).
- **Validation.** The k-mer emulator was checked against real vsearch runs on shuffled FASTAs with the same seeds: 301 FullDatabase runs per dataset (seed 0 = distributed order, seeds 1–300 shuffled) and 101 Inhibitory runs. It reproduced **100% of per-feature assignments in every run** for all three datasets.
- The naive "uniform among max-identity ties" emulator is strongly biased. It matched the real assignment for only 45–57% of features with the same seed. Chi-square rejected it for 95–97% of multi-reference features (median TVD 0.75), and per-sample means correlated only r = 0.66 (meta), 0.98 (bullfrog) and 0.77 (SERDP). It was **not** used. For comparison, the k-mer emulator had chi-square p<0.05 in 2.7–5.1% of features, which is the nominal rate (median TVD 0.02, Monte-Carlo noise of 300 runs), and per-sample r ≥ 0.998.
- Figure: ![validation]({{artifact:4007b1e9-9213-4eca-8bb9-8903e7d30c99}})

## 2. Order-invariant vs order-dependent metrics (10,000 simulated orders, seeds 1–10,000)
- **Inhibitory-subset mode Proportion (the pipeline's Inhibitory_Proportion) is order-invariant.** In all 101 real shuffled Inhibitory runs and all 10,000 simulated orders the per-sample value had range 0. Match status is also order-invariant in Full mode.
- **maxrejects effect.** Every feature with a full-coverage ≥99% Inhibitory hit was matched in every order. A few features were matched only through partial-overlap hits, which iddef 2 does not penalise for terminal gaps: meta 12 features (0.01% of reads), bullfrog 3 (0.05%), SERDP 15 (0.24%). Many features whose only exhaustive hits were partial-overlap artefacts were never matched: meta 1,478, bullfrog 16, SERDP 69.
- Inhibitory-mode **richness** varies a little (bullfrog 15% of frogs by ≤1 reference; SERDP 5.6% by ≤2; meta 0.1%). Full richness varies in 30% of bullfrog frogs and 6.7% of SERDP frogs, and in 0% of meta samples.
- **Top-hit-label metrics are strongly order-dependent.** These are the Full-db top hit labelled Inhibitory, and its ratio to Full_Proportion.

| dataset | features with 0<P(inhibitory call)<1 | share of matched reads | median 95% width (top-hit inh. prop.) | % samples width >0.05 / >0.10 | % samples changing quartile (≥5% of orders) | ratio: median width; % >0.10 | % samples whose ">50% inhibitory" class flips (5–95%) |
|---|---|---|---|---|---|---|---|
| meta (n=2,411) | 76 | 28.8% | 0.09 | 68% / 47% | 75% | 0.18; 70% | 21% |
| bullfrog (n=40) | 46 | 33.6% | 0.03 | 18% / 8% | 68% | 0.25; 100% | 13% |
| SERDP (n=372) | 94 | 38.0% | 0.25 | 90% / 84% | 83% | 0.37; 94% | 69% |

- Inhibitory richness among the Full top hits shifts by ≥1 reference within the 95% interval in 97% (meta), 100% (bullfrog) and 100% (SERDP) of samples.
- The distributed file order is not a typical order. In SERDP it sits below the permutation median for most frogs.
- Figure: ![intervals]({{artifact:179fb6aa-0c73-4640-b7af-7bb32d64cf64}})

## 3. Propagation to downstream tests (metric: top-hit-label inhibitory proportion)
- **Meta-analysis contrasts (MWU).** Significance does not depend on order: p<0.05 in ≥99.99% of orders for all three contrasts. Effect size does depend on order:
  - Adult vs larvae: rank-biserial 0.41 in the distributed order, 95% range 0.19–0.49 across orders.
  - Madagascar vs USA: −0.28 (−0.52 to −0.17).
  - *A. boreas* vs *L. catesbeianus*: 0.33 (0.28–0.39).
- The order-invariant Inhibitory-mode metric reverses the sign of the Madagascar–USA contrast (+0.31). It also gives a non-significant boreas–catesbeianus contrast (0.06, p=0.28). So the conclusion depends far more on which metric is used than on order.
- The ratio metric is order-fragile for adult vs larvae: p<0.05 in 35% of orders, and the sign matches the default order in only 39%.
- **SERDP infection (conditional logit, site × season).** OR per SD was 0.86 in the distributed order (p=0.42). Across orders the 95% range was 0.76–1.00 (min 0.74, max 1.09), with p<0.05 in 0% of orders. The ratio metric reached p<0.05 in 3.0% of orders.
- **SERDP load among infected (OLS + stratum).** beta was +0.09 in the distributed order, but the across-order distribution is bimodal: 95% range −0.18 to +0.14, median −0.12, and the sign matches the default order in only 31% of orders. p<0.05 in 0% of orders (ratio metric 0.8%).
  - One 150-bp feature (`6a40d31c25c63fc295917d2a6e3c007b`, 6,187 rarefied reads) is called Inhibitory in 33% of orders and drives the split: mean beta +0.06 when it is called Inhibitory, −0.14 when it is not.
- **Bullfrog pond-mean log Bd (~ metric + Type, 10 ponds).** beta 0.13 in the distributed order (p=0.20), 95% range 0.11–0.17, p<0.05 in 0% of orders. Within-treatment Spearman correlations were never significant.
- Figure: ![tests]({{artifact:97c07914-6cdc-4d36-99a9-6a1d7f4bb069}})

## Interpretation
- A user who takes the label of the reported top hit inherits a hidden, irreproducible random component. vsearch output is fully deterministic for a given file. But any re-ordered, re-downloaded or re-subset FASTA can change up to 38% of matched reads' calls, and per-sample values by up to 0.7 (SERDP) or 1.0 (meta, single-feature samples).
- In these data, permutation uncertainty did not create or remove any of the nominal significance calls for the main metric. It does change effect sizes substantially, and in SERDP it changes the sign of the load association.
- The pipeline's default subset-mode Inhibitory_Proportion avoids order dependence entirely. It does so by counting any inhibitory isolate at ≥99%, which is why it tracks "has any cultured relative" (see earlier rounds).

## Files
- [AmphiBac_order_permutation_results.xlsx]({{artifact:7d523c2a-5932-4770-ae11-3f4fb79ff64e}}) — all tables
- [meta_per_sample_permutation_intervals.csv]({{artifact:78591eb5-a7fe-4994-a03b-8b82d31f3a65}}) — meta per-sample intervals + host metadata
- [bf_per_sample_permutation_intervals.csv]({{artifact:eb9cf079-a0b9-4ff7-8509-64ed5bb8f733}}) — bullfrog per-sample
- [serdp_per_sample_permutation_intervals.csv]({{artifact:cf0d2321-7305-4ea5-aa96-3dc38816cc48}}) — SERDP per-sample
- [meta_per_feature_call_probability.csv]({{artifact:af0f9f00-90b9-4614-b222-17a774a02709}}) — meta per-feature P(inhibitory call)
- [bf_per_feature_call_probability.csv]({{artifact:129fd148-a0be-4cf0-b036-632d177bec48}}) — bullfrog per-feature
- [serdp_per_feature_call_probability.csv]({{artifact:7fb1dd51-b55c-4d71-a001-259070edae39}}) — SERDP per-feature
- [emulator_validation.csv]({{artifact:7a149312-44da-491d-9e47-b3aca73dfa56}}) — emulator vs real runs
- [propagated_tests.csv]({{artifact:551ea260-09ca-4de1-a6ee-cd509f4e7923}}) — test distributions summary
- [tie_mechanics.csv]({{artifact:aa07745f-4169-489c-b59d-91b9a77c1c63}}) — tie mechanics
- [AmphiBac_order_permutation_scripts.tar.gz]({{artifact:a5a4d160-9968-4431-be43-33dcfdb81550}}) — scripts p01–p13
