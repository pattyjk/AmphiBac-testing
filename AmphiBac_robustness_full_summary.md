# Does the AmphiBac workflow predict Bd-inhibitory function? A summary of all tests to date

## Bottom line

The AmphiBac database is internally consistent. The common way it is used, however — closed-reference matching of amplicon features at 99% identity, then counting "inhibitory" reads — does not recover Bd-inhibitory function. The main reasons:

1. **Isolate labels are not shared by identical 16S sequences.** 69–74% of groups of identical fragments contain conflicting Bd labels.
2. **Prediction of a held-out isolate's label is near chance once whole studies are held out.** Cohen's κ is about 0.1 for the pipeline's rules, at every fragment length, region and threshold tested.
3. **Genus alone predicts as well as, or better than, any 16S-identity rule.** This holds for V4 up to a 1,100-bp V2–V8 window.
4. **The per-sample "inhibitory proportion" is indistinguishable from what random labels give.** It measures how much of the community has a cultured relative in the database (ρ = 0.92–0.97 with the shuffled-label value).
5. **No metric tracked Bd infection or load** in either Bd dataset, once site and season (or treatment) were accounted for. A naive analysis of the leopard-frog data does find a significant association, but in the *opposite* direction, and it is driven by site.
6. **The practice is common and inconsistent.** At least 57 of the 116 citing papers whose methods could be read matched their own sequences to the database. Their matching rules range from genus name only to 100% identity.
7. **Paired studies confirm it.** Against each study's own same-population isolates, with that study removed from AmphiBac, the pipeline gives κ ≈ 0.02 (subset mode) to 0.08 (top hit). Genus gives 0.15.

## 1. Workflow and data tested

**The workflow.** We tested the `tools/antifungal_predictions_pipeline.py` script and the older `scripts/Antifungal_Predictions_conda.sh` from github.com/AmphiBac/AmphiBac-Database. Both run:

```
qiime vsearch cluster-features-closed-reference --p-perc-identity 0.99 --p-strand both
```

against one of four subsets: Full (8,739 sequences), Inhibitory (2,831), InhibitoryStrict (2,273) and Facilitating (645). They then report per-sample Total, Proportion and Richness.

q2-vsearch sets no `maxaccepts`, `maxrejects` or coverage options, so vsearch's defaults apply:

- the first k-mer-ranked hit that passes the threshold is accepted;
- terminal gaps are ignored when computing identity;
- there is no minimum query coverage.

We reproduced this exactly with vsearch 2.32.0. QIIME itself was not run.

| Data | Use |
|---|---|
| AmphiBac 2025.1 (qza = FASTA, 8,739 seqs), plus releases 2020.6 and 2023.2 | Isolate-level tests and release stability |
| Amphibian meta-analysis: 2,411 samples, 48,431 deblur sOTUs (90 bp) | Real-data impact and permutation null |
| Bullfrog nitrogen experiment: 40 frogs, 10 ponds, 3,700 ASVs (120 bp), Bd loads | Real-data impact and Bd validation |
| SERDP leopard frogs: 372 frogs, 6 species-sites, 4 seasons, 12,105 sequenced features (150 bp), Bd loads (100 infected) | Large-scale Bd validation |
| In-silico fragments cut from isolates: 90 / 100 / 120 / 150 bp, V4, V3–V4, and ten *E. coli*-coordinate windows of 254–1,100 bp | Length, threshold and region tests |

**Ground truth.** Each isolate's `AntiFungal.Function.with.AP.function` label, which is the column that defines the subsets. Using `Bd.Inhibition` instead gave the same conclusions.

## 2. The database itself

![Database audit]({{artifact:b7464cfe-bbf8-4b68-8705-c7eff66252f4}})

- **Labels are a hard cutoff on a continuous score.** Inhibitory means a CFS score above 80, Facilitating means −20 or below, and everything else is NoEffect. As a result, 790 isolates that the original study called inhibitory are labelled NoEffect.
- **Identical sequences carry different functions.** Among identical 515F fragments with ≥2 labelled isolates, 74% of groups (90 bp), 73% (120 bp) and 69% (253 bp) contain conflicting labels.
- **InhibitoryStrict removes only Facilitating look-alikes.** It drops just the 558 Inhibitory sequences that are 100% identical to a Facilitating isolate, a rule we reproduced exactly. Yet 748 of the 2,273 Strict sequences are 100% identical to a NoEffect isolate.
- **Sequence-format problems:**
  - 5,451 sequences are stored reverse-complemented.
  - 1,296 have no 515F site.
  - 61 placeholder sequences of 2–8 bp can never be hit.
  - 1,355 contain IUPAC codes.
  - `Assay.Method` contains auto-fill errors (from "97 well" through "223 well CFS assay").
- **Coverage along 16S.** Most sequences are Sanger reads spanning *E. coli* positions about 90–870. The largest common window is V3–V4 (296–795, 500 bp), which 6,037 tested isolates (77%, from 27 of 30 studies) span without IUPAC codes. Beyond about 600 bp, coverage collapses and the surviving set becomes biased toward long-read studies (45–51% inhibitory vs 38%).

![Fragment coverage]({{artifact:16e859fe-f259-4d69-8ae2-3505be53950c}})

## 3. How vsearch picks a "top match"

![Top hit vs exhaustive]({{artifact:eaa04451-a001-4e1c-b14b-cf70ffb44660}})

- **The default never picks a lower-identity reference, but ties are the norm.** At the top identity, 61% of matched meta-analysis sOTUs (85% of matched reads) and 77% of bullfrog ASVs have ≥2 tied references. In the longer windows the tie rate is 79–99%.
- **Which tied reference wins depends on file order.** Across 10 shuffles of the reference FASTA:
  - the assigned reference changed for 34% of matched sOTUs (67% of matched reads);
  - the inhibitory yes/no call changed for 14% of sOTUs (26% of matched reads);
  - in the bullfrog data, the call changed for 18% of ASVs;
  - in longer windows, 7–39% of calls changed.

  Thread count changed nothing.
- **Mixed labels are common among qualifying references.** 40% of matched sOTUs (80% of matched reads) have both inhibitory and non-inhibitory references at ≥99%. The default label disagrees with the majority label for 104 of 553 meta-analysis hits.
- **Partial-coverage references pass at "100%".** Because terminal gaps are ignored, a reference covering only 20–40 bp of a 90-bp query passes the 99% threshold. In real data, 22 of 553 default hits are partial. In windows of 400 bp and longer, 36–77% of qualifying hits are partial.
- **Strand matters.** With q2-vsearch's own default (`--strand plus`), 294 of 553 meta-analysis hits are lost. The AmphiBac scripts do set `both`.
- **Ambiguity codes.** N and compatible IUPAC codes score as matches.
- **Missed matches.** The default `maxrejects 32` almost never misses a genuine full-coverage match (5 of 537).

![Edge cases]({{artifact:874b4b84-6fd0-48c4-818e-d07171d442a0}})

## 4. Thresholds and fragment length

- **The same percentage means different mismatches at different lengths.** 99% identity allows 0 / 1 / 1 / 1 / 2 / 4 mismatches at 90 / 100 / 120 / 150 / 253 / 429 bp, and up to 11 at 1,100 bp. **At 90 bp, 99% is identical to 100%**: every per-sample output was the same.
- **Leave-one-out prediction at 120 bp and 99%, holding out only the isolate itself:**
  - Inhibitory-subset mode: sensitivity 0.94, **specificity 0.15**, PPV 0.40 (prevalence 0.38), κ = 0.07.
  - Default top-hit label: κ = 0.15.
  - Majority vote of all hits: κ = 0.32.
- **Holding out whole studies, every rule falls to κ ≈ 0.09–0.12.** The apparent skill largely comes from isolates in the same study sharing labels.
- **Longer windows: ten regions on a fixed panel of 1,397 isolates, leave-study-out at 99%:**
  - default top hit κ 0.13–0.23;
  - subset rules 0.15–0.31;
  - majority of full-coverage hits 0.27–0.44.

  No region or length beats the **genus-only baseline (κ 0.45)**.

![Kappa curves]({{artifact:ad797c86-021b-4658-9095-64c039f08558}})

![Window synthesis]({{artifact:acecbd8d-553e-4fa7-9bc4-9cbb3e50d459}})

## 5. Does 16S similarity carry functional signal at all?

- **Genus baseline.** Predicting each held-out isolate from its genus's majority label (RDP v18 SINTAX) gives leave-study-out κ 0.21 (120 bp) and 0.25 (V4), against ≤0.12 for every 99% V4 rule. Genera are themselves mixed: 46 of 50 large genera are 10–90% inhibitory, and *Pseudomonas* ranges from 0% to 100% inhibitory between studies.
- **Agreement against sequence difference.**
  - Identical-V4 pairs agree mostly when they come from the same study: κ 0.30 within studies vs 0.09 between.
  - Between studies, P(partner inhibitory | focal inhibitory) is a flat 0.52–0.57 from 0 to 10 V4 differences, then drops to random.
  - On the fixed panel it is 0.77–0.83 at ≥99.5% identity in every window (random pairs 0.51). That is a genus-scale step, not a fine-scale decay.
- **Correction to an earlier result.** Round 2 found that near-full-length identity seemed to help. Round 3 showed this was an isolate-set effect: on the fixed panel, V4 already reaches the same agreement.
- **Phylogenetic signal is weak and the same in every window.** Fritz–Purvis D is 0.52–0.76 (typical conserved microbial traits average about 0.07). τD is not deeper than expected once labels are shuffled within studies. Blomberg's K on the score is about 3×10⁻⁵.

![Decay curves]({{artifact:3b809a47-3603-487e-a1a0-2995d18e7e88}})

![Genus baseline]({{artifact:b345b5be-def9-4292-912e-638073988986}})

## 6. Database releases and out-of-sample prediction

- **Stable outputs.** Full and Facilitating outputs are stable across releases (ρ ≥ 0.999), and Inhibitory output is nearly stable from 2023.2 to 2025.1 (ρ = 0.995).
- **InhibitoryStrict is not stable from 2020.6 to 2025.1.** 26% of samples shift by more than 0.05 because the subset was redefined; only 1 of 7,382 shared isolates was relabelled. The change flips the *A. boreas* vs *L. catesbeianus* contrast (2020.6: r = −0.20, p = 0.0004; later releases: n.s.), and essentially one reference drives the flip.
- **Top-hit labels drift between releases without new data.** Between 2020.6 and 2023.2, 11% of matched sOTUs changed their top-hit label, mostly through tie-breaking rather than new isolates.
- **Prospective prediction is at chance.** Predicting isolates added after 2020.6 and after 2023.2 from the older release gives κ ≤ 0.15 at 120 bp and V4, and −0.23 to 0.22 in the longer windows. These estimates rest on few studies.

## 7. Impact on real community data

![What inhibitory reads are]({{artifact:a6185c00-32ba-40a4-8f96-840217fe9300}})

- **Many "inhibitory" reads come from mostly non-inhibitory lineages.** In the meta-analysis, 56% of pipeline-"inhibitory" reads come from features whose ≥99%-identical isolates are mostly non-inhibitory, and 78% come from features that also hit the Facilitating subset. In the SERDP frogs, the median frog has 74% of reads "inhibitory" and 57% "facilitating".
- **The metric tracks cultured coverage.** Inhibitory_Proportion vs Full_Proportion: ρ = 0.947 in the meta-analysis and r = 0.996 in the bullfrogs.
- **Label-permutation null.** We shuffled the labels among database isolates 1,000 times while keeping every sequence match fixed. The observed metric correlates with the shuffled-label mean at ρ = 0.92–0.97. Only 1.9% of meta-analysis samples exceed the null 95th percentile, below the 5% expected by chance; for bullfrogs and SERDP frogs, none do.
- **Host comparisons appear with random labels too.** The adult vs larvae and Madagascar vs USA contrasts are significant in 100% of shuffles.
- **Method choices change conclusions.** *A. boreas* vs *L. catesbeianus* is n.s. under the pipeline but significant under consensus rules. Madagascar vs USA reverses direction under a ≥75%-of-isolates rule.
- **Minor sensitivities.** InhibitoryStrict vs Inhibitory barely matters (ρ = 0.996). Rarefaction seed and thread count are negligible.

![Permutation null]({{artifact:c97dafc4-3662-4743-acac-6b61fe7b2a71}})

![Contrasts by method]({{artifact:895214e2-5fa4-46c5-93ed-fb711a3e6ed2}})

## 8. Validation against Bd infection

**Bullfrog nitrogen experiment** (40 frogs, 10 ponds)

- Bd load is completely separated by treatment, so only within-treatment or pond-level tests are interpretable.
- No AmphiBac metric was associated with log Bd in any analysis: naive, within-treatment, mixed model with pond, or pond means with treatment (15 metrics; minimum BH q = 0.29).
- All 15 naive correlations are positive.
- Power is low: with 10 ponds, only |r| ≥ 0.80 is detectable.

![Bd bullfrog]({{artifact:a3094ba4-400d-4264-9553-942866638a1b}})

**SERDP leopard frogs** (372 frogs; 100 infected)

- **Naive test.** Infected frogs have *more* "inhibitory" reads (median 0.81 vs 0.70; MWU p = 0.002, q = 0.012).
- **Within site × season** (conditional logit; 344 frogs, 14 strata), no metric predicts infection. For pipeline Inhibitory at 99%, the OR is 1.02 per SD (95% CI 0.72–1.45). The minimum BH q is 0.79.
- **Bd load among infected frogs.** No metric is associated (all p ≥ 0.18).
- **Adjusting for cultured share** changes nothing.
- **Label-permutation null.** The real labels' effect is typical of shuffled labels (71% of shuffles are as or more extreme).
- **Why the naive test is significant.** Sites differ in both prevalence and "inhibitory" share: TN, LA and Panama are high on both, NM and PA low.

![SERDP Bd]({{artifact:a609d8f1-5a73-4b9e-bd23-412926a87c8a}})

## 8b. Is the score distribution multimodal?

Hartigan's dip test on the 6,245 CFS-scored isolates, with 2,000 bootstraps and a ±0.05 jitter to break ties:

- **Pooled scores are unimodal** (D = 0.0034, p = 0.96). The distribution is left-skewed, with a main mode near 93 and a shallow shoulder near 14; the trough between them is only 4% lower than the shoulder.
- **Within studies**, two of 14 are bimodal after BH correction: AEM.01486-15 (n = 719, modes near 6 and 87) and AEM.00186-17 (n = 119, modes near 15 and 96).
- **Within genera**, none of 29 is bimodal after correction.

The >80 / ≤−20 label cutoffs therefore split one skewed distribution rather than separating natural groups.

![Dip test]({{artifact:d5a3e1f4-a7d2-43cf-8c45-2a66adb33bdb}})

Tables: [score_diptest.xlsx]({{artifact:f2384e74-ce47-46c6-80ea-a687bbf6e2fd}}), [diptest_by_group.csv]({{artifact:fb33c488-299a-481f-8221-c62beeee46dd}}).

## 8c. Tree-based prediction (PICRUSt2-style placement)

We placed reads into a 7,279-isolate reference tree with EPA-ng and predicted traits with castor hidden-state methods, Brownian/OU score models and a study-aware phylogenetic Gaussian process.

- **Benchmark (30-fold leave-study-out).** No tree-based method beat genus majority or majority voting over ≥99% hits:
  - nearest tip and castor methods: V4 κ 0.11–0.12, AUC 0.57;
  - study-aware phylogenetic model (the only one to tie the baselines): AUC 0.63, against 0.64 for genus majority and 0.62 for majority voting.
- **Variance in the inhibition score:** phylogeny 0.26, study 0.25, Bd genotype 0.08, residual 0.42.
- **Leopard-frog Bd test.** No tree-based metric predicts infection within site × season (OR 0.89–1.27 per SD, all CIs include 1).
- **Shuffled labels.** The observed effects sit inside the null.

![Placement benchmark]({{artifact:b95974d0-c4fd-4253-8aee-a4be3b20e4ef}})

Report: [AmphiBac_placement_method_report.md]({{artifact:2ae9d06a-ec9d-4758-8c63-e4b657649c3b}}).

## 8d. Primer use in amphibian microbiome studies

Of 180 primary studies that state their region:

| Region | Share of studies |
|---|---|
| V4 | 57% |
| V3–V4 | 29% |
| V4–V5 | 5% |
| Others | 9% |

V4 fell from 80% of studies in 2015–17 to 42% in 2024–26. Gut studies favour V3–V4, mostly with 338F/806R.

Report: [amphibian_16S_primer_survey.md]({{artifact:1b572d29-5d04-4b2b-bcf0-310054fbd134}}).

## 8e. Host-matched databases and reference-order uncertainty

**Host- or location-matched databases.** Restricting AmphiBac to the study host or location makes predictions worse:
- **Coverage:** there are no Massachusetts bullfrog isolates, and no *R. pipiens* or *R. blairi* isolates. Restriction leaves 55–98% of previously matched reads unmatched in 9 of 11 sets.
- **Partial overlaps:** the remaining matches are often partial overlaps. Every Tennessee *L. sphenocephalus* match covers only part of the 150-bp read.
- **Predictability:** same-host references predict isolate labels worse than all references (κ 0.08 vs 0.28 on the same 551 isolates).
- **Bd:** no restricted metric predicts Bd.

**Reference-order permutations.** We used 10,000 orders, with an emulator of vsearch's k-mer, length and position tie-break validated at 100% agreement against more than 1,200 real shuffled runs.
- **Pipeline metric:** subset-mode Inhibitory_Proportion is exactly order-invariant.
- **Top-hit metrics:** the median 95% width of the top-hit inhibitory share is 0.09 (meta), 0.03 (bullfrog) and 0.25 (SERDP).
- **Conclusions:** no significance call flipped, but effect sizes shifted widely, e.g. adult vs larvae rank-biserial 0.19–0.49.

![Per-sample intervals]({{artifact:179fb6aa-0c73-4640-b7af-7bb32d64cf64}})

Report: [AmphiBac_hostmatch_order_addendum.md]({{artifact:678a8fea-f010-4664-ae66-5bf4d79e4b28}}).
## 8f. Worked example: one tied SERDP read through the order permutations

This traces one feature through 1,001 **real** vsearch runs (the distributed order plus 1,000 shuffles, q2-vsearch closed-reference settings): `6a40d31c…`, 150 bp, 6,187 rarefied reads. It is the feature that drives the bimodal SERDP Bd-load slope in 8e.
- **The tie.** The read matches 13 isolates at ≥99%. Six match at 100%, all *L. catesbeianus* isolates from one Virginia study (AEM.01486-15): three Inhibitory (scores 93, 85, 83) and three NoEffect (77, 17, 11). The isolate scoring 77 is NoEffect only because of the 80 cutoff. Seven more match at 99.3% (5 NoEffect, 2 Inhibitory).
- **Which isolate wins.** All six 100% matches share the same number of 8-mers with the read (143), so shorter targets are preferred. Only three ever won: AmphiBac_5440 (Inhibitory, 347/1,001), 5442 (NoEffect, 330) and 5563 (NoEffect, 324). The two 1,270-bp isolates are the shortest overall; 5563 is the shortest of the forward-stored ones. The longer Inhibitory isolates 5506 and 5507 never won.
- **Call rates.** The distributed file order picks 5440, so the default run calls the read Inhibitory. Across orders that happens 35% of the time (emulator: 33%). The pipeline's subset mode calls it Inhibitory 100% of the time.
- **Per frog.** 12 frogs carry the read: 10 infected Panama (Ngäbe-Buglé) frogs, where it makes up to 0.47 of reads, plus two others at ≤0.006. In the three heaviest carriers, the top-hit inhibitory share falls from 0.56–0.61 to 0.09–0.18 when the call switches to NoEffect.
- **Bd load.** Flipping only this call, with every other feature at its default-order call, moves the Bd-load slope among the 100 infected frogs from +0.09 per SD (p = 0.49) to −0.14 (p = 0.30). That single switch spans most of the 10,000-order 95% range (−0.18 to +0.14).

![Tied feature walkthrough]({{artifact:cc55e8c5-2a1e-40fb-a705-c3b233db73fd}})

Tables: [tied_feature_references.csv]({{artifact:4690de67-7760-4a8e-8c68-5c9242c8532a}}), [tied_feature_real_runs.csv]({{artifact:db06260b-4cca-4864-be1b-6313824e67f8}}), [tied_feature_frogs.csv]({{artifact:0454ee67-98ad-468e-8624-dffc9d694200}}). Script: [tied_feature_walkthrough.py]({{artifact:87cffd2a-7b89-4540-a616-23184683316b}}).

## 8g. How often the database is used to predict function

The survey covers works citing Woodhams et al. 2015 (*Ecology* 96:595, doi 10.1890/14-1837.1). OpenAlex lists 222 citing works. After removing 12 PeerJ peer-review reports and 6 preprint/article duplicates, 204 unique papers remain.
- **Text sources.** Full text came from Europe PMC, Unpaywall and publisher open access. Each paper was classified from its methods and citing passages by a language model, with a second pass on methods sections. The 25 weakest calls were re-run on the whole text. Texts under 10k characters were treated as abstracts.

| Use of the database | Methods read | Abstract only | Total |
|---|---|---|---|
| Predictive: matched own amplicon OTUs/ASVs | 54 | 20 | 74 |
| Predictive: matched own isolate 16S | 3 | 0 | 3 |
| Records used as data (meta-analysis, assays, probiotic choice) | 6 | 8 | 14 |
| Background citation only | 51 | 8 | 59 |
| Cannot tell from the abstract | 2 | 52 | 54 |
| **Total** | 116 | 88 | 204 |

- **Confirmed predictive use: 57 papers**, about half of the 116 whose methods could be read. A further 20 are probable from the abstract alone. The 52 unresolved papers are mostly closed access (Wiley, Springer, *Proc. B*), so the true count is likely 75–90.
- **Growth.** Confirmed predictive papers numbered 1–3 per year in 2016–2019, 5–7 per year in 2020–2022, and 7–11 per year in 2023–2025.
- **No common rule.** Of the 57:
  - 17 used 100% identity, 8 used 99% and 8 used 97–98%;
  - 2 used a lower threshold or only a BLAST e-value;
  - 5 matched by genus or family name without comparing sequences, and 2 used phylogenetic methods;
  - 15 did not state the rule or cited another paper for it.
  
  Only 11 used vsearch, as the current pipeline does. Given sections 3–4, "proportion inhibitory" values from these studies are not directly comparable.
- **Not verified.** The table also records whether each study linked predictions to Bd and whether it ran its own assays. Neither field was audited.
- **Coverage.** Studies that cite only the AmphiBac 2023/2025 releases or the GitHub repository are not counted.

Table: [woodhams2015_citers_classified.csv]({{artifact:f12dbce2-573d-40ae-a5e7-125f85c95d8e}}). Script: [woodhams2015_citation_survey.py]({{artifact:456d6a28-c6d6-40d7-bfa6-b92b6f30a555}}).

## 8h. Does a different Inhibitory cutoff (or an IC50) help?

**No IC50 is available.** AmphiBac assays are single-concentration CFS tests, and score = (1 − proportional Bd growth) × 100. A cutoff of 50 is therefore the single-dose analogue of an IC50.

**Test.** The Inhibitory cutoff was swept from 0 to 95 and every isolate was relabelled at each cutoff.
- **Isolates:** about 2,900–3,100 CFS-scored isolates from 17 studies; agar-plate isolates have no score.
- **Prediction:** from exhaustive ≥99% full-coverage hits to isolates in *other* studies, at 120, 150 and 254 bp of V4.
- **Comparisons:** pipeline-style rules, the mean score of hits, and the genus baseline. Intervals come from 300 study-cluster bootstraps.

**Results:**
- **No cutoff gets κ above about 0.17**, for any rule or length.
- **50 is no better than 80.** At 150 bp: any inhibitory hit 0.13 vs 0.10; majority of hits 0.09 vs 0.15; genus −0.02 vs 0.15. Intervals overlap across cutoffs.
- **Predicting the continuous score is also weak:** Spearman ρ = 0.22–0.25, whether from the hits' mean score or from genus.
- **The cutoff mainly changes prevalence:** the share of isolates called Inhibitory falls from 91% at 0 to 19% at 95.
- **Within-study agreement is roughly double.** Restricting hits to the *same* study gives κ 0.25–0.43 and ρ ≈ 0.48. Some of these hits may be re-isolates of the same strain, so this is an approximate ceiling. The limit is lab/assay heterogeneity, not the threshold.

![Cutoff sweep]({{artifact:92fb2e8a-8a52-4572-a30f-1883a1707dd4}})

Table: [cutoff_sweep_results.xlsx]({{artifact:968698c3-29c6-497b-9623-3900eb29db8f}}). Script: [cutoff_sweep.py]({{artifact:6a9688d3-51ad-44d0-808e-93139824f7ca}}).

## 8i. Paired studies: AmphiBac vs each study's own isolates

**Question.** In studies that both cultured and sequenced, the study's own Bd-assayed isolates give a label for some of its own ASVs. Does the pipeline's call for those ASVs agree once the study's isolates are removed from AmphiBac?

**Data.** Eight curated paired studies were considered ([paired_culture_amplicon_studies.xlsx]({{artifact:fbc12368-3f52-47fb-a4f7-857d92a4a5fe}})). Four were usable at sequence level:
- **Walke 2015/2017:** 719 scored isolates. The reads are one pooled file with no sample IDs. The supplied mapping file belongs to the separate bullfrog Bd experiment, so this study is pooled.
- **Rebollar 2019:** 354 scored isolates. Reads were demultiplexed by index read into 61 barcode groups, R01–R61, covering 56% of spots; there are no host labels.
- **Bletz 2018:** 694 isolates with Bd and Bsal assays from Dryad, from the same populations but different swabs.
- **Ellison 2019 (*Mol Ecol*):** 43 isolates from the same swabs. Scores were derived from the raw OD data and are provisional; per-sample Bd is available.

The other studies were not usable: Jiménez 2019 has no isolate sequences, the *PNAS* 2019 isolates have no antifungal assays, and Assis 2020 gives only genus-level calls.

**Methods.**
- **ASVs:** a uniform 150-nt V4 window after 515F, denoised with vsearch UNOISE3 and uchime3.
- **Own-isolate label:** the majority label (score >80) of the study's isolates covering the whole ASV at 100% identity.
- **Pipeline:** q2-vsearch emulation against AmphiBac 2025.1 with that study's records removed (736, 384, 694 and 0 records respectively). Ellison is not in AmphiBac.

**Results (ASV level, Bd).**

| Rule (study held out) | Walke | Rebollar | Pooled 4 studies |
|---|---|---|---|
| Subset mode (pipeline) | κ 0.02 (−0.09–0.13) | 0.04 (−0.07–0.14) | 0.02 (−0.02–0.07) |
| Top hit | −0.03 | 0.30 (0.06–0.52) | 0.08 |
| Majority of ≥99% hits | 0.06 | 0.10 | 0.04 |
| Genus majority | 0.13 | 0.13 | 0.15 (0.03–0.27) |
| *Subset, own isolates included* | *0.24* | *0.18* | *0.15* |
| *Top hit, own isolates included* | *0.35* | *0.41* | *0.32* |

- **Held out, the pipeline barely agrees with the study's own isolates.** Subset mode catches most own-Inhibitory ASVs (sensitivity 0.85 pooled) but calls most non-inhibitors Inhibitory too (specificity 0.20, PPV 0.20).
- **Including the study's own isolates roughly doubles κ** (pooled top hit 0.08 → 0.32). Much of the apparent agreement when a lab analyses its own data with AmphiBac is circular: the ASVs match the lab's own deposited isolates.
- **Genus does at least as well as any sequence rule**, matching the cross-study results in sections 6 and 8h.
- **Bletz 2018 is uninformative for Bd.** Only 3 of 694 isolates score >80, and no own-matched ASV is Bd-Inhibitory. For Bsal, held-out κ is −0.01 to 0.11, but only about 3 own-Inhibitory ASVs exist.
- **Ellison 2019 is too small to judge:** 13 ASVs, with κ ranging from −0.38 to 0.37 across rules.

**Coverage.** A study's own isolates cover a minority of its reads: 37% (Walke, pooled), 20% (Rebollar), 22% (Bletz) and 3% (Ellison). AmphiBac ≥99% hits from other studies cover 51–64%. So the pipeline mostly labels reads using isolates from other labs and other hosts.

**Per sample.** The pipeline's inhibitory proportion is unrelated to the inhibitory proportion from the study's own isolates:
- Rebollar: ρ = 0.00 (60 groups).
- Ellison: ρ = −0.12 (98 samples).
- Bletz, against the Bsal-based own proportion: ρ = −0.07.

**Bd (Ellison 2019, 87 frogs, 56 infected).**
- **Pipeline metric:** not related to infection status (OR per SD 0.97 after site and life stage, p = 0.93). Its slope on load among infected frogs is +0.34 per SD (p = 0.06).
- **Own-isolate metric:** inhibitory proportion and mean score were associated with *higher* load before adjustment (p = 0.005 and 0.03). Neither association remained after site and life stage (p = 0.72 and 0.97).
- **Bletz 2018:** in Experiment 1 (6 exposed vs 6 control, Mann–Whitney), the Bd-based metrics did not differ (p ≥ 0.70). The own-isolate Bsal-inhibitory proportion after exposure was nominally higher in exposed animals (median 0.051 vs 0.022, p = 0.065), but its change from before exposure did not differ (p = 0.82). With n = 12 this is weak evidence.

**Interpretation.**
- Even in the best test case (same animals or populations, same lab assay), AmphiBac calls for an ASV do not reproduce the study's own measurement once that study is held out.
- Agreement in published analyses is inflated when a study's isolates are in the database.

**Caveats.**
- ASVs are treated as independent in the bootstrap intervals.
- The Rebollar groups may include non-study samples, and Walke has no per-sample data.
- Ellison scores are provisional.
- All results use 150-nt V4; Ellison's V3–V4 reads were reduced to merged reads only (28%).

![Paired-study κ by rule]({{artifact:67003799-5d4f-46c0-b682-76227461a438}})

Files: [paired_asv_agreement.xlsx]({{artifact:8385f84a-000a-4ccf-8e41-2e23da109505}}) (per-ASV calls, agreement, per-sample, Bd tests, Rebollar barcodes, README), [paired_asv_tables.tar.gz]({{artifact:844e1d8c-bce4-4e57-9ed8-0f2f2f09ec2c}}), [paired_analysis_scripts.tar.gz]({{artifact:7d892d66-f3c3-47bc-bc7e-b06db650f9c2}}).

Report: [AmphiBac_paired_study_report.md]({{artifact:d9143200-974c-476e-b83d-fd1ba83f0be6}}).

## 9. Recommendations

1. **Report matches as mismatches over the read length, not percent identity.** At 90 bp, 99% means an exact match.
2. **Summarise all qualifying hits instead of the first one.** Use vsearch directly with `--maxaccepts 0 --maxrejects 0`, and report n, % inhibitory, % facilitating and median score per feature.
3. **Require full query coverage** (`--query_cov 1.0` or alnlen = qlen). This matters more as fragments get longer.
4. **Always use `--strand both`.** Most reference sequences are stored reverse-complemented.
5. **Do not interpret subset-mode "inhibitory proportion" as function.** Report the cultured share (Full_Proportion) and the Inhibitory/Full ratio alongside it, and test any claim against a label-permutation null run within genus and within study.
6. **Benchmark against a genus-only baseline**, validating with whole studies held out. A 16S rule that cannot beat genus adds nothing.
7. **Control for site and season (or treatment) in any Bd analysis.** Naive correlations are confounded.
8. **State the database release and subset.** InhibitoryStrict results from 2020.6 are not comparable with later releases.
9. **Use the continuous CFS score where possible** rather than the >80 / ≤−20 cutoffs.
10. **Exclude** NotTested/PossibleContamination isolates, references under 100 bp, and references that do not span the amplicon.
11. **Suggested database fixes:**
    - orient all sequences;
    - remove the placeholder sequences;
    - build InhibitoryStrict against NoEffect isolates as well;
    - ship consensus tables with scores;
    - fix `Assay.Method`.
12. **Collapse or report ties.** When ≥99% hits disagree, call the feature ambiguous, or report the range of results across shuffled reference orders (see 8f).
13. **State the matching rule exactly:** tool, identity, coverage, strand, database release and subset. 15 of 57 predictive papers did not.
14. **Don't tune the 80 cutoff.** Model the continuous score with study (assay) as a random effect, or standardise scores within study before pooling. New assays with dilution series would allow true IC50s.
15. **Hold out the study's own isolates when validating.** In paired studies, including them roughly doubles κ (0.08 → 0.32 pooled, top hit), so agreement within a lab's own data is partly circular.

## 10. Caveats

- **Ground truth.** The labels are in-vitro assays from many labs and Bd strains, so part of the disagreement is assay noise rather than biology.
- **Study-level dependence.** Isolates and pairs within a study are not independent, so we used study-cluster bootstraps. Prospective and some window cells rest on 1–5 studies.
- **Genus baseline.** It uses RDP v18 SINTAX on full isolate sequences. It is optimistic relative to calling genus from an amplicon fragment.
- **Pipeline emulation.** QIIME was emulated with the identical vsearch call. Subset mode under permutation was approximated by "any ≥99% full-coverage hit is Inhibitory", which matches the pipeline for 276 of 291 inhibitory features.
- **Long fragments.** Only isolate-level tests could be done above 250 bp; all community datasets are 90–150 bp.
- **Bullfrog Bd test.** It has low power (10 ponds).
- **SERDP rep-set.** It was taxonomy-filtered before alignment: 3,490 features (about 11% of reads, up to 80% in some Vermont frogs) have no sequence and were treated as unmatched. Excluding the worst-affected frogs did not change the results.
- **Host contrasts.** These ignore study and site confounding. They show method sensitivity, not biology.

## 11. Open items

- **Bullfrog map (`nut_enrich_16S_map.txt`).** `NitrateuM` appears to be sorted independently of the rows and should be checked against the original. Six mapped frogs lack 16S data (P10F1, P2F1, P2F2, P2F4, P2F5, P8F5).
- **SERDP.** Rerun with the unfiltered rep-seqs if they are available.
- **Tests not yet run:**
  - leave-host-out or leave-country-out validation;
  - assay heterogeneity (Bd genotype, temperature, agar vs CFS);
  - upstream processing (primers, trimming, DADA2 vs deblur, 97% OTUs);
  - multi-copy 16S and Sanger error;
  - duplicate `Isolate.ID` entries;
  - splitting Bd from Bsal inhibition.
- **Citation survey.** 52 citing papers could be judged only from the abstract. Their full texts (via institutional access) would close the 57–90 range.

## 12. File index

**Reports**
- [Round 1: database audit, top-hit mechanics, length × threshold, real-data impact]({{artifact:68283292-5522-4f52-929b-09571edb7201}})
- [Round 2: permutation null, genus baseline, phylogenetic signal, releases, bullfrog Bd]({{artifact:8887dd6f-7443-4aef-bcb0-4245efda3b44}})
- [Round 3: ten 16S windows of 254–1,100 bp]({{artifact:20ea9098-1f41-46ab-8a4c-2f31e52fd496}})

**Key tables**
- Database audit: [audit_per_sequence.csv]({{artifact:76d77534-2987-4b44-99f0-acd51045b8ce}}), [audit_identical_V4_groups.csv]({{artifact:5e4d146f-ca89-4e3c-b0f7-36f9c727bfa9}})
- Top-hit mechanics: [tophit_summary.xlsx]({{artifact:4059a0e7-156e-4368-8fdd-b4b6c5ef185b}})
- Length × threshold: [lt_length_threshold_summary.xlsx]({{artifact:dca8f85d-9cb7-452e-a77c-26d9476052e0}}), [lt_mismatch_table.csv]({{artifact:12f3afda-f4b1-4728-ab84-a002a2f21761}})
- Real-data impact: [realdata_impact_summary.xlsx]({{artifact:c0bc7610-92fa-444e-9876-a604fe00dd95}}), [meta_per_sample_methods.csv]({{artifact:c99f3b59-0fb2-40d3-9194-db86479dee0c}}), [bf_per_sample_methods.csv]({{artifact:ff449495-bf24-4712-bb39-89fc00548481}})
- Null and genus: [round2_null_genus_summary.xlsx]({{artifact:425f654c-4028-4320-839e-98cb6910588d}}), [genus_assignments.csv]({{artifact:5f91e1d1-4bab-4e34-903c-fda11446522b}})
- Phylogenetic signal: [r2_phylo_signal_summary.xlsx]({{artifact:9ba920fb-3c5b-45ac-82ee-3582c75f5193}})
- Release stability: [AmphiBac_round2_release_stability_summary.xlsx]({{artifact:2fcd0bbb-7230-4184-bb96-6f6470fe1b1c}})
- Bullfrog Bd: [bd_validation_tests.xlsx]({{artifact:4cd17f3e-ee0e-486a-9346-e704ebf9c55d}})
- SERDP Bd: [SERDP_AmphiBac_Bd_tests.xlsx]({{artifact:c994ebbb-e60f-4167-9653-d84df84ef997}}), [serdp_per_sample_metrics.csv]({{artifact:2a986529-db6d-4fcf-8735-feae89a6747d}})
- Fragment coverage: [candidate_windows.csv]({{artifact:7c5c57f9-3120-4954-aeb4-21cd6e9af97f}})
- Round-3 windows: [r3_window_synthesis_panel.csv]({{artifact:e4665aa1-6e8a-4540-8839-8f34cdb3212b}}), [AmphiBac_round3_prediction_grid.xlsx]({{artifact:921588d3-1d6b-46a5-93b9-723bd664c464}}), [r3_decay_signal_tables.xlsx]({{artifact:0535790a-e640-43d0-a091-559c1079883f}}), [AmphiBac_round3_prospective_windows.xlsx]({{artifact:fbb2fab4-9ad3-4c89-a2da-afd9839f66a3}})

**Inputs and scripts**
- Inputs: [database bundle]({{artifact:635e3d25-3992-409c-9c87-382965c05cb1}}), [window fragments]({{artifact:305d0f1b-add2-4c0a-b447-f74675b90744}})
- Round 1: [audit, top-hit, length/threshold, real-data]({{artifact:33a0a17c-c8e2-46d1-9ed3-15d1899b5395}})
- Round 2: [null/genus]({{artifact:eeb11567-e054-4544-8bcf-41fe7a372658}}), [phylogenetic signal]({{artifact:90f4dade-a5dc-4dca-860c-a908783232f5}}), [releases]({{artifact:9bf9b060-cf92-4381-9fe8-84e588464d41}}), [bullfrog Bd]({{artifact:24d6f2d9-53cc-4dca-9b29-3e8b2ea48061}})
- SERDP: [SERDP]({{artifact:f07dba6f-e7ba-4831-9e23-1a4be319829c}})
- Round 3: [prediction grid]({{artifact:603fd98b-a8be-422d-b936-0de50d9d70fc}}), [decay/signal]({{artifact:2e7f4b77-7a9f-42f6-87a2-28aa43d4531e}}), [prospective]({{artifact:d0c3d1ff-d25b-47bf-9930-e3ea56bd17a6}})

**Worked example and citation survey (8f, 8g)**
- [fig_tied_feature_walkthrough.png]({{artifact:cc55e8c5-2a1e-40fb-a705-c3b233db73fd}}), [tied_feature_references.csv]({{artifact:4690de67-7760-4a8e-8c68-5c9242c8532a}}), [tied_feature_real_runs.csv]({{artifact:db06260b-4cca-4864-be1b-6313824e67f8}}), [tied_feature_frogs.csv]({{artifact:0454ee67-98ad-468e-8624-dffc9d694200}}), [tied_feature_walkthrough.py]({{artifact:87cffd2a-7b89-4540-a616-23184683316b}})
- [woodhams2015_citers_classified.csv]({{artifact:f12dbce2-573d-40ae-a5e7-125f85c95d8e}}), [woodhams2015_citation_survey.py]({{artifact:456d6a28-c6d6-40d7-bfa6-b92b6f30a555}})

**Cutoff sweep (8h)**
- [fig_cutoff_sweep.png]({{artifact:92fb2e8a-8a52-4572-a30f-1883a1707dd4}}), [cutoff_sweep_results.xlsx]({{artifact:968698c3-29c6-497b-9623-3900eb29db8f}}), [cutoff_sweep.py]({{artifact:6a9688d3-51ad-44d0-808e-93139824f7ca}})

**Paired studies (8i)**
- [AmphiBac_paired_study_report.md]({{artifact:d9143200-974c-476e-b83d-fd1ba83f0be6}}), [paired_asv_agreement.xlsx]({{artifact:8385f84a-000a-4ccf-8e41-2e23da109505}}), [paired_kappa_by_rule.png]({{artifact:67003799-5d4f-46c0-b682-76227461a438}}), [paired_asv_tables.tar.gz]({{artifact:844e1d8c-bce4-4e57-9ed8-0f2f2f09ec2c}}), [paired_analysis_scripts.tar.gz]({{artifact:7d892d66-f3c3-47bc-bc7e-b06db650f9c2}}), [paired_culture_amplicon_studies.xlsx]({{artifact:fbc12368-3f52-47fb-a4f7-857d92a4a5fe}})
