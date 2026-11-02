# AmphiBac addendum: do host- or location-matched reference databases help?

**Question.** Users often restrict AmphiBac to isolates from their own host species or region before running the q2-vsearch closed-reference step. Does this restriction give metrics that better predict Bd outcomes? Does a same-host reference predict an isolate's function better than the whole database does?

**Short answer.** No. Restriction mostly removes matches. On the panel where same-host prediction is possible at all, it predicts isolate function worse than the full database. It also introduces an alignment artefact: in several restricted sets, a large share of the "matches" (all of them for Tennessee) are partial overlaps that vsearch scores as 99-100% identity because terminal gaps are ignored.

## Methods (brief)
- **Pipeline emulation.** Every restricted Full / Inhibitory / InhibitoryStrict / Facilitating subset was searched with the exact q2-vsearch call: `vsearch --usearch_global --id 0.99 --strand both --qmask none --minseqlength 1` with default maxaccepts 1 and maxrejects 32 (vsearch 2.32.0, 2 threads).
  - Per sample: Proportion = matched reads / rarefying depth; Richness = number of distinct reference IDs matched.
  - Full subsets keep tested and untested isolates. Inhibitory and Facilitating membership was taken from the released FASTAs, which are identical to the Func labels (checked).
- **Meta-analysis.** All 169 meta sets were run as real restricted vsearch searches (83 same-host, 82 same-host+country, 4 same-country). Each search used the features present in that set's samples, so the results are exact and no filtering emulation was needed.
  - The exhaustive Full-DB search (maxaccepts 0, maxrejects 0, id ≥ 0.99) was filtered to each set's references to give a full-query-coverage "any-hit" Inhibitory proportion. This was compared with the real runs; see the confirmation section below.
- **Permutation intervals.** Func labels were shuffled 200 times among the tested isolates of each restricted set. The per-sample null mean and 2.5–97.5% interval of the any-hit Inhibitory proportion are in the per-sample tables.
- **Bd tests.**
  - SERDP: conditional logit of infection within site × season strata, and OLS of log10 Bd + stratum among infected frogs. Fitted pooled over the sites where each set applies, and per site, with the same frogs scored with the full database for comparison.
  - Bullfrog: within-treatment Spearman correlation, and pond-mean OLS with Type as a covariate.
  - BH-FDR correction within dataset × database × outcome.
- **Isolate level.**
  - Each tested isolate's V4 region (after 515F, before the 806R reverse-complement site, primers excluded; 5,555 of 7,873 tested isolates) was searched all-vs-all with vsearch `--iddef 1` (terminal gaps count, so a match needs ≥ 99% identity over the whole region).
  - Leave-study-out prediction used reference sets (i) all, (ii) same host species, (iii) same host genus, (iv) same country.
  - Two predictors: majority Inhibitory share of the ≥ 99% matches, and the bacterial-genus majority among reference-set isolates (genus from round-2 SINTAX/RDP on full-length sequences).
  - Kappa and AUC with 1,000 study-cluster bootstrap replicates, on common panels.

## 1. Restricted subsets
Inventory: [subset_inventory.csv]({{artifact:efdb8b7a-dc17-4ff9-9df6-5d9fdb94165f}}). Total = Full subset size; Tested = isolates with a Bd assay (Inhibitory / NoEffect / Facilitating).

| Set | Total | Tested | Inhibitory | Inh. strict | Facilitating | Studies (tested) |
|---|---|---|---|---|---|---|
| Unrestricted | 8,739 | 7,873 | 2,831 | 2,273 | 645 | 33 |
| Bullfrog: *L. catesbeianus* (Virginia, Colorado) | 295 | 269 | 97 | 63 | 39 | 2 |
| Bullfrog: Massachusetts ('MA'; all *N. viridescens*) | 137 | 134 | 34 | 24 | 14 | 1 |
| Bullfrog: species + MA | **0** | – | – | – | – | – |
| Genus fallback *Lithobates* + *Rana* | 528 | 502 | 221 | 146 | 40 | 8 |
| US Northeast (present: MA, New Jersey) | 157 | 154 | 50 | 36 | 14 | 2 |
| SERDP: *L. sphenocephalus* = Tennessee = species + TN | 73 | 73 | 6 | 5 | 0 | 1 |
| SERDP: Panama (country) | 1,180 | 1,150 | 404 | 384 | 104 | 4 |
| SERDP: PA, VT, NM, LA (state); *R. pipiens* / *R. blairi* (species) | **0** | – | – | – | – | – |

- The Panamanian leopard frog is not in AmphiBac. *R. pipiens*, *R. blairi* and the Panama site therefore fall back to *Lithobates* + *Rana*, which contains no Panama isolates.
- All 73 Tennessee / *L. sphenocephalus* isolates lack the 515F site. They start inside V4 and overlap the 150-bp SERDP reads by only about 111 bp.

**Meta-analysis coverage.**

| Restriction | Samples with a restricted metric (of 2,411) | Sets |
|---|---|---|
| Same host species | 1,681 | 83 |
| Same host + country | 1,492 | 82 |
| Same country | 2,113 | 4 |

- Same-host sets are small: median 23 tested isolates (range 1–269).
- 76 of the 83 host sets draw their tested isolates from a single study.

## 2. Restricted metrics versus full-database metrics
Tables: [restricted_vs_full_comparison.csv]({{artifact:089aff51-7450-42be-8d7b-5be3dea25ed6}}), [bullfrog_per_sample_restricted.csv]({{artifact:15e8b8d2-60b0-4c2d-9fd5-f0a28fdb6393}}), [serdp_per_sample_restricted.csv]({{artifact:964e6667-db5f-4cbf-85d1-4b960715371f}}), [meta_per_sample_restricted.csv.gz]({{artifact:30ba1506-ab0b-4c9a-9003-42e480a3c06a}}).

![Coverage lost by restriction]({{artifact:c924bf2e-d780-463f-aa19-13bf1a717d8f}})

**Coverage lost** is the share of reads matched with the full DB that become unmatched with the restricted DB. **Partial share** is the share of restricted matches that are partial overlaps.

| Set | n | Median coverage lost | Mean partial share | Spearman ρ, Inhibitory Prop. vs full |
|---|---|---|---|---|
| Bullfrog *L. catesbeianus* | 40 | 0.57 | 0.25 | 0.94 |
| Bullfrog MA | 40 | 0.90 | 0.00 | 0.57 |
| Bullfrog *Lithobates* + *Rana* | 40 | 0.55 | 0.26 | 0.93 |
| Bullfrog US Northeast | 40 | 0.88 | 0.00 | 0.75 |
| SERDP *L. sphenocephalus* (LA + TN) | 174 | 0.70 | 1.00 | 0.12 |
| SERDP *Lithobates* + *Rana* (all sites) | 372 | 0.37 | 0.14 | 0.73 |
| SERDP Tennessee (TN) | 80 | 0.70 | 1.00 | 0.33 |
| SERDP Panama | 13 | 0.64 | 0.00 | 0.40 |
| Meta same host | 1,681 | 0.97 | 0.18 | 0.03 |
| Meta same host + country | 1,492 | 0.98 | 0.19 | 0.04 |
| Meta same country | 2,113 | 0.13 | 0.02 | 0.80 |

- Meta same-host: median restricted cultured share 0.015, against 0.685 with the full database. 25% of covered samples have no restricted match at all.
- Inhibitory / Full ratio, Spearman ρ against the full-DB ratio: bullfrog −0.09 to 0.10; SERDP *L. sphenocephalus* −0.31; SERDP Tennessee −0.13; meta same host 0.11.

![Restricted vs full Inhibitory proportion]({{artifact:060d5e31-e85f-4787-b91b-b720c495c91e}})

**Partial-overlap artefact.** Because iddef 2 ignores terminal gaps, an isolate that covers only part of the read scores 100% identity on the overlap.
- With the full database, a full-coverage hit usually wins the k-mer ranking.
- In a restricted database, partial overlaps are often the only candidates left.
- For Tennessee, 100% of restricted matches are 61–124-bp overlaps of the 150-bp reads. The full-coverage any-hit Inhibitory proportion is exactly 0 for every Tennessee frog, but the pipeline reports a mean of 0.088.

**Confirmation of the exhaustive-filter derivation.** The full-coverage any-hit proportion matched the real restricted vsearch output exactly in:
- 80% of meta samples for same host (89% within 0.01);
- 84% for same host + country;
- 88% for same country;
- 100% of bullfrog samples for the MA and Northeast sets.

Where they differ, the real run has the higher proportion, because of partial-overlap hits the filter excludes (mean difference 0.007 for meta same host). The derivation is therefore a lower bound on what the pipeline reports. All per-sample metrics reported here come from the real runs.

**Permutation intervals.** Very few samples exceed the 95th percentile of their own restricted set's shuffled-label null:
- meta: 0.5% (same host), 0.4% (same host + country), 0.2% (same country);
- bullfrog: 0% for every set;
- SERDP: 5.1% for the *Lithobates* + *Rana* set, 0% for the others.

## 3. Bd tests
Results: [bd_tests_restricted.xlsx]({{artifact:c2920a2e-2335-4550-9ca4-c38e8f02d778}}) and [bd_tests_restricted.csv]({{artifact:207ed5de-4da4-4160-b946-43d5aa2210d6}}); 841 rows, 674 estimable.

**SERDP, best available match per site, pooled.**
- Infection, restricted Inhibitory proportion: OR 0.89 per SD (0.53–1.48), n = 344. Same frogs with the full DB: OR 1.02 (0.72–1.45).
- Load among infected: β 0.26 (−0.08 to 0.59), n = 100.

**SERDP, sites with a true species or location match.**
- Tennessee (species + location), infection: OR 0.80 (0.41–1.57), n = 80, 4 strata. This metric consists entirely of partial-overlap hits.
- Panama (country), load: β 0.28 (−0.73 to 1.30), n = 13. All 13 Panama frogs are infected, so no infection test is possible there.

**Bullfrog.**
- Within-treatment Spearman correlations for restricted Inhibitory proportion, all p ≥ 0.14:
  - High N, n = 25: −0.04 to −0.01.
  - Low N, n = 15: 0.33 to 0.40.
- Pond-mean OLS, MA set: β 0.22 per SD (0.02–0.42), p = 0.036, q = 0.20. The sign is positive: more "inhibitory" reads go with higher Bd.

**Multiple testing.** No restricted-metric test survives FDR correction (minimum q = 0.20 bullfrog, 0.24 SERDP). Nominal p < 0.05 hits are 16 of 340 SERDP rows and 5 of 94 bullfrog rows, about what chance would produce. The two nominally significant richness and Panama effects point in the "more inhibitory = more Bd" direction.

## 4. Isolate level: same-host references do not predict function better
Tables: [isolate_level_kappa_auc.csv]({{artifact:99f4293a-791b-43e2-9206-3ed6e722fbc4}}), [isolate_level_predictability.csv]({{artifact:e62d7397-4657-40eb-add4-38098f647a47}}), [isolate_level_predictions.csv.gz]({{artifact:fdbcaee1-b0da-488e-ad2c-8625a934d777}}).

![Isolate-level kappa]({{artifact:cd53e426-57f7-4cad-8aef-f562a9a0c17b}})

**How many isolates can be predicted at all?** Isolates with at least one ≥ 99% V4 match in another study, of 5,555 with V4:

| Reference set | Predictable isolates | % of V4 isolates |
|---|---|---|
| All | 4,881 | 88% |
| Same host genus | 874 | 16% |
| Same host species | 551 | 9.9% |
| Same country | 1,838 | 33% |

**Accuracy on the same 551 isolates** (those with a same-host match; 16 studies, 51% Inhibitory). 95% study-cluster bootstrap CIs in parentheses; Δκ is the paired difference from the all-isolates reference.

| Reference set | V4-match κ | V4-match AUC | Δκ vs all | Genus-majority κ |
|---|---|---|---|---|
| All | 0.28 (0.11–0.42) | 0.67 | – | 0.38 (0.22–0.51) |
| Same host species | 0.08 (−0.15 to 0.29) | 0.51 | −0.19 (−0.32 to −0.03) | 0.23 |
| Same host genus | 0.09 | 0.53 | −0.18 (−0.32 to −0.02) | 0.24 |
| Same country | 0.13 | 0.59 | −0.15 (−0.24 to −0.00) | 0.28 |

- Same-host genus-majority κ differs from all by Δκ −0.15 (−0.28 to −0.04).
- On the 535-isolate panel with a match under all four reference sets, the results are the same.

## Conclusion
Restricting AmphiBac to the study host or location does not fix the pipeline.
- It removes most of the cultured matches: median coverage lost is 0.55–0.98 across 9 of 11 sets, 0.37 for SERDP *Lithobates* + *Rana*, and 0.13 for meta same-country.
- It makes the Inhibitory metric depend on partial-overlap hits that vsearch counts as 99–100% identity.
- No Bd association appears.
- Among other-study isolates, same-host references predict function worse than the full database (κ 0.08 vs 0.28 on the same isolates).

Practical advice: if host-matched metrics are reported, require full query coverage (for example `--iddef 1`, or a check on `qcov`) and report how much coverage the restriction lost.

## Caveats
- Location and species metadata for isolates are sparse: there are no PA, VT, NM, LA or Massachusetts-bullfrog isolates. "Host-matched" therefore often means genus *Lithobates* + *Rana* from Virginia, Colorado, California and New Jersey.
- Most same-host sets come from a single isolate study, so host and study effects cannot be separated.
- Meta host-name matching added new synonym and spelling fixes (listed under deviations).
- Panel C of the isolate table (all isolates with any match) compares different isolate subsets across reference sets; use panels A and B for paired comparisons.

Scripts: [hostmatched_scripts.tar.gz]({{artifact:f92206f9-be9c-4c0c-b21e-702bf0b3d84c}}). Order: h00 prepare → h01 subsets → h02a/b/c vsearch → h04 metrics → h05 isolate level → h06 Bd tests → h07 figures.
