# How robust are AmphiBac-based Bd-function predictions?

Stress tests of the AmphiBac 2025.1 closed-reference workflow: fragment length, identity threshold, top-hit selection, database subsets and edge cases.

## 1. What the pipeline does

The recommended workflow is `tools/antifungal_predictions_pipeline.py`; the older equivalent is `scripts/Antifungal_Predictions_conda.sh`. It runs

```
qiime vsearch cluster-features-closed-reference --p-perc-identity 0.99 --p-strand both
```

against one of four subsets: FullDatabase (8,739 sequences), Inhibitory (2,831), InhibitoryStrict (2,273) or Facilitating (645). Per sample it reports:

- `<Mode>_Total`: matched reads.
- `<Mode>_Proportion`: Total divided by the rarefying depth.
- `<Mode>_Richness`: the number of matched reference IDs.

q2-vsearch passes vsearch these options: `--usearch_global --id 0.99 --strand both --qmask none --minseqlength 1`. It sets no `--maxaccepts`, `--maxrejects` or coverage options, so vsearch's defaults apply:

- `maxaccepts 1`: stop at the first candidate that passes the threshold.
- `maxrejects 32`: give up after 32 failed candidates.
- `iddef 2`: terminal gaps are excluded from the identity calculation.
- No minimum query coverage.

All tests below call vsearch 2.32.0 directly with exactly these arguments. QIIME 2 itself was not run.

**Test data:**

- The AmphiBac 2025.1 qza files. Their contents are identical to the 2025.1 FASTA.
- The amphibian meta-analysis: 48,431 deblur sOTUs of 90 bp, in 2,411 samples with ≥1,000 reads.
- The bullfrog data: 3,700 ASVs of 120 bp, in 40 usable samples.
- In-silico amplicons cut from AmphiBac isolates at 90, 100, 120 and 150 bp, full V4 (~253 bp) and V3–V4 (~429 bp).

**Ground truth:** the isolate's `AntiFungal.Function.with.AP.function` label. This is the column that defines the subsets.

## 2. Headline findings

| # | Issue | Magnitude | Consequence |
|---|---|---|---|
| 1 | **Database labels are a hard cutoff on a continuous score.** Inhibitory means CFS score >80, Facilitating ≤−20, NoEffect otherwise. | 790 isolates that the source study called inhibitory are NoEffect. | Isolates scoring 79 and 81 get opposite labels. |
| 2 | **Identical V4 sequences do not share a function.** | Among identical 515F fragments with ≥2 labelled isolates, 74% of groups (90 bp), 73% (120 bp) and 69% (253 bp) have conflicting Bd labels. | No 16S-identity rule can be accurate. This is the same within-sOTU spread as in our panel A. |
| 3 | **Subset ("Inhibitory") mode counts any feature with one inhibitory relative as inhibitory.** | Meta-analysis: 143/291 Inhibitory-mode sOTUs (49%) have mostly (<50% inhibitory) non-inhibitory ≥99% isolates. Bullfrog: 133/191 (70%). 56% of meta-analysis "inhibitory" reads come from such features, and 78% come from features that *also* hit the Facilitating subset. | The inhibitory proportion is inflated and mixes inhibitory with facilitating. |
| 4 | **InhibitoryStrict does not fix this.** It removes only Inhibitory sequences with a 100% hit to a Facilitating isolate (558, reproduced exactly). | 748 of the 2,273 Strict sequences are 100% identical to a NoEffect isolate. 72% of Strict-matched meta-analysis sOTUs (89% of reads) also match NoEffect isolates. Strict vs Inhibitory: per-sample ρ = 0.996. | Strict gives almost the same answer as Inhibitory. |
| 5 | **Prediction accuracy is low.** | Leave-self-out at 120 bp/99%, Inhibitory mode: sensitivity 0.94, **specificity 0.15**, PPV 0.40 (inhibitory prevalence 0.38), κ = 0.07. Default top-hit label: κ = 0.15. Majority vote over all ≥99% hits: κ = 0.32. **Leave-study-out:** κ = 0.09–0.12 for every rule (95% CIs reach ≈0.03). | The best rule's apparent skill mostly reflects isolates from the same study sharing labels. |
| 6 | **The default "top hit" is arbitrary among ties.** | The default never picks a lower-identity reference (0/553 meta, 0/259 bullfrog), but 61% of matched meta-analysis sOTUs (85% of matched reads) and 77% of bullfrog ASVs tie with ≥2 references at the top identity. Across 10 shuffles of reference order, the assigned reference changed for 34% of meta sOTUs (67% of matched reads) and the inhibitory yes/no call for 14% (26% of matched reads). Bullfrog: 51% and 18%. | Function and Richness depend on the order of the FASTA file. Thread count changed nothing (0 changes). |
| 7 | **Identity thresholds mean different things at different lengths.** | 99% allows 0 / 1 / 1 / 1 / 2 / 4 mismatches at 90 / 100 / 120 / 150 / 253 / 429 bp. **At 90 bp, 99% is the same as 100%**: every per-sample output is identical. 97% allows 2 / 3 / 3 / 4 / 7 / 12 mismatches. | "99% identity" is not comparable across studies with different read lengths. |
| 8 | **Terminal-gap exclusion lets partial references pass at "100%".** | A reference covering only 20–40 bp of a 90-bp query reports 100% identity and passes 0.99 (24/25 and 25/25 synthetic tests). In real data, 22/553 meta-analysis default hits align over only 31–87 bp of the 90-bp sOTU, and 12 hits are to references with <90 bp after 515F. An exhaustive search adds 4,284 meta sOTUs whose only "hits" are partial overlaps (median 19 bp). | Partial-coverage hits will grow if anyone raises `maxrejects`. Coverage should be required explicitly. |
| 9 | **Strand.** 5,451 of the 7,443 references that carry a 515F site are stored reverse-complemented. | With q2-vsearch's own default (`--strand plus`), 294/553 meta-analysis hits and 85/259 bullfrog hits are lost, and the label changes for 333/553. | The AmphiBac scripts do set `both`. Anyone running q2-vsearch or another tool manually must too. |
| 10 | **Other edge cases.** | (a) N and compatible IUPAC codes score as matches; 7 meta default hits involve such codes. (b) 61 references of 2–8 bp contain no 8-mer and can never be hit. (c) 11% of default top hits at 90–150 bp are NotTested/PossibleContamination isolates. (d) Richness counts reference IDs, not features: 285 → 277 (meta) and 190 → 145 (bullfrog) for Inhibitory; 270 → 164 for meta Facilitating. | Small individually, but each biases the outputs. |
| 11 | **The inhibitory proportion mostly measures the "cultured share".** | Meta-analysis: per-sample Inhibitory_Proportion vs Full_Proportion ρ = 0.947. A median 92% of AmphiBac-matched reads per sample are "inhibitory" under the pipeline. | The metric mostly tracks how much of a community has a cultured relative in AmphiBac. |
| 12 | **Biological conclusions can change with the labelling rule.** | The adult > larvae difference holds under every metric. *A. boreas* vs *L. catesbeianus* (the two most-sampled named species; "Boophis sp" n=170 skipped as genus-level) is n.s. under the pipeline and every threshold or subset variant (p = 0.25–0.36). It is significant under majority vote (p = 4e-6), ≥75% consensus (p = 4e-12) and score weighting (p = 1e-6). Madagascar > USA under the pipeline (p = 3e-28) **reverses** under the ≥75% consensus rule (p = 9e-7). | Per-sample inhibitory metrics are not robust to reasonable method choices. These contrasts ignore study and site confounding; they demonstrate method sensitivity, not biology. |

Things that turned out **not** to matter much:

- Rarefaction seed: median |Δ| = 0.010, 95th percentile 0.036, ρ = 0.998.
- Thread count.
- `maxrejects 32` missing genuine full-coverage matches: 5/537 meta-analysis sOTUs, 0 bullfrog.
- Keeping IUPAC-containing fragments: κ changed by <0.01.
- Longer fragments help modestly when isolates from the same study are allowed. For the top-hit rule, κ rises from 0.10 (90 bp) to 0.22 (V4) to 0.31 (V3–V4) at 99%. Leave-study-out, length adds at most 0.06–0.09 κ for any rule.

## 3. Figures

**Database audit.** (a) Sequence available after 515F. (b) The cutoff that defines the function labels. (c) Label conflicts among identical V4 fragments. (d) What InhibitoryStrict removes.

![Database audit]({{artifact:b7464cfe-bbf8-4b68-8705-c7eff66252f4}})

**Default top hit vs exhaustive search.** Meta-analysis sOTUs (top row) and bullfrog ASVs (bottom row).

![Top hit vs all qualifying hits]({{artifact:eaa04451-a001-4e1c-b14b-cf70ffb44660}})

**Reference-order, thread and subset sensitivity.**

![Order and subset sensitivity]({{artifact:78bc355b-db86-4f64-b7f2-d3a75b5aaf01}})

**Edge cases:** terminal gaps, IUPAC codes, short references, strand.

![Edge cases]({{artifact:874b4b84-6fd0-48c4-818e-d07171d442a0}})

**Leave-one-out agreement (Cohen's κ)** by length × threshold × rule, with 95% study-cluster bootstrap CIs.

![Kappa curves]({{artifact:ad797c86-021b-4658-9095-64c039f08558}})

**Sensitivity, specificity and PPV at 99%.**

![Sens/spec/PPV]({{artifact:5cbc7b69-e8e9-40ec-adf9-886afa3086c5}})

Heatmaps: [leave-self-out]({{artifact:a5ed28ba-1cf8-4121-a7f0-831e6bc6739d}}), [leave-study-out]({{artifact:a1af3ef8-a9ce-4ae0-9c3f-97c53e9c5a01}}).

**What "inhibitory" reads actually are.** Pipeline-inhibitory reads binned by the inhibitory share of all ≥99%-identical tested isolates.

![Inhibitory read composition]({{artifact:a6185c00-32ba-40a4-8f96-840217fe9300}})

**Per-sample agreement with the pipeline output under alternative settings.**

![Rank agreement]({{artifact:6ac98f33-4db9-4c6e-b9ef-84b4ddbacd3d}})

**Biological contrasts under each metric** (meta-analysis).

![Contrasts by method]({{artifact:895214e2-5fa4-46c5-93ed-fb711a3e6ed2}})

Per-sample method scatter: [fig_rd1_method_scatter.png]({{artifact:6551f4fd-0950-4339-b8fb-074195f4e36d}}).

## 4. Recommendations for users and for the database

1. **Report a match in mismatches, not percent identity, for the read length used.** For example, "≤1 mismatch over 120 bp" rather than "99%". At 90 bp, 99% means an exact match.
2. **Use all qualifying hits, not the first one.** Run vsearch directly with `--maxaccepts 0 --maxrejects 0` and `--userout`; q2-vsearch closed-reference does not expose these options. Then summarise every isolate ≥ threshold: n, % inhibitory, % facilitating, median score. Report the conflict rather than hiding it.
3. **Require full query coverage.** Add `--query_cov 1.0` or filter `alnlen == qlen`, because iddef 2 ignores terminal gaps.
4. **Do not use subset mode as a function call.** At minimum, also report Full_Proportion (the cultured share) and the Inhibitory/Full ratio. A feature should count as inhibitory only with a stated consensus rule, e.g. ≥75% of ≥3 isolates from ≥2 studies.
5. **Use the continuous score where possible**, rather than the >80 / ≤−20 cutoffs.
6. **Exclude** NotTested/PossibleContamination isolates, references <100 bp, and references that do not span the amplicon.
7. **Always use `--strand both`** with this database, because most sequences are stored reverse-complemented.
8. **Validate across studies.** A leave-study-out κ of about 0.1 is the realistic expectation for a binary inhibitory call from 16S identity, at any length or threshold tested.
9. **Suggested database improvements:**
   - Orient all sequences.
   - Drop placeholder sequences of 2–8 bp.
   - Build InhibitoryStrict against NoEffect isolates as well as Facilitating ones.
   - Ship per-sequence consensus tables, with the continuous scores.
   - Fix the auto-filled `Assay.Method` strings ("97 well CFS assay" … "223 well CFS assay").

## 5. Methods notes and caveats

- **Matching:** vsearch 2.32.0 with the exact q2-vsearch call. Exhaustive searches used `--maxaccepts 0 --maxrejects 0` at 0.97 and were filtered to exact identity (alnlen − mism − gaps)/alnlen. "Qualifying" hits are ≥ threshold and span the full query, unless stated otherwise. Exhaustive searches still pass through vsearch's k-mer prefilter.
- **Leave-one-out:**
  - Self-exclusion used `--self`.
  - Study exclusion used per-study databases with the default call.
  - Unmatched isolates were scored as non-inhibitory, which is how the pipeline behaves; matched-only metrics are in the workbook.
  - Ties in majority vote went to non-inhibitory.
  - CIs come from a study-cluster bootstrap (B = 1,000, 28 studies). The V3–V4 fragment used the Klindworth 341F primer (CCTACGGGNGGCWGCAG).
- **Real-data runs:**
  - The meta-analysis was rarefied to 1,000 reads (seeds 1 and 2); the bullfrog data to 1,883 reads, the minimum depth after dropping P10F3 and P4F1.
  - The hit map from unrarefied features was applied to the rarefied tables (verified identical, 285/285).
  - Samples with no match are reported as 0; the real pipeline's left-merge may give NaN.
  - In `meta_contrasts.csv`, the Kruskal–Wallis rows store their p-values in the column named `p_MWU`.
- **Ground truth:** isolate labels are in-vitro CFS or agar results from many labs and Bd strains. Disagreement partly reflects assay noise, not only strain-level biology. Results using `Bd.Inhibition` as truth were similar (e.g. Inhibitory mode, 120 bp/99%: κ = 0.084).
- **Meta-analysis sOTUs:** 232 sOTUs with lowercase sequence IDs were upper-cased and kept, giving 48,431 sOTUs.

## 6. Files

- Data and tables:
  - [audit_per_sequence.csv]({{artifact:76d77534-2987-4b44-99f0-acd51045b8ce}})
  - [audit_identical_V4_groups.csv]({{artifact:5e4d146f-ca89-4e3c-b0f7-36f9c727bfa9}})
  - [tophit_summary.xlsx]({{artifact:4059a0e7-156e-4368-8fdd-b4b6c5ef185b}})
  - [meta_per_query_tophit.csv]({{artifact:9046d550-57de-43c9-812b-76a0df897fdc}})
  - [bull_per_query_tophit.csv]({{artifact:db7919a0-7aa2-46fa-a570-b16dd323ee6a}})
  - [lt_length_threshold_summary.xlsx]({{artifact:dca8f85d-9cb7-452e-a77c-26d9476052e0}})
  - [lt_metrics_long.csv]({{artifact:8c242b36-35d5-4e77-9ea4-c2b73691cabd}})
  - [lt_mismatch_table.csv]({{artifact:12f3afda-f4b1-4728-ab84-a002a2f21761}})
  - [realdata_impact_summary.xlsx]({{artifact:c0bc7610-92fa-444e-9876-a604fe00dd95}})
  - [meta_per_sample_methods.csv]({{artifact:c99f3b59-0fb2-40d3-9194-db86479dee0c}})
  - [bf_per_sample_methods.csv]({{artifact:ff449495-bf24-4712-bb39-89fc00548481}})
  - [meta_feature_labels.csv]({{artifact:bc52136c-455b-4227-a907-5ac426ce676f}})
  - [bf_feature_labels.csv]({{artifact:8da5e7fb-c613-4550-9d83-07f5f6448a7d}})
  - [meta_contrasts.csv]({{artifact:5f48020a-0ee0-47be-85ab-6c980dd04b0c}})
- Inputs: [amphibac_2025.1_testbundle.tar.gz]({{artifact:635e3d25-3992-409c-9c87-382965c05cb1}}) holds the database FASTAs, metadata and audit table.
- Scripts: [AmphiBac_robustness_scripts.tar.gz]({{artifact:33a0a17c-c8e2-46d1-9ed3-15d1899b5395}}) contains the audit, top-hit, length/threshold and real-data scripts.
