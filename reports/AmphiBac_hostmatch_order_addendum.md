# Host-matched databases and reference-order uncertainty

## Summary

1. **Restricting AmphiBac to the sampled host species or location makes predictions worse, not better.** Most matched reads are lost, the remaining matches are often partial overlaps, and isolate labels become less predictable.
2. **For the pipeline's inhibitory metric, reference order adds no uncertainty.** Subset-mode Inhibitory_Proportion is exactly order-invariant across 10,000 reference orders.
3. **For metrics that depend on which reference is assigned, reference order adds large uncertainty.** This covers the top-hit label share and richness. The median 95% width for the top-hit inhibitory share is 0.09 in the meta-analysis and 0.25 in the leopard frogs.

## 1. Host- and location-matched databases

### What exists in AmphiBac

| Dataset | Species match | Location match | Species + location |
|---|---|---|---|
| Bullfrog, Massachusetts | *L. catesbeianus*: 269 isolates, from Virginia (142) and Colorado (127) | Massachusetts: 134 isolates, all *Notophthalmus viridescens* | none |
| SERDP *R. sphenocephala* (LA, TN) | *L. sphenocephalus*: 73 isolates, all from Tennessee | Tennessee: the same 73 | Tennessee only |
| SERDP *R. pipiens* (PA, VT), *R. blairi* (NM), Panama leopard frog | none; genus fallback is *Lithobates* + *Rana* (528 isolates) | Panama: 1,150 isolates from other hosts; none from PA, VT, NM or LA | none |
| Meta-analysis | 83 host species with isolates; 76 of these 83 sets come from a single study | per-country sets | per host × country |

### Coverage and artefacts

- **Most previously matched reads are lost.** The median share of previously matched reads that become unmatched is 0.55–0.98 in 9 of 11 restricted sets. For example: bullfrog species-only 0.57, Massachusetts-only 0.90; SERDP Tennessee 0.70, Panama 0.64; meta same host 0.97, same host + country 0.98. Only *Lithobates* + *Rana* for SERDP (0.37) and same-country for the meta-analysis (0.13) keep most reads.
- **Small restricted databases expose the partial-overlap artefact.** In the default q2-vsearch call, terminal gaps are ignored, so a reference covering only part of a read can still count as a 99% match.
  - All 73 Tennessee *L. sphenocephalus* isolates lack the 515F primer site, so every Tennessee match is a 61–124 bp overlap of a 150-bp read.
  - The pipeline nonetheless reports a mean inhibitory proportion of 0.088 for Tennessee frogs. Requiring full coverage gives 0 for every frog.
  - Partial overlaps are a median 25% of matched reads for bullfrog species-only matching.
- **Restricted metrics bear little relation to the full-database ones in the meta-analysis.** Spearman ρ with the full-database metric is 0.03 for same host, 0.04 for host + country and 0.80 for same country. With same-host matching, the median cultured share falls from 0.685 to 0.015, and 25% of covered samples have no match at all.

![Coverage lost by restriction]({{artifact:c924bf2e-d780-463f-aa19-13bf1a717d8f}})

### Bd tests

No restricted metric predicts Bd after FDR correction (minimum q 0.20 bullfrog, 0.24 SERDP).
- **SERDP, best available match per site, pooled:** infection OR 0.89 (0.53–1.48), n = 344.
- **SERDP, Tennessee species + location:** OR 0.80 (0.41–1.57), n = 80.
- **SERDP, Panama:** all 13 frogs are infected, so there is only a load test: β 0.28 (−0.73 to 1.30).
- **Bullfrog, Massachusetts database:** the pond-level slope is +0.22 (p = 0.036, q = 0.20). That is the wrong direction for protection.

### Is function more predictable within a host?

No. On the same 551 isolates (16 studies), with labels predicted from other studies only:

| Reference set | Majority of ≥99% V4 matches, κ | Majority of ≥99% V4 matches, AUC | Genus majority, κ |
|---|---|---|---|
| All other-study isolates | 0.28 | 0.67 | 0.38 |
| Same host species | 0.08 (Δκ −0.19, CI −0.32 to −0.03) | 0.51 | 0.23 |
| Same host genus | 0.09 | — | — |
| Same country | 0.13 | — | — |

Only 551 of 5,555 isolates with V4 (9.9%) even have a same-host match from another study. Restriction removes the many cross-host relatives that carry most of the (weak) information, and no host-specific functional signal replaces it.

![Same-host predictability]({{artifact:cd53e426-57f7-4cad-8aef-f562a9a0c17b}})

## 2. Uncertainty from reference order (10,000 permutations)

### How vsearch actually picks

We reverse-engineered vsearch's default choice. Candidates are ranked by number of shared unique 8-mers, then by shorter target length, then by file position, keeping 33 candidates and allowing 32 rejects.

An emulator built on this rule reproduced 100% of feature assignments in 301 real shuffled runs against the full database and 101 against the Inhibitory subset, for each dataset. A naive model that picks uniformly among references tied at the top identity is badly wrong (45–57% agreement), because the length tie-break makes most identity ties unreachable. The 10,000 permutations therefore used the validated emulator.

### Which metrics depend on order

- **Order-invariant:** the pipeline's subset-mode Inhibitory_Proportion and whether each feature matches at all. Both have range 0 across 101 real runs and 10,000 simulated orders.
- **Order-dependent:** anything that depends on which reference is assigned. This covers the top-hit label share, Full and Inhibitory richness, and the Inhibitory/Full ratio. 29–38% of matched reads have an order-dependent label.

| Dataset | Samples | Median 95% width, top-hit inhibitory share | Samples with width > 0.10 | Samples changing quartile in ≥5% of orders |
|---|---|---|---|---|
| Meta-analysis (90 bp) | 2,411 | 0.09 | 47% | 75% |
| Bullfrog (120 bp) | 40 | 0.03 | 8% | 68% |
| SERDP (150 bp) | 372 | 0.25 | 84% | 83% |

The value from the distributed file order lies inside the permutation interval for more than 99% of samples. In SERDP it is slightly low on average (median 0.029 below the permutation median).

![Per-sample intervals]({{artifact:179fb6aa-0c73-4640-b7af-7bb32d64cf64}})

### Effect on conclusions

No significance call flipped for this metric. Meta-analysis contrasts were p < 0.05 in at least 99.99% of orders; the SERDP and bullfrog tests were significant in 0%. Effect sizes did move substantially:

| Test | Default order | 95% range across orders |
|---|---|---|
| Meta, adult vs larvae (rank-biserial) | 0.41 | 0.19–0.49 |
| SERDP infection OR | 0.86 | 0.77–1.00 |
| SERDP load β | 0.09 | −0.18 to 0.14 (bimodal) |
| Bullfrog pond-level slope | 0.13 | 0.10–0.17 |

The SERDP load slope has the default order's sign in only 31% of orders. That instability comes from one abundant feature, which is called Inhibitory in 33% of orders.

![Propagated tests]({{artifact:97c07914-6cdc-4d36-99a9-6a1d7f4bb069}})

These intervals capture only the arbitrary tie-breaking. Uncertainty from labels (identical sequences with conflicting labels) and from assays is additional and larger.

## Implications

- **Do not restrict AmphiBac to the study host or location.** The database is too thin per host (mostly one study per host), and restriction amplifies the partial-coverage artefact. If restriction is used anyway, require full query coverage.
- **Report order uncertainty for any metric based on the top hit, or avoid such metrics.** Running the pipeline across shuffled reference orders gives a cheap interval. The subset-mode inhibitory proportion is immune to order, but that is because it only records whether any inhibitory isolate matched, which is why it behaves like cultured share.

## Caveats

- **Meta restricted metrics.** These came from real restricted vsearch runs for all 169 host and country sets. The host-name harmonisation added several genus spelling fixes, giving 83 host sets.
- **Isolate-level comparison.** Matching used ≥99% identity with terminal gaps counted. The genus predictor uses the SINTAX genus of each isolate's full sequence.
- **Order permutations.**
  - The Inhibitory-subset orders were permuted independently of the full-database orders.
  - Real-run validation used 301 shuffled full-database and 101 Inhibitory runs per dataset; the remaining orders were emulated.
  - For the meta-analysis, real runs used the 3,773 sOTUs with any ≥99% hit, since sOTUs without one can never match.
- **Rarefaction.** Single seeds: meta-analysis 1,000 reads, bullfrog 1,883 reads.

## Files

**Host-matched databases**
- [AmphiBac_hostmatched_addendum.md]({{artifact:ab0d8dcf-b211-46c7-80aa-9fec4af1162a}}) — detailed write-up
- [subset_inventory.csv]({{artifact:efdb8b7a-dc17-4ff9-9df6-5d9fdb94165f}})
- [restricted_vs_full_comparison.csv]({{artifact:089aff51-7450-42be-8d7b-5be3dea25ed6}})
- Per-sample restricted metrics: [bullfrog]({{artifact:15e8b8d2-60b0-4c2d-9fd5-f0a28fdb6393}}), [SERDP]({{artifact:964e6667-db5f-4cbf-85d1-4b960715371f}}), [meta]({{artifact:30ba1506-ab0b-4c9a-9003-42e480a3c06a}})
- [bd_tests_restricted.xlsx]({{artifact:c2920a2e-2335-4550-9ca4-c38e8f02d778}})
- [isolate_level_kappa_auc.csv]({{artifact:99f4293a-791b-43e2-9206-3ed6e722fbc4}})
- [restricted vs full figure]({{artifact:060d5e31-e85f-4787-b91b-b720c495c91e}})
- [hostmatched_scripts.tar.gz]({{artifact:f92206f9-be9c-4c0c-b21e-702bf0b3d84c}})

**Order permutations**
- [AmphiBac_order_permutation_addendum.md]({{artifact:ae58e62f-34d8-4d11-af03-e88e1f026ea2}}) — detailed write-up
- [AmphiBac_order_permutation_results.xlsx]({{artifact:7d523c2a-5932-4770-ae11-3f4fb79ff64e}})
- Per-sample intervals: [meta]({{artifact:78591eb5-a7fe-4994-a03b-8b82d31f3a65}}), [bullfrog]({{artifact:eb9cf079-a0b9-4ff7-8509-64ed5bb8f733}}), [SERDP]({{artifact:cf0d2321-7305-4ea5-aa96-3dc38816cc48}})
- Per-feature call probabilities: [meta]({{artifact:af0f9f00-90b9-4614-b222-17a774a02709}}), [bullfrog]({{artifact:129fd148-a0be-4cf0-b036-632d177bec48}}), [SERDP]({{artifact:7fb1dd51-b55c-4d71-a001-259070edae39}})
- [emulator validation figure]({{artifact:4007b1e9-9213-4eca-8bb9-8903e7d30c99}}), [emulator_validation.csv]({{artifact:7a149312-44da-491d-9e47-b3aca73dfa56}})
- [propagated_tests.csv]({{artifact:551ea260-09ca-4de1-a6ee-cd509f4e7923}})
- [AmphiBac_order_permutation_scripts.tar.gz]({{artifact:a5a4d160-9968-4431-be43-33dcfdb81550}})
