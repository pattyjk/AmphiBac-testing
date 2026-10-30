# AmphiBac robustness, round 3: longer fragments and different 16S regions

This addendum follows [round 1]({{artifact:68283292-5522-4f52-929b-09571edb7201}}) and [round 2]({{artifact:8887dd6f-7443-4aef-bcb0-4245efda3b44}}). It reruns the isolate-level tests from rounds 1–2 on ten 16S windows longer than 250 bp.

Tests that need community amplicons cannot be rerun at these lengths, because the available amplicon data are only 90–150 bp. This excludes:

- the per-sample pipeline outputs;
- the label-permutation null;
- the Bd-load tests.

## Design

All 8,739 database sequences were aligned to *E. coli* 16S (J01859). Each window was then cut on *E. coli* coordinates, so every window is defined the same way for every isolate. Window choice is described in the [coverage analysis]({{artifact:16e859fe-f259-4d69-8ae2-3505be53950c}}).

| Window | *E. coli* coordinates | Isolates with the fragment |
|---|---|---|
| V4 | 533–786 | 6,372 |
| V5–V6 | 780–1060 | 2,393 |
| V4–V5 | 533–906 | 2,120 |
| V2–V3 | 101–510 | 3,992 |
| V3–V4 (429 bp) | 357–785 | 6,148 |
| V6–V8 | 960–1390 | 2,270 |
| V3–V4 (500 bp) | 296–795 | 5,959 |
| V2–V4 | 201–800 | 5,095 |
| V1–V4 | 121–820 | 4,101 |
| V2–V8 | 236–1335 | 1,849 |

**Fixed panel.** Longer and more 3′ windows keep different studies and a different label mix. To separate the effect of the window from the effect of which isolates are included, every metric is also reported on a fixed panel. The panel is the 1,397 isolates present in all ten windows: 17 studies, 52% inhibitory. Window comparisons below use the panel unless stated otherwise.

## Results

![Window synthesis]({{artifact:acecbd8d-553e-4fa7-9bc4-9cbb3e50d459}})

| Test (fixed panel) | Across the ten windows | Compared with V4 |
|---|---|---|
| Leave-study-out κ at 99%, R1 (default top hit) | 0.13–0.23 | V4 0.154; no window differs with a paired CI excluding 0 |
| R2/R3 (any hit to Inhibitory / InhibitoryStrict) | κ 0.15–0.31; specificity 0.23–0.44 | Still mainly sensitivity with low specificity |
| R4 (majority of ≥99% full-coverage hits) | κ 0.27–0.44 | V4 0.32. A modest CI-supported gain for V1–V4, V2–V4, V4–V5 and V2–V8 (0.44, +0.11) |
| Genus majority (no sequence matching; window-independent upper bound) | κ 0.45 | Matches or beats every sequence rule in every window; V2–V8 R4 is within CI (−0.01) |
| Between-study P(partner inhibitory \| focal inhibitory) for pairs ≥99.5% identical | 0.77–0.83 (random pairs 0.51) | The same in every window |
| Between-study κ for identical fragments | 0.29–0.49 | No clear gain with length or region |
| Fritz–Purvis D (isolate tips) | 0.52–0.58 (within-study null ~0.93–1.00) | Weak clumping, identical across windows |
| consenTRAIT τD | No consistent deep clustering | Significant against the within-study null only for V2–V8 (p = 0.023) |
| Prospective prediction (isolates added after 2020.6 / 2023.2) | κ −0.23 to 0.22 (cells with n ≥ 30) | No window rescues it |

**Mechanics**

- **Mismatches allowed at 99%.** The threshold allows 2 mismatches (V4, V5–V6) up to 11 (V2–V8), verified from the alignments.
- **Partial-coverage hits.** Because terminal gaps are ignored, 36–77% of qualifying hits only partly cover the fragment in the V2–V3, V2–V4, V1–V4, V4–V5, V5–V6 and V2–V8 windows. Requiring full coverage raises R2 κ by up to +0.22 (V2–V8) but barely changes R1 (≤ +0.05).
- **Ties.** 79–99% of matched queries tie at the top identity. Reference-order shuffles change the default top-hit inhibitory call for 7–39% of isolates; V4 is the worst (35–39%).

## Interpretation

1. **Longer fragments and other regions do not change the round-1/2 conclusion.** On a fixed set of isolates, near-identical 16S sequences agree on Bd inhibition about 80% of the time when one is inhibitory, regardless of region or length. That is the same as genus-level agreement. Knowing the genus predicts as well as any 99% matching rule in any window.
2. **The apparent "full-length helps" result from round 2 is an isolate-set effect.** On all isolates, between-study agreement for identical fragments rises from κ 0.03 (V4, n = 6,372) to 0.41 (V2–V8, n = 1,849). On the fixed panel, V4 already gives 0.42. The long-read studies simply agree with each other more.
3. **Majority vote over all full-coverage hits (R4) is the only matching rule that approaches the genus baseline**, and it improves modestly with longer windows. The default pipeline rules (R1 and the subset rules) stay at κ ≈ 0.13–0.31 everywhere.
4. **Long windows make the coverage problem worse.** Partial-coverage references dominate qualifying hits in long windows, so requiring full query coverage becomes essential.

## Caveats and declared limits

- **Genus baseline.** It uses each isolate's full-sequence genus (RDP v18 SINTAX), so it is the same for every window and is an optimistic upper bound.
- **Full-coverage variants of R1–R3.** These are post-hoc filters on the default top hit. The exhaustive full-coverage variants are exact.
- **V4 window.** It is defined by *E. coli* coordinates (533–786), so its isolate set differs slightly from the round-1 515F-anchored V4.
- **Bootstrap and metrics.** Bootstrap B = 500, study-cluster resampling; metrics are binary only.
- **Trees and signal statistics.** Single FastTree ML trees (−nt −gtr −gamma). D comes from a vectorised re-implementation of caper::phylo.d, verified against caper. Unique-tip nulls used 199 permutations × 5 draws.
- **Prospective test.**
  - The fixed panel contains no added isolates, so a substitute panel was used: added isolates present in all six 5′-side windows (198 after 2020.6, 70 after 2023.2).
  - V4–V5, V5–V6 and V2–V8 have fewer than 30 added isolates.
  - Cells rest on 1–5 studies, so study-cluster CIs are uninformative.

## Files

- Synthesis: [fig_r3_window_synthesis.png]({{artifact:acecbd8d-553e-4fa7-9bc4-9cbb3e50d459}}), [r3_window_synthesis_panel.csv]({{artifact:e4665aa1-6e8a-4540-8839-8f34cdb3212b}}), [fragment bundle]({{artifact:305d0f1b-add2-4c0a-b447-f74675b90744}}).
- Prediction grid:
  - reports and figures: [addendum]({{artifact:44ac3b79-1300-486d-a40d-c8cbbdded4f4}}), [κ by window figure]({{artifact:dd89cb72-8f34-4fad-935a-482f1a737216}}), [mechanics figure]({{artifact:330f9358-36e6-43cd-a131-bb746bb66439}});
  - tables: [tidy summary]({{artifact:10e47e16-f1ae-4b3a-8bfb-bb059b2a63bf}}), [long metrics]({{artifact:b9f24575-6135-4b28-bb05-15cf515554f5}}), [paired contrasts]({{artifact:05a5f1cc-d416-4a32-8c80-21beb691e0e0}}), [xlsx]({{artifact:921588d3-1d6b-46a5-93b9-723bd664c464}});
  - [scripts]({{artifact:603fd98b-a8be-422d-b936-0de50d9d70fc}}).
- Decay and signal:
  - reports and figures: [addendum]({{artifact:3d42e016-7434-4e60-8cce-1aff8546649d}}), [decay κ]({{artifact:64aa7d93-e628-481b-83a4-c423f4a570dc}}), [P(partner) by identity]({{artifact:c80ee6ec-5f05-47ca-a243-190e8a240650}}), [D/τD]({{artifact:324b6feb-6d1b-4e17-af02-56b1f476dc33}});
  - tables and trees: [tidy summary]({{artifact:ee6c4bc6-37fb-4a9f-a070-e714b5085bf2}}), [xlsx]({{artifact:0535790a-e640-43d0-a091-559c1079883f}}), [trees]({{artifact:090c193c-5368-4041-899a-6283e91a79a3}});
  - [scripts]({{artifact:2e7f4b77-7a9f-42f6-87a2-28aa43d4531e}}).
- Prospective: [summary]({{artifact:095c2ceb-b157-4b1d-bae4-9013e1f8000f}}), [figure 99%]({{artifact:bb0a5051-f21f-4b30-9964-135c5a651fef}}), [xlsx]({{artifact:fbb2fab4-9ad3-4c89-a2da-afd9839f66a3}}), [scripts]({{artifact:d0c3d1ff-d25b-47bf-9930-e3ea56bd17a6}}).
