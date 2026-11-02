# 16S primer and region use in amphibian microbiome studies

**V4 is the most common target, but not as dominant as assumed, and its share is falling.** We screened 524 Europe PMC records, of which 180 primary studies state the 16S region they amplified:

| Region | Studies | Share |
|---|---|---|
| V4 | 102 | 57% |
| V3–V4 | 52 | 29% |
| V4–V5 | 9 | 5% |
| Full-length (V1–V9) | 5 | 3% |
| V1–V3 | 4 | 2% |
| Seven other regions (V3–V5, V1–V2, V8–V9, V5–V6, V3, and two multi-region studies) | 8 | 4% |

V4 accounted for 80% of studies in 2015–17 but only 42% in 2024–26. Over the same period V3–V4 rose from 0% to 41%, and in gut and fecal studies it is already the most common region (55%).

![Regions targeted by amphibian microbiome studies]({{artifact:d75f4fc4-856e-469b-9635-f60a78330b5a}})

*Figure. (a) Amplified 16S spans on E. coli coordinates. The label beside each span gives the number and share of the 180 studies; bar colour picks out V4 and V3–V4. (b) Region shares by publication period; the two studies from 2014 or earlier are not shown. (c) Region shares by body site. Body site comes from the extraction; studies sampling neither skin nor gut are omitted.*

## V4 and the Earth Microbiome Project primers

V4 became the default through the 515F/806R primers of the Earth Microbiome Project. It entered amphibian work through early Illumina skin studies ([Loudon et al. 2013](https://doi.org/10.1038/ismej.2013.200); [Longo et al. 2015](https://doi.org/10.1098/rsos.140377); [Walke et al. 2015](https://doi.org/10.1371/journal.pone.0139848)). V4 still dominates skin studies (64% of 141), in line with the amphibian skin meta-analyses that standardised on it.

Reporting of the primer version is poor. Of the 102 V4 studies, only 18 identify the original Caporaso 515F/806R. Only 5 identify the degenerate Parada 515F-Y and/or Apprill 806RB primers, which improve coverage of SAR11 and Crenarchaeota. The other 78 do not say which version they used. This matters when comparing studies at the level of exact sequence variants, because the 515F-Y change alters which taxa amplify.

## V3–V4: two distinct primer traditions

The rise of V3–V4 reflects two separate traditions.

The first is the Klindworth/Herlemann pair, 341F with 785R or 805R. It is the default of Illumina's 16S library protocol and is used almost entirely in skin studies: 15 of 15 papers with this pair are skin studies (e.g. [Demircan et al. 2018](https://doi.org/10.1038/s41598-018-29373-y); [Ujszegi et al. 2023](https://doi.org/10.1186/s13568-023-01585-0)).

The second is 338F/806R, the standard of several commercial sequencing providers. It dominates gut work: 14 of its 24 papers are gut studies. Among papers with full text, V3–V4 papers are much more likely to have Chinese affiliations or providers (33 of 52) than V4 papers (21 of 100). The remaining V3–V4 papers use 341F/806R (4) or do not name the primers (9). Recent examples include [Dong et al. 2026](https://doi.org/10.3389/fvets.2026.1900113) and [Xue-Dong et al. 2026](https://doi.org/10.3389/fmicb.2026.1751143).

The two V3–V4 pairs share 806R-type reverse primers, so their reads overlap with V4 amplicons only in the V4 portion. Both also run on 2 × 250–300 bp chemistry; 10 of the 16 NovaSeq studies with a stated region are V3–V4.

## Less common regions

**V4–V5 (515F with 907R, 909R, 926R or 939R).** Nine studies, for example [Kruger 2020](https://doi.org/10.3389/fmicb.2020.01330), [Campos et al. 2024](https://doi.org/10.1016/j.crmicr.2023.100213) and [Zhang et al. 2024](https://doi.org/10.3389/fmicb.2024.1356161). A cluster of US studies uses 515F-Y/939R but describes it as "V3–V5", although both primers flank V4–V5 ([Jiménez et al. 2022](https://doi.org/10.1128/aem.01818-21); [Hauser et al. 2024](https://doi.org/10.7554/elife.92168); [Osborne et al. 2024](https://doi.org/10.1093/ismejo/wrae104); [Madison et al. 2025](https://doi.org/10.1128/iai.00569-24)).

**Full-length 16S.** PacBio or 27F/1492R amplicons have appeared only recently, in skin ([Zhu et al. 2024](https://doi.org/10.1016/j.isci.2024.109907); [Schuck et al. 2026](https://doi.org/10.1186/s42523-026-00593-2)) and gut ([Zhu et al. 2024](https://doi.org/10.1093/ismejo/wrad009)).

**Early 454 regions.** The pyrosequencing studies used the 5′ end of the gene:
- V1–V2, 27F/388R ([McKenzie et al. 2011](https://doi.org/10.1038/ismej.2011.129));
- V1–V3, 27F/533R ([Chang et al. 2016](https://doi.org/10.1186/s12866-016-0660-4); [Huang et al. 2018](https://doi.org/10.3389/fmicb.2017.02670)). Both papers call this the "V4 region", but those primers amplify V1–V3.

**Single studies.** V5–V6, 783F/1027R ([Federici et al. 2015](https://doi.org/10.1264/jsme2.me15041)), and V8–V9, V8f/1510r ([Bouchali et al. 2026](https://doi.org/10.1128/msystems.01436-25)).

**Other markers alongside 16S.** Fungal ITS is the most common (48 mentions across ITS1, ITS2 and unspecified ITS), followed by shotgun metagenomes (12) and 18S (11).

## Sequencing platforms

Illumina MiSeq is used in 130 studies, followed by NovaSeq (23) and HiSeq (15). Roche 454 (8 studies) belongs to the pre-2017 literature, and PacBio (3) is recent.

## Implications for AmphiBac-style matching

For reference matching against isolate databases such as AmphiBac, the region mix matters. A V4-only reference set covers about 57% of studies. Sequences from the other common regions sit in different parts of the gene:
- **V3–V4** shares only its V4 portion with V4 data;
- **V4–V5** extends past 806R into V5;
- **V1–V3** shares no sequence with V4 at all.

In our coverage analysis, the V3–V4 window was well represented in AmphiBac's Sanger isolate sequences, but V5 and beyond were poorly covered. So V4–V5 and the 3′ regions can only be matched against a minority of isolates.

## Methods and limits

- **Search.** We queried Europe PMC titles and abstracts for amphibian taxon terms together with skin, microbiome or microbiota terms and 16S, amplicon or sequencing-platform terms. This returned 524 records (2002–2026).
- **Full texts.** Open-access XML was retrieved for 383 records. Methods-section sentences mentioning primers, regions or platforms were extracted, and a language model classified each record: primary sequencing study or not, 16S region, primers, primer version, platform, and other markers. It also returned a supporting quote, used for checking.
- **Corrections.** Where both primers were named, the region was re-derived from their E. coli positions. This corrected 12 records, including 8 papers whose own stated region contradicts their primers.
- **Validation.** The extraction agreed with independent regex detection of 515F and 806R in 78 of 79 papers, and spot checks of quoted evidence matched.
- **Records without a stated region.** Of the 301 primary studies:
  - 98 could be assessed only from the abstract (no open full text), so their region is not assessed rather than absent;
  - 15 used only ITS, 18S or shotgun sequencing;
  - 8 did not state a region.
- **Coverage limits.** Records whose title or abstract lacks the search terms are missed, as are paywalled methods sections. Commercially sequenced gut studies are likely over-represented relative to skin studies, which are more often open access. The skin and gut split comes from the extraction, not from manual review.

## Files

- Per-study table (180 studies: region, primers, primer version, platform, evidence quote): [amphibian_microbiome_16S_regions.csv]({{artifact:6b17702e-5fee-41df-b90d-7c97942be428}})
- All 301 screened primary studies, with assessment scope: [screened_primary_studies.csv]({{artifact:1a84fc8a-bbdc-4cb5-bc3b-f8a03e5bf30f}})
- Search, download and regex scripts: [primer_survey_scripts.tar.gz]({{artifact:d9b6d077-6bc7-4ef7-a6d6-10af807facb6}}). The language-model classification step ran inside the analysis session, so it is not included.
