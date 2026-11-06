AmphiBac paired-study export: reads and metadata used in the paired culture + amplicon test (summary section 8i)
Built 2026-11-06. Total size 10.8 GB, 471 files. MANIFEST.tsv lists every file with byte size, plus md5 for read files.

Studies used in the analysis
1. walke2015/ - Walke et al. 2015 AEM (isolates reused by Walke 2017; see walke2017/ for its extra tables)
   reads/SRR2162774.fastq.gz: single pooled run from SRA (SRP062395), R1 only, ~250 nt, V4 after 515F, primers and barcodes absent.
   Walke2015_global_mapping_file_user_supplied.csv: your mapping file (63 samples with barcodes). The SRA reads carry no barcodes, so this run
   cannot be split by sample; it was analysed pooled.
   isolates.fasta / isolates_genbank.gb: KM186980-KM187640 plus AmphiBac-only isolates; isolate_assays.csv: per-isolate CFS score/label from AmphiBac.
2. rebollar2019/ - Rebollar et al. 2019 ISME J (PRJNA521543)
   reads/: raw multiplexed SRR8570210 (R1), SRR8570211 (R2), SRR8570212 (8-nt index).
   reads_demultiplexed/: R01..R61 _R1/_R2 fastq.gz, split here by exact index match with all 8 index bases Q>=30 (barcodes in rebollar_barcodes.csv;
   61 barcodes above the largest abundance break). Assigned 8,908,591 of 15,789,745 spots (56.4%); counts in demultiplex_counts.csv.
   No public mapping file: R01..R61 have NO host labels and may include non-study samples. The analysis used R1 only.
3. bletz2018/ - Bletz et al. 2018 Proc B (PRJNA477390)
   reads/: 266 single-end R1 fastq.gz (only R1 deposited); reads_manifest.csv and samples.csv map runs to samples (242 wild + 24 Experiment 1).
   RSPB-2018-0758_Bletz_etal_Dryad_AllData.xlsx (Dryad), isolate_assays.csv (Bd and Bsal proportional growth), isolates*.fasta / .gb.
4. molecol2019/ - Ellison et al. 2019 Mol Ecol (Dryad doi:10.5061/dryad.8b83gn8)
   reads/2018_Ellison_Rsierrae_microbiome_Illumina_sequences.fna.gz: NOT fastq. Dryad provides only a demultiplexed QIIME split_libraries FASTA
   without quality scores (V3-V4 341F/805R). The analysis used only merged reads (>=400 nt, 28%).
   mapping file, Bd results (per-sample Bd status/ZE), CFS raw OD results, isolate_assays.csv (scores derived here from OD; provisional), isolates.fasta.

Also included
- asv_tables/tables_pkg/<study>/: ASV sequences (asvs.fa), ASV x sample tables, per-ASV calls (own-isolate label + every pipeline rule), logs.
- analysis/: paired_asv_agreement.xlsx (results), paired_culture_amplicon_studies.xlsx (study curation), report, analysis scripts.
Each study folder has README.txt (sources, primers, what is missing), samples.csv and ena_filereport.tsv where applicable.
Not included: Jimenez 2019, PNAS 2019 and Assis 2020 reads (downloaded but not usable for the test, so not used).
