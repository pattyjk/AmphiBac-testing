AmphiBac robustness tests - export package
Start with AmphiBac_robustness_full_summary.md; all of its links are relative paths inside this folder.
reports/  per-round reports (round 1, round 2, round 3). Their internal links point to the original artifact store; the same files are in this package by filename.
figures/  all figures (PNG, 300 dpi)
tables/   result tables (csv/xlsx)
scripts/  analysis scripts (tar.gz archives per track, plus score_diptest.py)
inputs/   amphibac_2025.1_testbundle.tar.gz (AmphiBac 2025.1 subset FASTAs, metadata, audit table) and amphibac_window_fragments.tar.gz (16S window fragments)
Not included: raw community data (meta-analysis, bullfrog and SERDP tables), which are the user's own inputs.
Tools: vsearch 2.32.0, Python 3.11 (pandas, numpy, scipy, statsmodels, scikit-learn, diptest), FastTree/mafft and R for phylogenetic signal.
