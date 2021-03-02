#!/usr/bin/env bash

# python bin/extract_files.py data/BOVERI-568_commercial.csv
# python bin/dump_variants.py data/BOVERI-568_commercial.csv -m CG001v5.1_Amplicon_Manifest_Panel5.1.12_20200911.tsv
# python bin/extract_files.py data/BOVERI-568_clinical.csv
# python bin/dump_variants.py data/BOVERI-568_clinical.csv -m CG001v5.1_Amplicon_Manifest_Panel5.1.12_20200911.tsv

python bin/check_indels.py data/v51NextSeq_commercial_samples_expected_indels_NO_WT.tsv \
  data/BOVERI-568_NextSeq_blacklist_2602.tsv BOVERI-568_commercial-report.yaml \
  > data/BOVERI-568_NextSeq_blacklist_2602_commercial_stats.tsv

python bin/analyze_variants.py data/v51NextSeq_commercial_samples_expected_indels_NO_WT.tsv BOVERI-568_commercial-report.yaml
python bin/analyze_variants.py data/v51NextSeq_commercial_samples_expected_indels_WT.tsv BOVERI-568_commercial-report.yaml
python bin/analyze_variants.py data/v51NextSeq_commercial_Horizon_samples_expected_indels_NO_WT.tsv BOVERI-568_commercial_Horizon-report.yaml
python bin/analyze_variants.py data/v51NextSeq_commercial_SeraSeq_samples_expected_indels_NO_WT.tsv BOVERI-568_commercial_SeraSeq-report.yaml
python bin/extract_qmrs_indels.py data/BOVERI-568_NextSeq_commercial.csv 0.9 0.5 0.25

python bin/check_indels.py data/v51NextSeq_clinical_samples_expected_indels_NO_WT.tsv \
  data/BOVERI-568_NextSeq_blacklist_2602.tsv BOVERI-568_clinical-report.yaml \
  > data/BOVERI-568_NextSeq_blacklist_2602_clinical_stats.tsv

python bin/analyze_variants_clinical.py data/v51NextSeq_clinical_samples_expected_indels_NO_WT.tsv BOVERI-568_clinical-report.yaml
python bin/analyze_variants_clinical.py data/v51NextSeq_clinical_samples_expected_indels_WT.tsv BOVERI-568_clinical-report.yaml
python bin/extract_qmrs_indels.py data/BOVERI-568_NextSeq_clinical.csv 0.9 0.5 0.25

python bin/analyze_variants.py data/v51NextSeq_commercial_samples_expected_indels_NO_WT.tsv BOVERI-568_commercial-report-grid.yaml
