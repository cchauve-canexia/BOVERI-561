#!/usr/bin/env bash

# echo "Extract results for commercial samples"
# python bin/extract_files.py data/BOVERI-568_commercial.csv
# python bin/dump_variants.py data/BOVERI-568_commercial.csv -m CG001v5.1_Amplicon_Manifest_Panel5.1.12_20200911.tsv

echo "Commercial samples: stats on blacklist"
python bin/check_indels.py data/v51NextSeq_commercial_samples_expected_indels_NO_WT.tsv \
  data/BOVERI-568_NextSeq_blacklist_2602.tsv BOVERI-568_commercial-report.yaml \
  > data/BOVERI-568_NextSeq_blacklist_2602_commercial_stats.tsv

echo "Commercial samples: non-wildtype samples"
python bin/analyze_variants.py data/v51NextSeq_commercial_samples_expected_indels_NO_WT.tsv BOVERI-568_commercial-report.yaml
echo "Commercial samples: wildtype samples"
python bin/analyze_variants.py data/v51NextSeq_commercial_samples_expected_indels_WT.tsv BOVERI-568_commercial-report.yaml
echo "Commercial samples: non-wildtype Horizon samples"
python bin/analyze_variants.py data/v51NextSeq_commercial_Horizon_samples_expected_indels_NO_WT.tsv BOVERI-568_commercial_Horizon-report.yaml
echo "Commercial samples: non-wildtype SeraSeq samples"
python bin/analyze_variants.py data/v51NextSeq_commercial_SeraSeq_samples_expected_indels_NO_WT.tsv BOVERI-568_commercial_SeraSeq-report.yaml
echo "Commercial samples: QMRS samples"
python bin/extract_qmrs_indels.py data/BOVERI-568_NextSeq_commercial.csv 0.9 0.5 0.25

# echo "Extract results for clinical samples"
# python bin/extract_files.py data/BOVERI-568_clinical.csv
# python bin/dump_variants.py data/BOVERI-568_clinical.csv -m CG001v5.1_Amplicon_Manifest_Panel5.1.12_20200911.tsv

echo "Clinical samples: stats on blacklist"
python bin/check_indels.py data/v51NextSeq_clinical_samples_expected_indels_NO_WT.tsv \
  data/BOVERI-568_NextSeq_blacklist_2602.tsv BOVERI-568_clinical-report.yaml \
  > data/BOVERI-568_NextSeq_blacklist_2602_clinical_stats.tsv

echo "Clinical samples: non-wildtype samples"
python bin/analyze_variants_clinical.py data/v51NextSeq_clinical_samples_expected_indels_NO_WT.tsv BOVERI-568_clinical-report.yaml
echo "Clinical samples: wildtype samples"
python bin/analyze_variants_clinical.py data/v51NextSeq_clinical_samples_expected_indels_WT.tsv BOVERI-568_clinical-report.yaml
echo "Clinical samples: QMRS samples"
python bin/extract_qmrs_indels.py data/BOVERI-568_NextSeq_clinical.csv 0.9 0.5 0.25

echo "Commercial samples: thresholds/parameters grid exploration"
python bin/analyze_variants.py data/v51NextSeq_commercial_samples_expected_indels_NO_WT.tsv BOVERI-568_commercial-report-grid.yaml
