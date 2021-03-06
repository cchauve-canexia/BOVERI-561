#!/usr/bin/env bash

# echo "Extract results for QMRS FFPE runs"
# python bin/extract_files.py data/QMRS_FFPE.csv
# python bin/dump_variants.py data/QMRS_FFPE.csv -m CG001v5.1_Amplicon_Manifest_Panel5.1.12_20200911.tsv

echo "QMRS FFPE: QMRS samples 0.25"
python bin/extract_qmrs_indels.py data/QMRS_FFPE.csv 0.9 0.5 0.25 \
  -q data/QMRS_FFPE_samples_list.txt  \
  -f data/fingerprint_genes_map.tsv \
  -m data/CG001v5.1_Amplicon_Manifest_Panel5.1.12_20200911.tsv \
  -b data/BOVERI-568_NextSeq_blacklist_0303.tsv
