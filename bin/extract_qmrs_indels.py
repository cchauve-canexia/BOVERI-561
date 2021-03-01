#!/usr/bin/env python3

"""
Extract QMRS called indels for a given score thresholds
"""
import argparse
from collections import defaultdict
import csv
import os
import pandas as pd
import yaml

# Dataframes column keys: should be imported from vcf_utils.py
SAMPLE = 'sample'
RUN_ID = 'run_id'
VAF = 'vaf'
CHR = 'chr'
POS = 'pos'
REF = 'ref'
ALT = 'alt'
NG = 'ng_est'
W_SCORE = 'w_score'
W_COMP = 'w_comp'
SCORE = 'score'
COMPLEXITY = 'complexity'
SUPPORT = 'support'
OVERLAP = 'overlap'
CONTROL = 'control'

# ------------------------------------------------------------------------------
# Global input

## Analysis: Features defining an indel
## Basic indel features
INDEL_FEATURES = [SAMPLE, CHR, POS, REF, ALT]
# Detected indels features
INDEL_FEATURES_VAF_1 = INDEL_FEATURES + [VAF]
INDEL_FEATURES_VAF_2 = [W_SCORE, SCORE, COMPLEXITY, SUPPORT, OVERLAP, CONTROL]
INDEL_FEATURES_VAF = INDEL_FEATURES_VAF_1 + INDEL_FEATURES_VAF_2

# ------------------------------------------------------------------------------
# Auxiliary functions

def get_runs_qmrs_data(run_id_list):
    """
    Returns the dataframes of observed indels for run in run_id_list and
    sample from sample_list from
    results files
    :param: run_id_list (list(str)): ID of the runs to consider
    :return DataFrame: dataframe of observed indels in QMRS samples
    """
    observed_indels_df_list = []
    for run_id in run_id_list:
        # Reading the pipeline results
        run_indels_file = os.path.join('results', run_id, f"{run_id}_indels.tsv")
        observed_df_aux = pd.read_csv(run_indels_file, sep='\t')
        # Reformatting the sample column to match the true indels format
        observed_df_aux[SAMPLE] = observed_df_aux.apply(
            lambda row: row[SAMPLE].split('_S')[0], axis=1
        )
        observed_indels_df = observed_df_aux.loc[
            observed_df_aux[SAMPLE].str.startswith('QMRS')
        ].round(3)
        observed_indels_df_list.append(observed_indels_df)
        # Excluding indels in the blasklist
    qmrs_observed_indels_df = pd.concat(observed_indels_df_list)
    qmrs_observed_indels_df.reset_index(drop=True, inplace=True)
    return qmrs_observed_indels_df

def add_weighted_score(indels_df, weight):
    """
    Add a weighted score column to an observed indels dataframe
    """
    result_indels_df = indels_df.copy()
    # Reducing the weight of complexity penalty by factor w
    result_indels_df[W_SCORE] = result_indels_df.apply(
        lambda row: (
            weight * row[COMPLEXITY] +
            row[SUPPORT] +
            row[OVERLAP] +
            row[CONTROL]
        ),
        axis=1
    )
    return result_indels_df

def export_indels(qmrs_observed_indels_df, score_max, min_vaf, out_prefix, out_file):
    # Output of FP indels
    header_1 = [SCORE, W_COMP]
    header_4 = [RUN_ID, SAMPLE, CHR, POS, REF, ALT, VAF]
    header_5 = [W_SCORE, SCORE, COMPLEXITY, SUPPORT, OVERLAP, CONTROL]
    out_file.write('\t'.join(header_1 + header_4 + header_5))
    for _, indel in qmrs_observed_indels_df.iterrows():
        if indel[W_SCORE] <= score_max and indel[VAF] >= min_vaf:
            indel_info = (
                out_prefix +
                [indel[RUN_ID]] +
                [indel[x] for x in INDEL_FEATURES] +
                [round(indel[VAF], 3)] +
                [round(indel[x], 3) for x in INDEL_FEATURES_VAF_2]
            )
            indel_str = '\t'.join([str(x) for x in indel_info])
            out_file.write('\n' + indel_str)

# ------------------------------------------------------------------------------
if __name__ == "__main__":
    # Analysis parameters
    # Input file
    ARGS_RUNS_FILE = ['runs_file', None, 'Runs file']
    # Scoring scheme
    ARGS_SCORE = ['score_max', None, 'Score max']
    ARGS_W_COMP = ['w_comp', None, 'Complexity weight']
    ARGS_MIN_VAF = ['min_vaf', None, 'Minimum calling VAF']
    parser = argparse.ArgumentParser(description='Indels testing: report')
    parser.add_argument(ARGS_RUNS_FILE[0], type=str, help=ARGS_RUNS_FILE[2])
    parser.add_argument(ARGS_SCORE[0], type=float, help=ARGS_SCORE[2])
    parser.add_argument(ARGS_W_COMP[0], type=float, help=ARGS_W_COMP[2])
    parser.add_argument(ARGS_MIN_VAF[0], type=float, help=ARGS_MIN_VAF[2])

    args = parser.parse_args()
    # Reading parameters
    RUNS_FILE = open(args.runs_file, 'r').readlines()
    RUNS_LIST = []
    for run_data in RUNS_FILE:
        RUNS_LIST.append(run_data.rstrip().split(',')[1])
    QMRS_INDELS_DF = add_weighted_score(get_runs_qmrs_data(RUNS_LIST), args.w_comp)
    OUT_PREFIX = [args.score_max, args.w_comp]
    OUT_FILE_NAME = args.runs_file.replace('data', 'results').replace('.csv', f"_QMRS_{args.score_max}_{args.w_comp}_{args.min_vaf}.tsv")
    OUT_FILE = open(OUT_FILE_NAME, 'w')
    export_indels(QMRS_INDELS_DF, args.score_max, args.min_vaf, OUT_PREFIX, OUT_FILE)
    OUT_FILE.close()
