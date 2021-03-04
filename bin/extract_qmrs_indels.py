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

from analyze_variants_utils import (
    SAMPLE,
    RUN_ID,
    VAF,
    CHR,
    POS,
    REF,
    ALT,
    NG,
    W_SCORE,
    W_COMP,
    SCORE,
    COMPLEXITY,
    SUPPORT,
    OVERLAP,
    CONTROL,
    add_weighted_score,
    augment_fingerprints,
    filter_blacklist,
    filter_fingerprints,
    read_blacklist,
)

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

def find_indel(df, run, sample, chr, pos, ref, alt):
    return list(df.loc[(df[RUN_ID]==run) & (df[SAMPLE]==sample) & (df[CHR]==chr) & (df[POS]==pos) & (df[REF]==ref) & (df[ALT]==alt)].index)

def get_runs_qmrs_data(
    run_id_list,
    qmrs_samples_list=None,
    blacklist=[],
    fingerprints_df=None
):
    """
    Returns the dataframes of observed indels for run in run_id_list and
    sample from sample_list from
    results files
    :param: run_id_list (list(str)): ID of the runs to consider
    :param: qmrs_samples_list (list(str)): list of the IDs of QMRS samples
    :param: blacklist (list(dict(str, str, int, str, str, str)): indexed by
    CHR, POS, REF, ALT, BL_ORIGIN
    :param: filter_artifacts (bool): if True filter artifacts from blakclist
    :param: fingerprints_df (DataFrame): dataframe of fingerprints characterized
    by columns AMP_CHR_COL, AMP_START_COL and AMP_END_COL
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
        if qmrs_samples_list is None:
            observed_indels_df = observed_df_aux.loc[
                observed_df_aux[SAMPLE].str.startswith('QMRS')
            ].round(3)
        else:
            observed_indels_df = observed_df_aux.loc[
                observed_df_aux[SAMPLE].isin(qmrs_samples_list)
            ].round(3)
        observed_indels_df_list.append(observed_indels_df)
        # Excluding indels in the blasklist
    qmrs_all_observed_indels_df = pd.concat(observed_indels_df_list)
    # qmrs_all_observed_indels_df[POS] = qmrs_all_observed_indels_df[POS].astype(int)
    qmrs_observed_indels_df_1 = filter_blacklist(
        qmrs_all_observed_indels_df, blacklist, False
    )
    qmrs_observed_indels_df_1.reset_index(drop=True, inplace=True)
    if fingerprints_df is not None:
        qmrs_observed_indels_df = filter_fingerprints(
            qmrs_observed_indels_df_1, fingerprints_df
        )
    else:
        qmrs_observed_indels_df = qmrs_observed_indels_df_1
    qmrs_observed_indels_df.reset_index(drop=True, inplace=True)
    return qmrs_observed_indels_df

def export_indels(
    qmrs_observed_indels_df, score_max, min_vaf, out_prefix, out_file
):
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
    ARGS_QMRS_LIST = ['-q', '--qmrs_samples', 'QMRS samples list file']
    ARGS_BLACKLIST = ['-b', '--blacklist', 'Blacklist file']
    ARGS_FINGERPRINTS = ['-f', '--fingerprints', 'Fingreprint genes file']
    ARGS_MANIFEST = ['-m', '--manifest', 'Manifest file']
    parser = argparse.ArgumentParser(description='Indels testing: report')
    parser.add_argument(ARGS_RUNS_FILE[0], type=str, help=ARGS_RUNS_FILE[2])
    parser.add_argument(ARGS_SCORE[0], type=float, help=ARGS_SCORE[2])
    parser.add_argument(ARGS_W_COMP[0], type=float, help=ARGS_W_COMP[2])
    parser.add_argument(ARGS_MIN_VAF[0], type=float, help=ARGS_MIN_VAF[2])
    parser.add_argument(
        ARGS_QMRS_LIST[0], ARGS_QMRS_LIST[1], type=str, help=ARGS_QMRS_LIST[2]
    )
    parser.add_argument(
        ARGS_BLACKLIST[0], ARGS_BLACKLIST[1], type=str, help=ARGS_BLACKLIST[2]
    )
    parser.add_argument(
        ARGS_FINGERPRINTS[0], ARGS_FINGERPRINTS[1], type=str,
        help=ARGS_FINGERPRINTS[2]
    )
    parser.add_argument(
        ARGS_MANIFEST[0], ARGS_MANIFEST[1], type=str, help=ARGS_MANIFEST[2]
    )
    args = parser.parse_args()
    # Reading parameters
    RUNS_FILE = open(args.runs_file, 'r').readlines()
    RUNS_LIST = []
    for run_data in RUNS_FILE:
        RUNS_LIST.append(run_data.rstrip().split(',')[1])
    if args.qmrs_samples is not None:
        QMRS_SAMPLES_LIST = [
            x.rstrip() for x in open(args.qmrs_samples, 'r').readlines()
        ]
    else:
        QMRS_SAMPLES_LIST = None
    MANIFEST_DF = pd.read_csv(args.manifest, sep='\t')
    FINGERPRINTS_DF = augment_fingerprints(
        pd.read_csv(args.fingerprints, sep='\t'), MANIFEST_DF
    )
    if args.blacklist is not None:
        BLACKLIST = read_blacklist(args.blacklist, MANIFEST_DF)
    else:
        BLACKLIST = []
    QMRS_INDELS_DF = add_weighted_score(
        get_runs_qmrs_data(
            RUNS_LIST,
            qmrs_samples_list=QMRS_SAMPLES_LIST,
            blacklist=BLACKLIST,
            fingerprints_df=FINGERPRINTS_DF
            ),
        args.w_comp
    )
    OUT_PREFIX = [args.score_max, args.w_comp]
    OUT_FILE_NAME = args.runs_file.replace(
        'data', 'results'
    ).replace(
        '.csv', f"_QMRS_{args.score_max}_{args.w_comp}_{args.min_vaf}.tsv"
    )
    OUT_FILE = open(OUT_FILE_NAME, 'w')
    export_indels(
        QMRS_INDELS_DF, args.score_max, args.min_vaf, OUT_PREFIX, OUT_FILE
    )
    OUT_FILE.close()
