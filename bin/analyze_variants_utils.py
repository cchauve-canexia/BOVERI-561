#!/usr/bin/env python3

"""
Common functions and variables for analyzing results in commercial and clinical
samples
"""
from collections import defaultdict
import os
import pandas as pd
import yaml

# Dataframes column keys: should be imported from vcf_utils.py
SAMPLE = 'sample'
RUN_ID = 'run_id'
VAF = 'vaf'
EXP_VAF = 'exp_vaf'
CHR = 'chr'
POS = 'pos'
REF = 'ref'
ALT = 'alt'
NG = 'ng_est'
W_SCORE = 'w_score'
SCORE = 'score'
COMPLEXITY = 'complexity'
SUPPORT = 'support'
OVERLAP = 'overlap'
CONTROL = 'control'

LVE_KEY = 'min_exp_vaf'
LPI_KEY = 'lower_pi'

# Indel calls status
FP = 'FP'
TP = 'TP'
FN = 'FN'
FN_U = 'FN_u'
FN_O = 'FN_o'
TN = 'TN'
STATUS_KEYS = [FP, TP, FN_U, FN_O, TN]
EXPECTED = 'E'
OBSERVED = 'O'

# Settings features
W_COMP = 'w_comp'
VAF_VAL = 'vaf_values'
NG_RANGE = 'ng_range'

## Analysis: Features defining an indel
## Basic indel features
INDEL_FEATURES = [SAMPLE, CHR, POS, REF, ALT]
# Expected indels features
INDEL_FEATURES_EXPVAF = INDEL_FEATURES + [EXP_VAF]
# Detected indels features
INDEL_FEATURES_VAF_1 = INDEL_FEATURES + [VAF]
INDEL_FEATURES_VAF_2 = [W_SCORE, SCORE, COMPLEXITY, SUPPORT, OVERLAP, CONTROL]
INDEL_FEATURES_VAF = INDEL_FEATURES_VAF_1 + INDEL_FEATURES_VAF_2

# Paremeters keys for the YAML file
GRID_KEY = 'penalty_grid'
VAF_KEY = 'vaf_range'
NG_KEY = 'ng_range'
SCORES_KEY = 'scores'
W_COMP_KEY = 'w_comp'
BLACKLIST_KEY = 'blacklist_files'
OUT_SUFFIX_KEY = 'out_suffix'
MANIFEST_KEY = 'manifest_file'
FILTER_ARTIFACTS_KEY = 'filter_artifacts'
FINGERPRINT_KEY = 'fingerprint_file'

# Amplicons manifest keys
AMP_NAME_COL = 'Amplicon_name'
AMP_START_COL = 'Start'
AMP_END_COL = 'End'
AMP_CHR_COL = 'Chr'

# Origin of a blacklisted indel
BL_ARTIFACT = 'Artifact'
BL_BLACKLIST = 'blacklisted'
# Used columns of blacklist files
BL_INDEL = 'id'
BL_ORIGIN = 'origin'

# ------------------------------------------------------------------------------
# Auxiliary functions

def row_to_tuple(row):
    """
    Dataframe row encoding an indel to a tuple for indexing purpose
    """
    return (row[RUN_ID], row[SAMPLE], row[CHR], row[POS], row[REF], row[ALT])

def read_parameters(yaml_file_path):
    """
    Read a YAML parameters file and returns a dictionary of parameters values
    """
    parameters = {BLACKLIST_KEY: [], GRID_KEY: None, OUT_SUFFIX_KEY: ''}
    features_dict, blacklist_files = {}, []
    with open(yaml_file_path) as c_file:
        parameters_dict = yaml.safe_load(c_file)
        for key, value in parameters_dict.items():
            if key == GRID_KEY:
                parameters[GRID_KEY] = [
                    (float(x.split('_')[0]), float(x.split('_')[1]))
                    for x in value.split()
                ]
            elif key in [SCORES_KEY, W_COMP_KEY]:
                value_split = str(value).split()
                step = float(value_split[2])
                [r_start, r_end] = [int(x) for x in value_split[0:2]]
                features_dict[key] = [
                    round(x * step, 2) for x in range(r_start, r_end)
                ]
            elif key in [VAF_KEY, NG_KEY]:
                parameters[key] = [
                    [float(x) for x in y.split('_')] for y in str(value).split()
                ]
            elif key == BLACKLIST_KEY:
                blacklist_files = value.split()
            elif key == OUT_SUFFIX_KEY:
                parameters[key] = value
            elif key in [MANIFEST_KEY, FINGERPRINT_KEY]:
                parameters[key] = pd.read_csv(value, sep='\t')
            elif key == FILTER_ARTIFACTS_KEY:
                parameters[key] = value
    for blacklist_file in blacklist_files:
        parameters[BLACKLIST_KEY] += read_blacklist(
            blacklist_file, parameters[MANIFEST_KEY]
        )
    if parameters[GRID_KEY] is None:
        parameters[GRID_KEY] = [
            (s, w)
            for s in features_dict[SCORES_KEY]
            for w in features_dict[W_COMP_KEY]
        ]
    return parameters

def augment_fingerprints(fingerprints_df, manifest_df):
    """
    Associate to each fingerprint its coordinates
    """
    fg_df = fingerprints_df.copy()
    fg_df[AMP_CHR_COL] = ''
    fg_df[AMP_START_COL] = 0
    fg_df[AMP_END_COL] = 0
    for index, fg in fingerprints_df.iterrows():
        amplicon_name = fg[AMP_NAME_COL]
        amplicon_index = list(
            manifest_df.loc[manifest_df[AMP_NAME_COL]==amplicon_name].index
        )[0]
        amplicon = manifest_df.loc[[amplicon_index]]
        fg_df.at[index, AMP_START_COL] = int(amplicon.at[amplicon_index, AMP_START_COL])
        fg_df.at[index, AMP_END_COL] = int(amplicon.at[amplicon_index, AMP_END_COL])
        fg_df.at[index, AMP_CHR_COL] = amplicon.at[amplicon_index, AMP_CHR_COL]
    return fg_df

def coord_to_del(chr, start, end, manifest_df):
    """
    Given an interval chr:start-end, returns the corresponding sequence
    extended by one base to the left
    """
    for _, amplicon in manifest_df.iterrows():
        test_1 = amplicon['Chr'] == chr
        test_2 = int(amplicon['Start']) < start
        test_3 = int(amplicon['End']) >= end
        if test_1 and test_2 and test_3:
            pos = start - int(amplicon['Start']) - 1
            ref_len = end - start + 1
            ref = amplicon['Amplicon'][pos:pos + ref_len + 2]
            alt = amplicon['Amplicon'][pos:pos + 1]
            return {CHR: chr, POS: pos, REF: ref, ALT: alt}

def read_blacklist(blacklist_tsv_file, manifest_df, sep='\t'):
    """
    Extracts from a CSV file a list of (aliquot, chr, pos, ref, alt) describing
    indels
    :param: blacklist_csv_file (str): path to a CSV file describing indels in
    the format chr,pos,ref,alt
    :param: manifest_df (DataFrame): amplicon manifest
    :param: sep (str): separator in CSV file
    :return: list(dict(str, str, int, str, str, str)): indexed by
    CHR, POS, REF, ALT, BL_ORIGIN
    """
    blacklist = []
    blacklist_df = pd.read_csv(blacklist_tsv_file, sep='\t')
    for _, row in blacklist_df.iterrows():
        indel_str = row[BL_INDEL].split(':')
        indel = {
            BL_ORIGIN: row[BL_ORIGIN],
            CHR: indel_str[0]
        }
        if len(indel_str) == 4:
            indel[POS] = int(indel_str[1])
            indel[REF], indel[ALT] = indel_str[2], indel_str[3]
        else:
            assert indel_str[2] == 'del'
            [start, end] = [int(x) for x in indel_str[1].split(',')]
            indel_dict = coord_to_del(indel_str[0], start, end, manifest_df)
            for key in [POS, REF, ALT]:
                indel[key] = indel_dict[key]
        blacklist.append(indel)
    return blacklist

def check_blacklist(expected_df, blacklist):
    """
    Check that no expected indel is in the blacklist
    """
    index = []
    for indel in blacklist:
        index += list(
            expected_df.loc[
                (expected_df[CHR]==indel[CHR]) &
                (expected_df[POS]==indel[POS]) &
                (expected_df[REF]==indel[REF]) &
                (expected_df[ALT]==indel[ALT])
            ].index
        )
    return len(index) == 0


def filter_blacklist(df, blacklist, filter_artifacts):
    """
    Returns a dataframe of variants from df by filtering all calls that match
    the features (chr, pos, ref, alt) in blacklist
    :param: df (DataFrame): dataframe of indels with a column SAMPLE
    :param: blacklist list(dict(str, str, int, str, str, str)): indexed by
    CHR, POS, REF, ALT, BL_ORIGIN
    :param: filter_artifacts (bool): if True, indels in the blacklist whose
    origin is ARTIFACT are considered, otherwise they are not
    :return: DataFrame: input dataframe without any row where the entries
    CHR, POS, REF, ALT match a blacklisted indel from blacklist
    """
    index = []
    if not filter_artifacts:
        blacklist_aux = [x for x in blacklist if x[BL_ORIGIN] != BL_ARTIFACT]
    else:
        blacklist_aux = blacklist
    for indel in blacklist_aux:
        index += list(
            df.loc[
                (df[CHR]==indel[CHR]) &
                (df[POS]==indel[POS]) &
                (df[REF]==indel[REF]) &
                (df[ALT]==indel[ALT])
            ].index
        )
    return df.loc[~df.index.isin(index)]

def filter_fingerprints(df, fingerprints_df):
    """
    Returns a dataframe of variants from df by filtering all calls that falls
    within a fingerprint
    :param: df (DataFrame): dataframe of indels with a column SAMPLE
    :param: fingerprints_df (DataFrame): dataframe of fingerprints characterized
    by columns AMP_CHR_COL, AMP_START_COL and AMP_END_COL
    """
    index_list = []
    for _, fg in fingerprints_df.iterrows():
        index_list += list(
            df.loc[
                (df[CHR]==fg[AMP_CHR_COL]) &
                (df[POS].between(fg[AMP_START_COL], fg[AMP_END_COL]))
            ].index
        )
    return df.drop(index=index_list)

def get_runs_data(
    run_id_list,
    samples_list,
    blacklist=[],
    filter_artifacts=False,
    fingerprints_df=None
):
    """
    Returns the dataframes of observed indels for run in run_id_list and
    sample from sample_list from
    results files
    :param: run_id_list (list(str)): ID of the runs to consider
    :param: samples_list (list(str)): list of samples in the run for which there
    is at least one expected indel
    :param: blacklist (list(dict(str, str, int, str, str, str)): indexed by
    CHR, POS, REF, ALT, BL_ORIGIN
    :param: filter_artifacts (bool): if True filter artifacts from blakclist
    :param: fingerprints_df (DataFrame): dataframe of fingerprints characterized
    by columns AMP_CHR_COL, AMP_START_COL and AMP_END_COL
    :return DataFrame: dataframe of observed indels
    """
    observed_indels_df_list = []
    for run_id in run_id_list:
        # Reading the pipeline results
        run_indels_file = os.path.join(
            'results', run_id, f"{run_id}_indels.tsv"
        )
        observed_df_aux = pd.read_csv(run_indels_file, sep='\t')
        # Reformatting the sample column to match the true indels format
        observed_df_aux[SAMPLE] = observed_df_aux.apply(
            lambda row: row[SAMPLE].split('_S')[0], axis=1
        )
        observed_indels_df = observed_df_aux.loc[
            observed_df_aux[SAMPLE].isin(samples_list)
        ].round(3)
        observed_indels_df_list.append(observed_indels_df)
        # Excluding indels in the blasklist
    all_observed_indels_df = pd.concat(observed_indels_df_list)
    observed_indels_df_1 = filter_blacklist(
        all_observed_indels_df, blacklist, filter_artifacts
    )
    observed_indels_df_1.reset_index(drop=True, inplace=True)
    if fingerprints_df is not None:
        observed_indels_df = filter_fingerprints(
            observed_indels_df_1, fingerprints_df
        )
    else:
        print('No fingerprint filtering')
        observed_indels_df = observed_indels_df_1
    observed_indels_df.reset_index(drop=True, inplace=True)
    return observed_indels_df

# ------------------------------------------------------------------------------
STAT_PRECISION = 3

def compute_ratio(x, y, precision=STAT_PRECISION):
    """
    Computes the ratio x/y with precision precision, but if y=0 in whih case
    returns 0.0
    """
    return (0.0 if y == 0 else  round(x / y, precision))

OUT_MAIN = 'out_main'
OUT_ERRORS = 'out_errors'
OUT_VAF = 'out_vaf'

def open_out_files(parameters, exp_indels_file):
    """
    Open and write headers for the three output files.
    """
    out_suffix = parameters[OUT_SUFFIX_KEY]
    _, dataset_name = os.path.split(exp_indels_file)
    out_file_name = dataset_name.replace('.tsv', f"{out_suffix}_out.tsv")
    out_file = open(os.path.join('results', out_file_name), 'w')
    header_1 = [VAF_VAL, NG_RANGE, SCORE, W_COMP]
    header_2 = [TP, FP, TN, FN, FN_U, FN_O]
    header_3 = ['sens.', 'spec.', 'acc.', 'youden', 'prec.', 'recall', 'F1', 'FDR']
    out_file.write('\t'.join(header_1 + header_2 + header_3))
    out_file_errors_name = out_file_name.replace('_out.tsv', '_errors.txt')
    out_file_errors = open(os.path.join('results', out_file_errors_name), 'w')
    out_file_vaf_name = out_file_name.replace('_out.tsv', '_vaf.tsv')
    out_file_vaf = open(os.path.join('results', out_file_vaf_name), 'w')
    header_4 = ['error', RUN_ID, SAMPLE, CHR, POS, REF, ALT, VAF, EXP_VAF]
    header_5 = [W_SCORE, SCORE, COMPLEXITY, SUPPORT, OVERLAP, CONTROL]
    out_file_errors.write('\t'.join(header_1 + header_4 + header_5))
    out_file_vaf.write('\t'.join(
            header_1 + INDEL_FEATURES + [VAF, EXP_VAF, 'status']
        )
    )
    return {
        OUT_MAIN: out_file, OUT_ERRORS: out_file_errors, OUT_VAF: out_file_vaf
    }


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

def export_statistics(index_dict, out_prefix, out_file):
    # Computing statistics
    nb_tp = len(index_dict[TP])
    nb_fn_unobserved = len(index_dict[FN_U])
    nb_fn_observed = len(index_dict[FN_O])
    nb_fn = nb_fn_unobserved + nb_fn_observed
    nb_fp = len(index_dict[FP])
    nb_tn = len(index_dict[TN])
    nb_all = nb_tp + nb_tn + nb_fp + nb_fn
    sensitivity = compute_ratio(nb_tp, nb_tp + nb_fn)
    specificity = compute_ratio(nb_tn, nb_tn + nb_fp)
    youden = round(sensitivity + specificity - 1.0, STAT_PRECISION)
    accuracy = compute_ratio(nb_tp + nb_tn, nb_all)
    precision = compute_ratio(nb_tp, nb_tp + nb_fp) # = PPV
    recall = sensitivity
    F1 = compute_ratio(2.0 * precision * recall, precision + recall)
    FDR = round(1.0 - precision, STAT_PRECISION)
    # Output of statistics
    stats = [nb_tp, nb_fp, nb_tn, nb_fn, nb_fn_unobserved, nb_fn_observed]
    stats += [sensitivity, specificity, accuracy, youden, precision, recall]
    stats += [F1, FDR]
    results = [str(x) for x in out_prefix + stats]
    out_file.write('\n' + '\t'.join(results))

def export_errors(
    expected_indels_df, observed_indels_df, index_dict, out_prefix, out_file
):
    # Output of FP indels
    for index in index_dict[FP]:
        indel_info = (
            out_prefix +
            [FP,  observed_indels_df.at[index, RUN_ID]] +
            [observed_indels_df.at[index, x] for x in INDEL_FEATURES] +
            [round(observed_indels_df.at[index, VAF], 3), 'nan'] +
            [round(observed_indels_df.at[index, x], 3) for x in INDEL_FEATURES_VAF_2]
        )
        indel_str = '\t'.join([str(x) for x in indel_info])
        out_file.write('\n' + indel_str)
    # Output of observed FN indels
    for (index_exp, index_obs) in index_dict[FN_O]:
        indel_info = (
            out_prefix +
            [FN_O,  observed_indels_df.at[index_obs, RUN_ID]] +
            [observed_indels_df.at[index_obs, x] for x in INDEL_FEATURES] +
            [round(observed_indels_df.at[index_obs, VAF], 3)] +
            [round(expected_indels_df.at[index_exp, EXP_VAF], 3)] +
            [round(observed_indels_df.at[index_obs, x], 3) for x in INDEL_FEATURES_VAF_2]
        )
        indel_str = '\t'.join([str(x) for x in indel_info])
        out_file.write('\n' + indel_str)
    # Output of unobserved FN
    for index in index_dict[FN_U]:
        indel_info = (
            out_prefix +
            [FN_U,  expected_indels_df.at[index, RUN_ID]] +
            [expected_indels_df.at[index, x] for x in INDEL_FEATURES_EXPVAF] +
            ['nan'] +
            ['nan' for x in INDEL_FEATURES_VAF_2]
        )
        indel_str = '\t'.join([str(x) for x in indel_info])
        out_file.write('\n' + indel_str)

def export_vafs(
    expected_indels_df, observed_indels_df, index_dict, out_prefix, out_file
):
    index_status_list = []
    # VAF detected and expected for TP
    for (index_exp, index_obs) in index_dict[TP]:
        index_status_list.append((index_exp, index_obs, TP))
    # Output of expected and detected VAF for detected FN
    for (index_exp, index_obs) in index_dict[FN_O]:
        index_status_list.append((index_exp, index_obs, FN))
    # Writing
    for (index_exp, index_obs, status) in index_status_list:
        observed_vaf = observed_indels_df.at[index_obs, VAF]
        expected_vaf = expected_indels_df.at[index_exp, EXP_VAF]
        vaf_info = (
            out_prefix +
            [observed_indels_df.at[index_obs, x] for x in INDEL_FEATURES] +
            [observed_vaf, expected_vaf, status]
        )
        vaf_str = '\t'.join([str(x) for x in vaf_info])
        out_file.write('\n' + vaf_str)

def compute_lpi_lve(expected_indels_df):
    lpi, lve = {}, {}
    for _, row in expected_indels_df.iterrows():
        lpi[row[SAMPLE]] = float(row[LPI_KEY])
        lve[row[SAMPLE]] = float(row[LVE_KEY])
    return lpi, lve


def compute_lve(expected_indels_df):
    lpi, lve = {}, {}
    for _, row in expected_indels_df.iterrows():
        lve[row[SAMPLE]] = float(row[LVE_KEY])
    return lve
