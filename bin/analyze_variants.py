#!/usr/bin/env python3

"""
Analyzing pipeline calls compared to expected calls to generate accuracy
statistics and identify False Positive (FP) and False Negative (FN) cases
"""
import argparse
from collections import defaultdict
import pandas as pd

from analyze_variants_utils import (
    STATUS_KEYS,
    GRID_KEY,
    VAF_KEY,
    NG_KEY,
    BLACKLIST_KEY,
    FILTER_ARTIFACTS_KEY,
    FINGERPRINT_KEY,
    MANIFEST_KEY,
    EXPECTED,
    OBSERVED,
    CHR,
    W_SCORE,
    OVERLAP,
    CONTROL,
    SAMPLE,
    RUN_ID,
    VAF,
    EXP_VAF,
    NG,
    TP,
    FP,
    TN,
    FN,
    FN_O,
    FN_U,
    OUT_MAIN,
    OUT_ERRORS,
    OUT_VAF,
    row_to_tuple,
    compute_lpi_lve,
    add_weighted_score,
    export_statistics,
    export_errors,
    export_vafs,
    read_parameters,
    augment_fingerprints,
    open_out_files,
    check_blacklist,
    get_runs_data,
)


def compare_indels(
    expected_df, observed_df, indel_2_idx, score_max, lpi, lve
):
    """
    Computes the index of true positive, false positive, false negative
    :param: expected_df (DataFrame): dataframe of all expected indels
    :param: observed_df (DataFrame): dataframe of all dected indels
    :param: indel_2_idx (dict): dictionary indexed by the tuple representation
    of indels and linking an indel to a dictionary indexed by EXPECTED and
    OBSERVED each linking to a list of index in expected_df and observed_df
    corresponding to entries in these dataframes for the indel
    :param: score_max (float): max penalty to call a detected indel
    :param: lpi (dict(str, float)): dictionary of lower prediction interval
    value for each sample
    :param: lve (dict(str, float)):  dictionary of lower value expected for
    each sample
    :return: dictionary index by STATUS_KEYS
    - index of TP calls in expected_df and observed_df
    - index of observed FN indels in expected_df and in observed_df
    - index of unobserved FN indels in expected_df
    - index of FP calls in observed_df
    - index of TN calls in observed_df
    Method:
    https://docs.google.com/document/d/15hFggOFMKXyzYGcWnSeLeUNNNA-p92Wy3pTCQyvjwF4
    """
    index_dict = {status: [] for status in STATUS_KEYS}
    for index, row in expected_df.iterrows():
        indel_tuple = row_to_tuple(row)
        indel_keys = indel_2_idx[indel_tuple].keys()
        if EXPECTED in indel_keys and OBSERVED in indel_keys:
            indel_idx = indel_2_idx[indel_tuple][OBSERVED]
            score = observed_df.at[indel_idx, W_SCORE]
            overlap = observed_df.at[indel_idx, OVERLAP]
            control = observed_df.at[indel_idx, CONTROL]
            test_overlap_control = overlap < 1.0 and control < 1.0
            test_call = test_overlap_control and score <= score_max
            vaf = observed_df.at[indel_idx, VAF]
            sample = row[SAMPLE]
            if test_call and vaf >= lpi[sample]:
                index_dict[TP].append((index, indel_idx))
            if (not test_call) or vaf < lpi[sample]:
                index_dict[FN_O].append((index, indel_idx))
        if EXPECTED in indel_keys and OBSERVED not in indel_keys:
             index_dict[FN_U].append(index)
    index_fp, index_tn = [], []
    for index, row in observed_df.iterrows():
        indel_tuple = row_to_tuple(row)
        indel_keys = indel_2_idx[indel_tuple].keys()
        if EXPECTED not in indel_keys:
            vaf, score, sample = row[VAF], row[W_SCORE], row[SAMPLE]
            overlap, control = row[OVERLAP], row[CONTROL]
            test_overlap_control = overlap < 1.0 and control < 1.0
            test_call = test_overlap_control and score <= score_max
            if test_call and  vaf >= lve[sample]:
                index_dict[FP].append(index)
            if (not test_call) and vaf >= lve[sample]:
                index_dict[TN].append(index)
    return index_dict

def process_indels(
    all_expected_indels_df, all_runs_indels_df,
    penalty_grid, settings_grid, out_files
):
    """
    Computing sensitivity, specificity, precision, recall, F1, FDR,
    accuracy youden's index for a grid of thresholds and parameters
    :param: all_expected_indels_df (DataFrame): all non-WT expected indels
    :param: all_runs_indels_df (DataFrame): all observed indels
    :param: penalty_grid (list((float, float))): list of pairs (score, w_comp)
    where score is the maximum penalty score to call an indel and w_comp the
    weight given to sequence complexity
    :param: settings_grid (list((float, float))): list of pairs (vaf_values,
    ng_values) where vaf_values is a list of expected VAF and ng_values a list
    of two float representing a min ng and max ng to consider
    :param: out_files (dict(files)): dictionary of opened files indexed by
    [OUT_MAIN, OUT_ERRORS, OUT_VAF]
    Write three files; If exp_indels_file is NAME.tsv:
    - NAME_out.tsv: statistics
    - NAME_errors.tsv: FP and FN indels
    - NAME_vaf.tsv: observed and expected VAFs for TP/FN indel calls
    """
    # Dictionary of indels to index in both dataframes
    print('Creating indels index')
    # indel_2_index (dict): dictionary indexed by the tuple representation
    # of indels and linking an indel to a dictionary indexed by EXPECTED and
    # OBSERVED each linking to a list of index in all_expected_df and
    # all_runs_indels_df corresponding to entries in these dataframes for the
    # indel
    indel_2_index = defaultdict(dict)
    for index, row in all_expected_indels_df.iterrows():
        indel_2_index[row_to_tuple(row)][EXPECTED] = index
    LPI, LVE = compute_lpi_lve(all_expected_indels_df)
    for index, row in all_runs_indels_df.iterrows():
        indel_2_index[row_to_tuple(row)][OBSERVED] = index
    # Loop on scoring schemes
    for (score_max, w_comp) in penalty_grid:
        print(f"SCORE={score_max}\tW_COMP={w_comp}")
        # Reducing the weight of complexity penalty by factor w
        runs_indels_df = add_weighted_score(all_runs_indels_df, w_comp)
        for (vaf_values, ng_range) in settings_grid:
            vaf_values_str = "{:<10}".format('_'.join([str(x) for x in vaf_values]))
            ng_range_str = f"{ng_range[0]}_{ng_range[1]}"
            print(f"\tNG={ng_range_str}\tVAF={vaf_values_str}")
            out_prefix = [vaf_values_str, ng_range_str, score_max, w_comp]
            # Expected indels
            expected_indels_df = all_expected_indels_df.loc[
                (all_expected_indels_df[EXP_VAF].isin(vaf_values)) &
                (all_expected_indels_df[NG].between(ng_range[0], ng_range[1]))
            ].copy()
            # Samples containing expected indels
            samples_list = list(pd.unique(expected_indels_df[SAMPLE]))
            # Observed indels
            observed_indels_df = runs_indels_df.loc[
                runs_indels_df[SAMPLE].isin(samples_list)
            ].copy()
            # Computing TP, FP, TN, FN indexes
            index_dict = compare_indels(
                expected_indels_df, observed_indels_df, indel_2_index,
                score_max, LPI, LVE
            )
            # Computing and exporting statistics
            export_statistics(index_dict, out_prefix, out_files[OUT_MAIN])
            # Exporting errors
            export_errors(
                expected_indels_df, observed_indels_df,
                index_dict, out_prefix, out_files[OUT_ERRORS]
            )
            # Exporting VAFs
            export_vafs(
                expected_indels_df, observed_indels_df,
                index_dict, out_prefix, out_files[OUT_VAF]
            )

# ------------------------------------------------------------------------------
if __name__ == "__main__":
    # Analysis parameters
    # Input file
    ARGS_EXPECTED_FILE = ['exp_indels_file', None, 'Expected indels file']
    # Black list file
    ARGS_PARAMETERS_FILE = ['parameters_file', None, 'Parameters YAML file']
    parser = argparse.ArgumentParser(description='Indels testing: report')
    parser.add_argument(ARGS_EXPECTED_FILE[0],
                        type=str,
                        help=ARGS_EXPECTED_FILE[2])
    parser.add_argument(ARGS_PARAMETERS_FILE[0],
                        type=str,
                        help=ARGS_PARAMETERS_FILE[2])
    args = parser.parse_args()
    # Reading parameters
    PARAMETERS = read_parameters(args.parameters_file)
    PENALTY_GRID = PARAMETERS[GRID_KEY]
    VAF_RANGES = PARAMETERS[VAF_KEY]
    NG_RANGES = PARAMETERS[NG_KEY]
    SETTINGS_GRID = [(vaf, ng) for vaf in VAF_RANGES for ng in NG_RANGES]
    PARAMETERS[FINGERPRINT_KEY] = augment_fingerprints(
        PARAMETERS[FINGERPRINT_KEY],
        PARAMETERS[MANIFEST_KEY]
    )
    # Creating and opening output files
    OUT_FILES = open_out_files(PARAMETERS, args.exp_indels_file)
    # Reading all expected indels
    print('Reading indels')
    ALL_EXPECTED_INDELS_DF = pd.read_csv(args.exp_indels_file, sep='\t')
    ALL_EXPECTED_INDELS_DF.rename(columns={'chromosome': CHR}, inplace=True)
    # List of runs and samples with at least one expected indel
    RUNS_LIST = list(ALL_EXPECTED_INDELS_DF[RUN_ID].unique())
    SAMPLES_LIST = list(ALL_EXPECTED_INDELS_DF[SAMPLE].unique())
    # Extracting all expected indels and all detected indels for the
    # selected samples, but the black-listed ones
    EXPECTED_INDELS_DF = ALL_EXPECTED_INDELS_DF.loc[
        ALL_EXPECTED_INDELS_DF[SAMPLE].isin(SAMPLES_LIST)
    ]
    assert check_blacklist(EXPECTED_INDELS_DF, PARAMETERS[BLACKLIST_KEY])
    RUNS_INDELS_DF = get_runs_data(
        RUNS_LIST, SAMPLES_LIST,
        blacklist=PARAMETERS[BLACKLIST_KEY],
        filter_artifacts=PARAMETERS[FILTER_ARTIFACTS_KEY],
        fingerprints_df=PARAMETERS[FINGERPRINT_KEY]
    )
    # Processing indels for all parameters and threshold settings
    print('Processing indels')
    process_indels(
        EXPECTED_INDELS_DF.copy(), RUNS_INDELS_DF.copy(),
        PENALTY_GRID, SETTINGS_GRID, OUT_FILES
    )
    for out_file in OUT_FILES.values():
        out_file.close()
