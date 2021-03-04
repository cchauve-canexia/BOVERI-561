#!/usr/bin/env python3

"""
Check the presence of a set of indels in observed indels
"""
import argparse
import numpy as np
import os
import pandas as pd
import yaml

# Dataframes column keys: should be imported from vcf_utils.py
from analyze_variants_utils import (
    MANIFEST_KEY,
    SAMPLE,
    RUN_ID,
    VAF,
    EXP_VAF,
    CHR,
    POS,
    REF,
    ALT,
    NG,
    W_SCORE,
    SCORE,
    COMPLEXITY,
    SUPPORT,
    OVERLAP,
    CONTROL,
    INDEL_FEATURES,
    INDEL_FEATURES_EXPVAF,
    INDEL_FEATURES_VAF_1,
    INDEL_FEATURES_VAF_2,
    INDEL_FEATURES_VAF,
)

# Paremeters keys for the YAML file
M = 'manifest_file'
# Used columns of blacklist files
INDEL_ID = 'id'

# ------------------------------------------------------------------------------
# Auxiliary functions

def read_parameters(yaml_file_path):
    """
    Read a YAML parameters file and returns a dictionary of parameters values
    """
    parameters = {}
    with open(yaml_file_path) as c_file:
        parameters_dict = yaml.safe_load(c_file)
        for key, value in parameters_dict.items():
            if key == MANIFEST_KEY:
                parameters[key] = pd.read_csv(value, sep='\t')
    return parameters

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

def read_indels(indels_tsv_file, manifest_df, sep='\t'):
    """
    Extracts from a CSV file a list of (chr, pos, ref, alt) describing
    indels
    :param: blacklist_csv_file (str): path to a CSV file describing indels in
    the format chr,pos,ref,alt
    :param: manifest_df (DataFrame): amplicon manifest
    :param: sep (str): separator in CSV file
    :return: list(dict(str, str, int, str, str)): indexed by
    CHR, POS, REF, ALT
    """
    indels = []
    indels_df = pd.read_csv(indels_tsv_file, sep='\t')
    for _, row in indels_df.iterrows():
        indel_str = row[INDEL].split(':')
        indel = {
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
        indels.append(indel)
    return indels

def get_runs_data(run_id_list, samples_list):
    """
    Returns the dataframes of observed indels for run in run_id_list and
    sample from sample_list from
    results files
    :param: run_id_list (list(str)): ID of the runs to consider
    :param: samples_list (list(str)): list of samples in the run for which there
    is at least one expected indel
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
    return all_observed_indels_df


# ------------------------------------------------------------------------------
if __name__ == "__main__":
    # Analysis parameters
    # Input file
    ARGS_EXPECTED_FILE = ['exp_indels_file', None, 'Expected indels file']
    ARGS_TSV_INDELS_FILE = ['tsv_indels_file', None, 'Checked indels file']
    # Black list file
    ARGS_PARAMETERS_FILE = ['parameters_file', None, 'Parameters YAML file']
    parser = argparse.ArgumentParser(description='Indels testing: report')
    parser.add_argument(ARGS_EXPECTED_FILE[0],
                        type=str,
                        help=ARGS_EXPECTED_FILE[2])
    parser.add_argument(ARGS_TSV_INDELS_FILE[0],
                        type=str,
                        help=ARGS_TSV_INDELS_FILE[2])
    parser.add_argument(ARGS_PARAMETERS_FILE[0],
                        type=str,
                        help=ARGS_PARAMETERS_FILE[2])
    args = parser.parse_args()
    # Reading parameters
    PARAMETERS = read_parameters(args.parameters_file)
    # Reading all expected indels
    ALL_EXPECTED_INDELS_DF = pd.read_csv(args.exp_indels_file, sep='\t')
    ALL_EXPECTED_INDELS_DF.rename(columns={'chromosome': CHR}, inplace=True)
    # Reading checked indels file
    CHECKED_INDELS_DF = pd.read_csv(args.tsv_indels_file, sep='\t')
    # List of runs and samples with at least one expected indel
    RUNS_LIST = list(ALL_EXPECTED_INDELS_DF[RUN_ID].unique())
    SAMPLES_LIST = list(ALL_EXPECTED_INDELS_DF[SAMPLE].unique())
    RUNS_INDELS_DF = get_runs_data(
        RUNS_LIST, SAMPLES_LIST
    )
    # Processing indels for all parameters and threshold settings
    print('indel\tnb_occ\tmin_vaf\tmax_vaf\tmean_vaf\tmin_ctrl\tmax_ctrl\tmean_ctrl')
    for _, indel_row in CHECKED_INDELS_DF.iterrows():
        indel_1 = indel_row[INDEL_ID].split(':')
        indel = {
            CHR: indel_1[0],
            POS: int(indel_1[1]),
            REF: indel_1[2],
            ALT: indel_1[3]
        }
        indel_df = RUNS_INDELS_DF.loc[
            (RUNS_INDELS_DF[CHR]==indel[CHR]) &
            (RUNS_INDELS_DF[POS]==indel[POS]) &
            (RUNS_INDELS_DF[REF]==indel[REF]) &
            (RUNS_INDELS_DF[ALT]==indel[ALT])
        ]
        nb_occurrences = len(indel_df.index)
        max_vaf = round(np.max(indel_df[VAF]), 2)
        min_vaf = round(np.min(indel_df[VAF]), 2)
        mean_vaf = round(np.mean(indel_df[VAF]), 2)
        max_control = round(np.max(indel_df[CONTROL]), 2)
        min_control = round(np.min(indel_df[CONTROL]), 2)
        mean_control = round(np.mean(indel_df[CONTROL]), 2)
        print(f"{indel[CHR]}:{indel[POS]}:{indel[REF]}:{indel[ALT]}\t{nb_occurrences}\t{min_vaf}\t{max_vaf}\t{mean_vaf}\t{min_control}\t{max_control}\t{mean_control}")
