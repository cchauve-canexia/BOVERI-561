"""
Generate plots from a grid exploration
"""
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns
import sys

from analyze_variants_utils import (
    NG_RANGE,
    VAF_VAL,
    SCORE,
    W_COMP,
    VAF,
    EXP_VAF,
    TP,
    FN,
)

PREFIX = sys.argv[1]

X_VAL = 'x'
SPEC1 = '1-spec'
STATUS = 'status'

# Statistics data
DATA_FILE = os.path.join('results', f"{PREFIX}_out.tsv")
DATA_DF = pd.read_csv(DATA_FILE, sep='\t')
DATA_DF[X_VAL] = DATA_DF.apply(lambda row: 10.0 * row[SCORE] + row[W_COMP], axis=1)
DATA_DF[SPEC1] = DATA_DF.apply(lambda row: 1.0 - row['spec.'], axis=1)

# VAF data
VAF_FILE = os.path.join('results', f"{PREFIX}_vaf.tsv")
VAF_DF = pd.read_csv(VAF_FILE, sep='\t')

# Grid
NG_RANGE_LIST = [x for x in list(DATA_DF[NG_RANGE].unique())]
NB_RANGE = len(NG_RANGE_LIST)
VAF_VAL_LIST = [x for x in list(DATA_DF[VAF_VAL].unique())]
VAF_VAL_LIST.sort()
NB_VAF = len(VAF_VAL_LIST)
GRID = [(v, r) for v in VAF_VAL_LIST for r in NG_RANGE_LIST]

# Data organized by setting
DATA_GRID, VAF_GRID_LIST, VAF_GRID = {}, {}, {}
for (v, r) in GRID:
    DATA_GRID[(v, r)] = DATA_DF.loc[
        (DATA_DF[NG_RANGE]==r) & (DATA_DF[VAF_VAL]==v)
    ]
for r in NG_RANGE_LIST:
    VAF_GRID_LIST[r] = []
    for v in VAF_VAL_LIST:
        VAF_GRID_LIST[r].append(VAF_DF.loc[
            (VAF_DF[NG_RANGE]==r) & (VAF_DF[VAF_VAL]==v)
        ])
    VAF_GRID[r] = pd.concat(VAF_GRID_LIST[r])


# Plot of sensitivity, specificity and PPV per setting
out_file = DATA_FILE.replace('.tsv', '.png')
fig, axes = plt.subplots(nrows=NB_VAF, ncols=NB_RANGE, figsize=(NB_RANGE * 10, NB_VAF * 5))
i, j = 0, 0
for vaf in VAF_VAL_LIST:
    for range in NG_RANGE_LIST:
        j = j % NB_RANGE
        DATA_GRID[(vaf, range)].plot(
            x=X_VAL,
            y=['sens.', 'spec.', 'prec.'],
            title=f"Expected VAF={str(vaf).rstrip()} DNA input={str(range).rstrip()}",
            ax=axes[i, j],
            grid=True
        )
        j += 1
    i += 1
plt.savefig(out_file)

# ROC curve per setting
out_file = DATA_FILE.replace('.tsv', '_ROC.png')
fig, axes = plt.subplots(nrows=NB_VAF, ncols=NB_RANGE, figsize=(NB_RANGE * 10, NB_VAF * 5))
i, j = 0, 0
for vaf in VAF_VAL_LIST:
    for range in NG_RANGE_LIST:
        j = j % NB_RANGE
        DATA_GRID[(vaf, range)].plot(
            x=SPEC1,
            y=['sens.'],
            title=f"Expected VAF={str(vaf).rstrip()} DNA input={str(range).rstrip()}",
            ax=axes[i, j],
            grid=True,
            kind='scatter'
        )
        j += 1
    i += 1
plt.savefig(out_file)

# VAP scatter plots
out_file = VAF_FILE.replace('.tsv', '.png')
fig, axes = plt.subplots(nrows=NB_RANGE, ncols=1, figsize=(NB_RANGE * 10, NB_VAF * 5))
i = 0
for range in NG_RANGE_LIST:
    correlation = VAF_GRID[range][VAF].corr(
        VAF_GRID[range][EXP_VAF],
        method='pearson'
    )
    vaf_plot = sns.scatterplot(
        data=VAF_GRID[range],
        x=VAF, y=EXP_VAF, hue=STATUS,
        palette={TP: 'blue', FN: 'orange'},
        ax=axes[i]
    )
    title = f"DNA input={str(range).rstrip()}"
    axes[i].set_title(f"{title}, Pearson correlation={round(correlation, 3)}")
    i += 1
plt.savefig(out_file)
