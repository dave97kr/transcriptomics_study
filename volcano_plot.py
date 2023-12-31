#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# List of files
files = ['gse121411_ex1_withcolumnandindex.txt', 'gse121411_ex2_withcolumnandindex.txt','gse121411_ex3_withcolumnandindex.txt']
dfs = []

for file in files:
    # load your data
    df = pd.read_csv(file, sep='\t') # Change this if your separator is different

    # Find columns that partially match your identifiers
    control_cols = df.filter(regex='GFP|wt|D845A').columns
    experiment_cols = df.filter(regex='S653C|L755S|V777L|L869R').columns

    # Calculate mean TPM for control and experimental groups
    control_df = df[control_cols].mean(axis=1)
    experiment_df = df[experiment_cols].mean(axis=1)

    # Calculate log fold change
    df['logFC'] = np.log2(experiment_df / control_df + 1)

    # Let's assume you have the q-value column named as 'qvalue' in the dataframe. If not, replace 'qvalue' with the correct column name
    df['-log10(qvalue)'] = -np.log10(df['qvalue'])

    dfs.append(df)

# Concatenate all dataframes
df_concat = pd.concat(dfs)

# Generate the volcano plot
plt.figure(figsize=(10,8))
plt.scatter(df_concat['logFC'], df_concat['-log10(qvalue)'], color='gray')
plt.title('Volcano plot', fontsize=20)
plt.xlabel('Log2 Fold Change', fontsize=15)
plt.ylabel('-Log10(Q value)', fontsize=15)

# Highlight points with high fold change and low q value
plt.scatter(df_concat['logFC'][(df_concat['qvalue']<0.05) & (abs(df_concat['logFC'])>1)], df_concat['-log10(qvalue)'][(df_concat['qvalue']<0.05) & (abs(df_concat['logFC'])>1)], color='red')

plt.show()