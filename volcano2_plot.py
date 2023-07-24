# Import required libraries
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Read in the data
df = pd.read_csv('GFP_fulv_DEG_results.txt', sep='\t')

# Add a column for -log10(adjusted P-value)
df['-log10_adj_pval'] = -np.log10(df['adj.P.Val'])

# Create the plot
plt.figure(figsize=(10, 5))

# Plot all points in the same color
plt.scatter(df['logFC'], df['-log10_adj_pval'], s=10, color='blue', edgecolor='k')

# Label the axes
plt.xlabel('Log2 Fold Change', fontsize=14)
plt.ylabel('-Log10 Adjusted P-value', fontsize=14)
plt.title('Volcano Plot', fontsize=20)

plt.xlim([-8,8])
plt.ylim([0,25])

# Define a threshold for significant points that should be labeled
label_threshold_pval = 5
label_threshold_logFC = 4

# Add labels for significant points
for i in range(len(df)):
    if df.iloc[i]['-log10_adj_pval'] > label_threshold_pval and abs(df.iloc[i]['logFC']) > label_threshold_logFC:
        plt.text(df.iloc[i]['logFC'], df.iloc[i]['-log10_adj_pval'], df.index[i], fontsize=8)

# Display the plot
plt.show()
