import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
# Read the data from the text file
res = pd.read_csv("GFP_fulv_deseq_results.txt", sep="\t")

# Remove rows with NA in the 'padj' or 'log2FoldChange' column
res = res.dropna(subset=['padj', 'log2FoldChange'])

# Adjust p-values for plotting (avoid -inf after log10 transformation)
res.loc[res['padj'] == 0, 'padj'] = np.finfo(float).eps

# Create a basic volcano plot
plt.figure(figsize=(10, 5))
plt.scatter(res['log2FoldChange'], -np.log10(res['padj']), alpha=0.5)
plt.xlabel('Log2 Fold Change', fontsize=14)
plt.ylabel('-Log10 Adjusted p-value', fontsize=14)
plt.title('Volcano Plot', fontsize=20)
plt.grid(True)
plt.show()