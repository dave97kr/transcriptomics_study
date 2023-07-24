# Install the required libraries if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("limma")
BiocManager::install("fgsea")

# Load the required libraries
library(limma)
library(fgsea)

# Load the limma results
results <- read.table("GFP_fulv_DEG_results.txt", header=TRUE, sep="\t", row.names=1)

# Check if the table has the correct format
head(results)

# Prepare the gene list
results$P.Value[is.na(results$P.Value)] <- 1  # replace NA values with 1
results$logP <- -log10(results$P.Value)
results$signedLogP <- results$logP * sign(results$logFC)
results <- results[order(results$signedLogP, decreasing = TRUE),]  # order by signedLogP in decreasing order

geneList <- results$signedLogP
names(geneList) <- rownames(results)

# Function to read GMT file
read.gmt <- function(file){
  con <- file(file, "r")
  pathways <- list()
  while(TRUE) {
    line <- readLines(con, n = 1)
    if(length(line) == 0) {
      break
    }
    line <- strsplit(line, "\t")[[1]]
    pathways[[line[1]]] <- line[3:length(line)]
  }
  close(con)
  return(pathways)
}

# Use the function
pathways <- read.gmt("REACTOME_ONCOGENIC_MAPK_SIGNALING.v2023.1.Hs.gmt")

# Modify pathways to only contain genes that are present in geneList
pathways <- lapply(pathways, function(x) x[x %in% names(geneList)])

# Filter out empty pathways
pathways <- pathways[sapply(pathways, length) > 0]

# Run GSEA
fgseaRes <- fgsea(pathways, geneList, minSize=15, maxSize=500)

# Convert fgseaRes into a data frame
fgseaRes_df <- as.data.frame(fgseaRes)

# Convert the leadingEdge column to character format
fgseaRes_df$leadingEdge <- sapply(fgseaRes_df$leadingEdge, function(x) paste(x, collapse = ", "))

# Write the data frame to a file
write.table(fgseaRes_df, file = "fgsea_results.txt", sep = "\t", row.names = FALSE)

# Start the graphical device
png("fgsea_plot.png", width=800, height=600)

# Limit the number of pathways to a maximum of 10
numPathways <- min(nrow(fgseaRes), 10)

# Plot the results (top 10 pathways)
topPathways <- head(fgseaRes[order(fgseaRes$padj), ], numPathways)$pathway
plotGseaTable(pathways[topPathways], geneList, fgseaRes, gseaParam = 0.5)

# Close the graphical device
dev.off()
