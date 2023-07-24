# Load the library into the R environment
library(limma)

# Read in the expression data without assigning row names initially
data1 <- read.table('gse121411_ex1_withcolumnandindex.txt', sep = '\t', header = TRUE, row.names = NULL)
data2 <- read.table('gse121411_ex2_withcolumnandindex.txt', sep = '\t', header = TRUE, row.names = NULL)
data3 <- read.table('gse121411_ex3_withcolumnandindex.txt', sep = '\t', header = TRUE, row.names = NULL)

colnames(data1)[1]<- 'gene1'
colnames(data2)[1]<- 'gene2'
colnames(data3)[1]<- 'gene3'
# Remove duplicates based on the 'Gene' column
data1 <- data1[!duplicated(data1$gene1), ]
data2 <- data2[!duplicated(data2$gene2), ]
data3 <- data3[!duplicated(data3$gene3), ]

# Now you can assign 'Gene' column as row names and remove it from the data frame
row.names(data1) <- data1$gene1; data1$gene1 <- NULL
row.names(data2) <- data2$gene2; data2$gene2 <- NULL
row.names(data3) <- data3$gene3; data3$gene3 <- NULL

# Merge the three data frames by row names (gene identifiers)
data <- cbind(data1, data2, data3)

# Subsetting the data based on control and treatment conditions
control_cols <- grep("GFP\\.Fulv", colnames(data))
exclude_control_cols <- grepl("GFP\\.Fulv\\+Palbo|GFP\\.Fulv\\+Ner", colnames(data))

# Keep only columns that contain "GFP.Fulv" and do not contain "GFP.Fulv+Palbo" or "GFP.Fulv+Ner"
control_cols <- control_cols[!exclude_control_cols[control_cols]]

experiment_cols <- grep("S653C\\.Fulv|L755S\\.Fulv|V777L\\.Fulv|L869R\\.Fulv", colnames(data))

# Exclude columns with "+Palbo" or "+Ner" in their names in the experiment columns
exclude_experiment_cols <- grepl("\\+Palbo|\\+Ner", colnames(data))

# Keep only columns that do not contain "+Palbo" or "+Ner"
experiment_cols <- experiment_cols[!exclude_experiment_cols[experiment_cols]]

control_data <- data[,control_cols]
experiment_data <- data[,experiment_cols]

# Combine the control and experimental data
combined_data <- cbind(control_data, experiment_data)

# Create a vector of group labels
group <- factor(c(rep("control", length(control_cols)), rep("experiment", length(experiment_cols))))

# Create a design matrix
design <- model.matrix(~0+group)

# Here we assume that "batch" is a vector that represents the batch for each sample
# batch <- c(rep("batch1", length_of_batch1), rep("batch2", length_of_batch2), ...)

# Perform voom normalization
v <- voom(combined_data, design, plot=TRUE)

# Remove batch effect
# v <- removeBatchEffect(v, batch=batch)

# Fit the model
fit <- lmFit(v, design)

contrast.matrix <- makeContrasts(groupexperiment-groupcontrol, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Extract the results
results <- topTable(fit2, adjust="fdr", number=Inf)

# Print the results
print(results)

# Write results to a text file
write.table(results, file = "GFP_fulv_DEG_results.txt", sep = "\t", quote = FALSE, row.names = TRUE)
