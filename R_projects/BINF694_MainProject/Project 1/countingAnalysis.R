# *************************************
# *  A Manning Smith
# *  04/12/2025
# *  Interim Analysis Working Doc
# *  Working end-2-end
# *
# *
# *************************************

# Load Libraries
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(readxl)
library(tibble)
library(edgeR)
library(dplyr)
library(TissueEnrich)
library(tidyverse)
library(biomaRt)

# Functions
get_summary_stats <- function(data) {
  result <- data.frame(
    Metric = c("Mean", "Median", "Minimum", "Maximum", "Q1 (25%)", "Q3 (75%)", "Std Dev", "Count"),
    row.names = "Metric"
)
}

# Load Data


## Load Count Data
liver_a_counts <- read.table(file="Project 1/data/liver_a.s.count.txt", sep = "\t", header = T, row.names=1)
liver_c_counts <- read.table(file="Project 1/data/liver_c.s.count.txt", sep = "\t", header = T, row.names=1)
liver_d_counts <- read.table(file="Project 1/data/liver_d.s.count.txt", sep = "\t", header = T, row.names=1)
stomach_a_counts <- read.table(file="Project 1/data/stomach_a.s.count.txt", sep = "\t", header = T, row.names=1)
stomach_3b_counts <- read.table(file="Project 1/data/stomach_3b.s.count.txt", sep = "\t", header = T, row.names=1)
stomach_3a_counts <- read.table(file="Project 1/data/stomach_3a.s.count.txt", sep = "\t", header = T, row.names=1)
names(liver_a_counts)[1] <- "count"
names(liver_c_counts)[1] <- "count"
names(liver_d_counts)[1] <- "count"
names(stomach_a_counts)[1] <- "count"
names(stomach_3b_counts)[1] <- "count"
names(stomach_3a_counts)[1] <- "count"


## Load Sample Info
sample_info <- read_excel("Project 1/data/masterData.xlsx", sheet = "sample_info")
sample_info <- sample_info[, 1:2] 
sample_info <- column_to_rownames(sample_info, var = colnames(sample_info)[1])

## Load Gene Info
gene_info <- read_excel("Project 1/data/masterData.xlsx", sheet = "gene_info")
gene_info <- dplyr::filter(gene_info, filter_selection == "YES")
gene_list <- gene_info$gene ## Combine Gene Lists

## Create count matrix
count_matrix <- cbind(
  liver_a = liver_a_counts$count,liver_c = liver_c_counts$count,liver_d = liver_d_counts$count,
  stomach_a = stomach_a_counts$count,stomach_3a = stomach_3a_counts$count,stomach_3b = stomach_3b_counts$count
)
rownames(count_matrix) <- rownames(liver_a_counts)

## Creating a DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_info,
  design = ~ tissue
)
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Run DESeq2 analysis
dds <- DESeq(dds)

res <- results(dds, contrast=c("tissue", "liver", "stomach"))
# Sort results by adjusted p-value
res_ordered <- res[order(res$padj), ]
summary(res)

## PCA Plot
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup="tissue") + 
  theme_minimal() +
  ggtitle("PCA of Liver vs Stomach Samples")


## Heatmap
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,
         display_numbers=TRUE,
         main="Sample Distance Matrix")

## Top 50 Gene Heatmap
gene_list <- rownames(res_ordered)
top_50_genes <- rownames(res_ordered)[1:50]  # Top 50 DE genes
top_gene_counts <- assay(vsd)[top_50_genes, ]
pheatmap(top_gene_counts,
         cluster_rows=TRUE,
         cluster_cols=TRUE,
         show_rownames=TRUE,
         annotation_col=data.frame(Tissue=sample_info$tissue, row.names=rownames(sample_info)),
         main="Top 50 Differentially Expressed Genes")

## Top 25 Genes for Each Tissue
res_fil <- res[!is.na(res$padj), ]
significant_genes <- res_fil[res_fil$padj < 0.05, ]
liver_genes <- significant_genes[significant_genes$log2FoldChange > 0, ]
stomach_genes <- significant_genes[significant_genes$log2FoldChange < 0, ]
liver_genes_ordered <- liver_genes[order(liver_genes$log2FoldChange, decreasing = TRUE), ]
stomach_genes_ordered <- stomach_genes[order(stomach_genes$log2FoldChange, decreasing = FALSE), ]
liver_top <- rownames(liver_genes_ordered)[1:25]
stomach_top <- rownames(stomach_genes_ordered)[1:25]
selected_genes <- c(liver_top, stomach_top)
top_gene_counts <- assay(vsd)[selected_genes, ]
rownames(top_gene_counts) <- gene_info[rownames(top_gene_counts), "gene"]

gene_origin <- data.frame(
  Origin = ifelse(rownames(top_gene_counts) %in% liver_top, "Liver", "Stomach"),
  row.names = rownames(top_gene_counts)
)

pheatmap(top_gene_counts,
         cluster_rows=TRUE,
         cluster_cols=TRUE,
         show_rownames=TRUE,
         annotation_col=data.frame(Tissue=sample_info$tissue, row.names=rownames(sample_info)),
         annotation_row=gene_origin,
         cutree_rows = 2,
         gaps_row = 25,
         main="Top Differentially Expressed Genes (25 Liver, 25 Stomach)",
         fontsize_row=8)

gene_overlap <- intersect(selected_genes, gene_info)


# Count significant genes (padj < 0.05 and |log2FC| > 1)
sig_genes <- res[!is.na(res$padj) & res$padj < 0.01 & !is.na(res$log2FoldChange) & abs(res$log2FoldChange) > 1.5, ]
cnt_sig_genes <- sum(res$padj < 0.01 & abs(res$log2FoldChange) > 1.5, na.rm=TRUE)
cnt_sig_genes

summary(sig_genes@listData)
lapply(as.data.frame(sig_genes@listData), summary)
  

summary_df <- data.frame(
Metric = c("Mean", "Median", "Minimum", "Maximum", "Q1 (25%)", 
           "Q3 (75%)", "Std Dev", "Count")
)

# Calculate statistics for each column
for (col in names(sig_genes@listData)) {
  col_data <- sig_genes@listData[[col]]
  summary_df[[col]] <- c(
    mean(col_data, na.rm = TRUE),
    median(col_data, na.rm = TRUE),
    min(col_data, na.rm = TRUE),
    max(col_data, na.rm = TRUE),
    quantile(col_data, 0.25, na.rm = TRUE),
    quantile(col_data, 0.75, na.rm = TRUE),
    sd(col_data, na.rm = TRUE),
    sum(!is.na(col_data))
  )
}

# Set row names and display formatted table
rownames(summary_df) <- summary_df$Metric
summary_df$Metric <- NULL

#kable(summary_df, caption = "Signficant Genes | Table X")

# Get gene symbols from rownames
sig_genes$symbol <- rownames(sig_genes)

# Sort by absolute fold change
if(nrow(sig_genes) > 300) {
  ordered_genes <- sig_genes[order(abs(sig_genes$log2FoldChange), decreasing=TRUE), ]
  top_genes <- ordered_genes[1:300, ]
} else {
  top_genes <- sig_genes
}

# Convert to regular data frame with gene symbols
top_genes_df <- as.data.frame(top_genes)

# Add gene symbols as column
top_genes_df <- as.data.frame(top_genes)
top_genes_df$gene_symbol <- rownames(top_genes)

# Reorganize columns to put gene symbol first
cols <- c("gene_symbol", setdiff(names(top_genes_df), "gene_symbol"))
top_genes_df <- top_genes_df[, cols]

# Write out results with gene symbol column
write.csv(top_genes_df, "Project 1/results/top_genes_for_enrichment.csv", row.names=FALSE)

# Create text file with just gene symbols
gene_symbols <- rownames(top_genes)
writeLines(gene_symbols, "Project 1/results/gene_symbols.txt")


# Save results to file (for your project deliverable)
write.table(as.data.frame(res_ordered), file="Project 1/results/liver_vs_stomach_DESeq2_results.tsv", 
            sep="\t", quote=FALSE)

###### Compare Selected Genes
markers<- intersect(gene_list, rownames(count_matrix))
filtered_counts <- count_matrix[markers, ]
filtered_counts[is.na(filtered_counts)] <- 0
filtered_counts[is.infinite(filtered_counts)] <- 0
row_var <- apply(filtered_counts, 1, var, na.rm=TRUE)
filtered_counts <- filtered_counts[row_var > 0,]

# Create heatmap
pheatmap(filtered_counts,
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_colnames = TRUE,
         fontsize_row = 8,
         cutree_rows = 2,
         gaps_row = 32,
         annotation_col=data.frame(Tissue=sample_info$tissue, row.names=rownames(sample_info)),
         main = "Gene Expression by Selected Genes",
)

markers<- intersect(gene_list, top_genes)
filtered_top_counts <- count_matrix[markers, ]

plotMA(res, main="MA Plot: Liver vs Stomach", ylim=c(-10, 10))

#### Part 2 Stuff

y <- DGEList(counts = count_matrix, 
             group = factor(c(rep("liver", 3), rep("stomach", 3))))

# Filter low-expressed genes
keep <- filterByExpr(y)
y <- y[keep, ]
y <- calcNormFactors(y) # Calculate normalization factors
y <- estimateDisp(y) # Estimate dispersion


plotBCV(y, main="Biological Coefficient of Variation")
plotMDS(y, labels=colnames(y), col=c(rep("blue",3), rep("red",3)), 
        main="Multi-Dimensional Scaling Plot")

et <- exactTest(y)
plotSmear(et, de.tags=rownames(topTags(et, n=100)), 
          main="Smear Plot of Differential Expression")
abline(h=c(-1, 1), col="blue")

results_df <- as.data.frame(res)
colnames(results_df)[colnames(results_df) == "padj"] <- "FDR"
colnames(results_df)[colnames(results_df) == "log2FoldChange"] <- "logFC"
summary(results_df)
markers<- intersect(gene_list, rownames(results_df))
selected_results <- results_df[markers, ]
summary(selected_results)

write.csv(selected_results, "Project 1/data/results.csv")
















markers <- intersect(gene_info$gene, rownames(count_matrix))
markers <- data.frame(gene = markers)
filtered_counts <- count_matrix[markers$gene, ]
filtered_counts[is.na(filtered_counts)] <- 0
filtered_counts[is.infinite(filtered_counts)] <- 0
row_var <- apply(filtered_counts, 1, var, na.rm=TRUE)
filtered_counts <- filtered_counts[row_var > 0,]

filtered_counts <- left_join(
  transform(filtered_counts, gene = rownames(filtered_counts)),
  select(gene_info, tissue_association),
  by = "gene"
)

gene_origin <- data.frame(
  Origin = filtered_counts$tissue_association,
  row.names = rownames(filtered_counts)
)

# Create heatmap
pheatmap(filtered_counts[1:6],
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_colnames = TRUE,
         fontsize_row = 8,
         annotation_row=gene_origin,
         annotation_col=data.frame(Tissue=sample_info$tissue, row.names=rownames(sample_info)),
         main = "Gene Expression by Selected Genes",
)