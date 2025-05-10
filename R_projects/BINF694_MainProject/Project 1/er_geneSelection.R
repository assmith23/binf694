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
library(writexl)

# Functions
get_summary_stats <- function(data) {
  result <- data.frame(
    Metric = c("Mean", "Median", "Minimum", "Maximum", "Q1 (25%)", "Q3 (75%)", "Std Dev", "Count"),
    row.names = "Metric"
  )
}

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

vsd <- vst(dds, blind=FALSE)

# Get differentially expressed genes with adjusted p-value < 0.05
sig_genes <- res_ordered[!is.na(res_ordered$padj) & res_ordered$padj < 0.05,]

# Separate liver and stomach upregulated genes
liver_genes <- sig_genes[sig_genes$log2FoldChange > 0,]
stomach_genes <- sig_genes[sig_genes$log2FoldChange < 0,]

# Get top 250 from each tissue (sorted by padj)
top_250_liver <- head(liver_genes, 250)
top_250_stomach <- head(stomach_genes, 250)

# Create files for DAVID
write.table(rownames(top_250_liver), "Project 1/results/top_250_liver_DEGs.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

write.table(rownames(top_250_stomach), "Project 1/results/top_250_stomach_DEGs.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

# Also save complete results table as requested in the project
complete_results <- as.data.frame(res_ordered)
complete_results$gene <- rownames(complete_results)
complete_results$tissue <- ifelse(complete_results$log2FoldChange > 0, "liver", "stomach")
complete_results <- complete_results[, c("gene", "tissue", "log2FoldChange", "padj", "pvalue", "baseMean", "lfcSE", "stat")]

# Filter top 250 from each tissue
liver_top <- complete_results[complete_results$tissue == "liver",][1:250,]
stomach_top <- complete_results[complete_results$tissue == "stomach",][1:250,]
top_500 <- rbind(liver_top, stomach_top)

# Export both sheets to Excel
write_xlsx(list("All_Results" = complete_results, 
                "Top_500" = top_500), 
           "Project 1/results/liver_stomach_DEG_results.xlsx")


liver_top_genes <- rownames(liver_genes)[1:25]
stomach_top_genes <- rownames(stomach_genes)[1:25]
top_50_genes <- c(liver_top_genes, stomach_top_genes)

# Get normalized counts for these genes
top_gene_counts <- assay(vsd)[top_50_genes, ]

# Create annotation data frame
anno_col <- data.frame(Tissue = sample_info$tissue, row.names = rownames(sample_info))

# Create color palette
tissue_colors <- c("liver" = "#1B9E77", "stomach" = "#D95F02")
anno_colors <- list(Tissue = tissue_colors)

# Generate heatmap
pheatmap(top_gene_counts,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         fontsize_row = 8,
         annotation_col = anno_col,
         annotation_colors = anno_colors,
         main = "Top 50 Differentially Expressed Genes",
         scale = "row")  # Scale rows for better visualization