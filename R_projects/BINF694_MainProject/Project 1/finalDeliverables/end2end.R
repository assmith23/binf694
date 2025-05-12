# *************************************
# *  A Manning Smith
# *  05/16/2025
# *  End 2 End Code Analysis
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

# Load data

## Set Working Directory
setwd("~/Documents/binf694/R_projects/BINF694_MainProject")

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
gene_list <- gene_info$gene

# End 2 End

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

## Run DESeq2 analysis
dds <- DESeq(dds)

res <- results(dds, contrast=c("tissue", "liver", "stomach"))

## Sort results by adjusted p-value
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

## Save Top 50 Genes to txt
top_50_genes_df <- data.frame(Gene=top_50_genes, row.names=NULL)
write.table(top_50_genes_df, file="Project 1/finalDeliverables/results/top_50_genes.txt", sep="\t", row.names=FALSE, quote=FALSE)