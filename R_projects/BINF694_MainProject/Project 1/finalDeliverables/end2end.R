# *************************************
# *  A Manning Smith
# *  05/16/2025
# *  End 2 End Code Analysis
# *
# *
# *************************************

# Load Libraries
library(DESeq2) # Version 1.48.0
library(pheatmap) # Version 1.0.12
library(RColorBrewer) # Version 1.1.3
library(ggplot2) # Version 3.5.2
library(readxl) # Version 1.4.5
library(tibble) # Version 3.2.1
library(edgeR) # Version 4.6.1
library(dplyr) # Version 1.1.4
library(TissueEnrich) # Version 1.28.0
library(tidyverse) # Version 2.0.0
library(biomaRt) # Version 2.64.0
library(writexl) # Version 1.5.4
library(knitr) # Version 1.50
library(png) # Version 0.1-8
library(tinytex) # Version 0.57
library(imager) # Version 1.0.3
library(bookdown) # Version 0.43

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

## Load Enrichment Results
enrichment_results <- read_excel("Project 1/finalDeliverables/results/enrichmentResults.xlsx", sheet = "ER_results")

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

vsd <- vst(dds, blind=FALSE)

## Sort results by adjusted p-value
res_ordered <- res[order(res$padj), ]
summary(res)

## Sort Genes by significance
sig_genes <- res_ordered[!is.na(res_ordered$padj) & res_ordered$padj < 0.05,] # Get differentially expressed genes with adjusted p-value < 0.05
liver_genes <- sig_genes[sig_genes$log2FoldChange > 0,] # Separate liver and stomach upregulated genes
stomach_genes <- sig_genes[sig_genes$log2FoldChange < 0,]

## Sort Top 50 Genes
gene_list <- rownames(res_ordered)
top_50_genes <- rownames(res_ordered)[1:50]  # Top 50 DE genes
top_gene_counts <- assay(vsd)[top_50_genes, ]

## Save Top 50 Genes to txt
top_50_genes_df <- data.frame(Gene=top_50_genes, row.names=NULL)
write.table(top_50_genes_df, file="Project 1/finalDeliverables/results/top_50_genes.txt", sep="\t", row.names=FALSE, quote=FALSE)

# Separate liver and stomach upregulated genes
liver_genes <- sig_genes[sig_genes$log2FoldChange > 0,]
stomach_genes <- sig_genes[sig_genes$log2FoldChange < 0,]

# Get top 250 from each tissue (sorted by padj)
top_250_liver <- head(liver_genes, 250)
top_250_stomach <- head(stomach_genes, 250)
top_400_total <- head(sig_genes, 400)

# Create files for DAVID
write.table(rownames(top_250_liver), "Project 1/results/top_250_liver_DEGs.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(rownames(top_250_stomach), "Project 1/results/top_250_stomach_DEGs.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(rownames(top_400_total), "Project 1/results/top_400_total_DEGs.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

# Tables

## Table 1
trim_df <- read.csv("~/Documents/binf694/R_projects/BINF694_MainProject/Project 1/data/trimming_results.csv", header=TRUE)
kable(trim_df, caption = "Trimming Results | Table 1")

## Table 2
map_df <- read.csv("~/Documents/binf694/R_projects/BINF694_MainProject/Project 1/data/mapping_results.csv", header=TRUE)
kable(map_df, caption = "Mapping Results | Table 2")

# Figures

## Figure 1 - PCA Plot
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup="tissue") + 
  theme_minimal() +
  ggtitle("PCA of Liver vs Stomach Samples | Figure 1")

## Figure 2 - Sample Distance Matrix
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,
         display_numbers=TRUE,
         main="Sample Distance Matrix | Figure 2")

## Figure 3 - Top 50 Gene Heatmap
pheatmap(top_gene_counts,
         cluster_rows=TRUE,
         cluster_cols=TRUE,
         show_rownames=TRUE,
         annotation_col=data.frame(Tissue=sample_info$tissue, row.names=rownames(sample_info)),
         main="Top 50 Differentially Expressed Genes | Figure 3")

## Figure 4 - Heatmap of Liver and Stomach Genes
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
#rownames(top_gene_counts) <- selected_genes[rownames(top_gene_counts), "gene"]

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
         main="Differentially Expressed Genes by Tissue | Figure 4",
         fontsize_row=8)

## Figure 5 - Heatmap of Selected Genes
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
         annotation_col=data.frame(Tissue=sample_info$tissue, row.names=rownames(sample_info)),
         main = "Gene Expression by Selected Genes | Figure 5",
)

## Figure 6 - Volcano Plot
# Create a bar plot for the top 10 enriched biological processes
enrich_bp <- enrichment_results %>%
  filter(category == "GOTERM_BP_DIRECT") %>%
  filter(run_file == "top400_total_DEG-DAVID") %>%
  arrange(desc(Count)) %>%
  head(10)

# Bar Chart
ggplot(enrich_bp, aes(x = reorder(Term, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 10 Enriched Biological Processes | Figure 6",
       x = "Biological Process",
       y = "Gene Count") +
  theme_minimal()

## Figure 7 - Cellular Component
# Create a bar plot for the top 10 enriched biological processes
enrich_bp <- enrichment_results %>%
  filter(category == "GOTERM_CC_DIRECT") %>%
  filter(run_file == "top400_total_DEG-DAVID") %>%
  arrange(desc(Count)) %>%
  head(10)

# Bar Chart
ggplot(enrich_bp, aes(x = reorder(Term, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 10 Enriched Cellular Component | Figure 7",
       x = "Cellular Component",
       y = "Gene Count") +
  theme_minimal()

## Figure 8 - Molecular Function
# Create a bar plot for the top 10 enriched Molecular functions
enrich_bp <- enrichment_results %>%
  filter(category == "GOTERM_MF_DIRECT") %>%
  filter(run_file == "top400_total_DEG-DAVID") %>%
  arrange(desc(Count)) %>%
  head(10)

# Bar Chart
ggplot(enrich_bp, aes(x = reorder(Term, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 10 Enriched Molecular Function | Figure 8",
       x = "Molecular Function",
       y = "Gene Count") +
  theme_minimal()

## Figure 9 - KEGG Pathway
# Create a bar plot for the top 10 enriched KEGG pathways
enrich_kegg <- enrichment_results %>%
  filter(category == "KEGG_PATHWAY") %>%
  filter(run_file == "top400_total_DEG-DAVID") %>%
  arrange(desc(Count)) %>%
  head(10)

# Bar Chart
ggplot(enrich_kegg, aes(x = reorder(Term, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 10 Enriched KEGG Pathways | Figure 9",
       x = "KEGG Pathway",
       y = "Gene Count") +
  theme_minimal()

## Figure 10 - Interpro
# Create a bar plot for the top 10 enriched Interpro
enrich_interpro <- enrichment_results %>%
  filter(category == "INTERPRO") %>%
  filter(run_file == "top400_total_DEG-DAVID") %>%
  arrange(desc(Count)) %>%
  head(10)
# Bar Chart
ggplot(enrich_interpro, aes(x = reorder(Term, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 10 Enriched Interpro | Figure 10",
       x = "Interpro",
       y = "Gene Count") +
  theme_minimal()

## Figure 11 - String Diagram
img <- load.image("Project 1/finalDeliverables/results/string_normal_image400.png")
plot(img, axes = FALSE, main = "STRING Diagram | Figure 11")

## Figure 12 - Tissue Results
img <- load.image("Project 1/finalDeliverables/results/tissueResults.png")
plot(img, axes = FALSE, main = "Tissue Results | Figure 12")

## Figure 13 - Cluster 1
img <- load.image("Project 1/finalDeliverables/results/Cluster1.png")
plot(img, axes = FALSE, main = "Cluster 1 | Figure 13")

## Figure 14 - Cluster 2
img <- load.image("Project 1/finalDeliverables/results/Cluster2.png")
plot(img, axes = FALSE, main = "Cluster 2 | Figure 14")