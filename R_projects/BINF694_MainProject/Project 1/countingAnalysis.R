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

## Gene Heatmap
top_genes <- rownames(res_ordered)[1:50]  # Top 50 DE genes
top_gene_counts <- assay(vsd)[top_genes, ]
pheatmap(top_gene_counts,
         cluster_rows=TRUE,
         cluster_cols=TRUE,
         show_rownames=TRUE,
         annotation_col=data.frame(Tissue=sample_info$tissue, row.names=rownames(sample_info)),
         main="Top 50 Differentially Expressed Genes")

# Count significant genes (padj < 0.05 and |log2FC| > 1)
sig_genes <- sum(res$padj < 0.05 & abs(res$log2FoldChange) > 1, na.rm=TRUE)

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