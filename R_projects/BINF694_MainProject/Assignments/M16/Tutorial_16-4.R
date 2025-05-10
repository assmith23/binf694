# *************************************
# *  A Manning Smith
# *  04/05/2025
# *  Tutorial 16.4
# *
# *
# *
# *************************************

# Load Libraries
#library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(edgeR)
library(GO.db)
library(org.Hs.eg.db)

## Loading and Preparing the Dataset
cancer_counts <- read.table(file="M15/oral_carcinoma_counts.txt", sep = "\t", header = T)
head(cancer_counts)
y <- DGEList(counts=cancer_counts[,3:8], genes=cancer_counts[,1:2])
rownames(y$counts) <- rownames(y$genes) <- y$genes$EntrezGene
y$genes$EntrezGene <- NULL
y <- calcNormFactors(y)
y$samples$group = c("N", "T", "N", "T", "N", "T")
y <- estimateDisp(y)
et <- exactTest(y, pair=c("N","T"))

## GO Biological Process Enrichment Analysis Using goana
goExact <- goana(et)
goExact.sort_values <- topGO(goExact, ont="BP", sort="Up", n=30)
#topGO(goExact, ont="BP", n=30)

top_5_upregulated = goExact[order(goExact$P.Up), ][1:5, ]
top_5_downregulated = goExact[order(goExact$P.Down), ][1:5, ]

top_5_upregulated
top_5_downregulated
