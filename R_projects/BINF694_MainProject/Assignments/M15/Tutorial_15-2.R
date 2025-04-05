# *************************************
# *  A Manning Smith
# *  04/05/2025
# *  Tutorial 15.2
# *
# *
# *
# *************************************

# Load Libraries
#library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(edgeR)


## Loading and Preparing the Dataset
cancer_counts <- read.table(file="oral_carcinoma_counts.txt", sep = "\t", header = T)
head(cancer_counts)
y <- DGEList(counts=cancer_counts[,3:8], genes=cancer_counts[,1:2])


## Editing Row Labels
rownames(y$counts) <- rownames(y$genes) <- y$genes$EntrezGene
y$genes$EntrezGene <- NULL
y$samples


## Calculating Library Scaling Factors
y <- calcNormFactors(y)
norm.factor


## Organizing Samples into Groups
y$samples$group = c("N", "T", "N", "T", "N", "T")
## Creating a Multi-Dimensional Scaling (MDS) Plot
plotMDS(y)


## Differential Expression Analysis in edgeR Using the Exact Test
### A. Calculating the Dispersion
y <- estimateDisp(y)

### B. Creating a BCV Plot
plotBCV(y)

### C. Differential Expression Analysis
et <- exactTest(y, pair=c("N","T"))
summary(de <- decideTests(et))

### D. Visualizing Differentially Expressed Genes
detags <- rownames(y)[as.logical(de)]
plotSmear(et, de.tags=detags)
abline(h=c(-1, 1), col="blue")

### E. Viewing a List of Differentially Expressed Genes
diffExpGenes <- topTags(et, n=6, p.value = 0.05)
diffExpGenes$table

### F. Exporting Differential Expression Results
write.table(diffExpGenes$table, file="tumor_v_normal_exactTest.txt", sep = "\t", row.names=TRUE, col.names=NA)


## Differential Expression Analysis in EdgeR Using the Generalized Linear Model (GLM)
### A. Constructing the Design Matrix
Patient <- factor(c(8,8,33,33,51,51))
Tissue <- factor(c("N","T","N","T","N","T"))
design <- model.matrix(~Patient+Tissue)
rownames(design) <- colnames(y)
design

### B. Calculating the Dispersion
y2 <- y
y2 <- estimateDisp(y2, design, robust=TRUE)