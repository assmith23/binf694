# *************************************
# *  A Manning Smith
# *  04/05/2025
# *  Tutorial 14.3
# *
# *
# *
# *************************************

# Load Libraries
library(DESeq2)
library(pheatmap)
library(RColorBrewer)


## Viewing and Loading the Datasets
  # Viewed in csv


## Load Data
exp_counts <- read.table(file="Assignments/M14/mouse_mammary_counts2.txt", sep = "\t", header = T, row.names=1)
sample_data <- read.table(file="Assignments/M14/sample_data.txt", sep = "\t", header = T, row.names=1)
head(exp_counts)
sample_data


## Creating a DESeqDataSet from the Input Data
# Q3
data_deseq <- DESeqDataSetFromMatrix(countData = exp_counts, colData = sample_data, design = ~ 1)
head(counts(data_deseq))

## Filtering the Dataset
# Q4
nrow(data_deseq)

# Q5
data_deseq <- data_deseq[ rowSums(counts(data_deseq)) > 1, ]
nrow(data_deseq)


## Transforming the Dataset
rld <- rlog(data_deseq, blind=FALSE)


## Creating a Sample Distance Heatmap
# Q6 --> 0

sampleDists <- dist(t(assay(rld)))
sampleDists


## Construct the heatmap
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste( rld$cell_type, rld$dev_stage, rld$replicate, sep="-" )
colnames(sampleDistMatrix) <- paste( rld$cell_type, rld$dev_stage, rld$replicate, sep="-" )
colors <- colorRampPalette( rev(brewer.pal(9, "OrRd")) )(255)
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists,
         col=custom_colors, display_numbers=TRUE, breaks = custom_breaks)


### Get mean of each Cell Type
lum_idx <- 7:12
basal_idx <- 1:6

within_lum <- sampleDistMatrix[lum_idx, lum_idx]
within_basal <- sampleDistMatrix[basal_idx, basal_idx]

get_upper <- function(mat) mat[upper.tri(mat)]

mean_within_lum <- mean(get_upper(within_lum))
mean_within_basal <- mean(get_upper(within_basal))


## Principal Component Analysis (PCA)
plotPCA(rld, intgroup = c("cell_type", "dev_stage"))


## Gene Heatmap
geneVars <- rowVars(assay(rld))
geneVarsOrdered <- order(geneVars, decreasing = TRUE)
topVarGenes <- head(geneVarsOrdered, 50)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c("cell_type","dev_stage")])
clear_col_names <- paste( rld$cell_type, rld$dev_stage, rld$replicate, sep="." )
topGenesHeatmap <- pheatmap(mat, annotation_col=df, labels_col = clear_col_names)