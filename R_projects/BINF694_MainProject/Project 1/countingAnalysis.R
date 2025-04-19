# *************************************
# *  A Manning Smith
# *  04/12/2025
# *  Interim Analysis Working Doc
# *
# *
# *
# *************************************

# Load Libraries
library(DESeq2)
library(pheatmap)
library(RColorBrewer)


## Load Data
exp_counts <- read.table(file="", sep = "\t", header = T, row.names=1)
sample_data <- read.table(file="sample_data.txt", sep = "\t", header = T, row.names=1)
head(exp_counts)
sample_data