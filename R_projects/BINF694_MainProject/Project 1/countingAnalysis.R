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
liver_counts <- read.table(file="Project 1/data/liver_a.s.count.txt", sep = "\t", header = T, row.names=1)
stomach_counts <- read.table(file="Project 1/data/stomach_3b.s.count.txt", sep = "\t", header = T, row.names=1)