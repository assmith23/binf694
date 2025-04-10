# *************************************
# *  A Manning Smith
# *  04/09/2025
# *  Tutorial 16.3
# *
# *
# *
# *************************************


# Load Libraries
library(readxl)
library(dplyr)


geneData <- read_excel("Assignments/M16/tnbc_vs_normal_up_rnaseq.xlsx")
View(geneData)

summary(geneData)

subset_temp <- geneData %>%
  filter(PValue >= 0.013 & PValue <= 0.481)
summary(subset_temp)

cor(subset2[, c("logFC", "FDR", "PValue")], use = "complete.obs")


# Specs: LogFC | 0.890-1.077 -- FDR | 0.056-0.668
subset1 <- geneData %>%
  filter(logFC >= 0.890 & logFC <= 1.077) %>%
  filter(FDR >= 0.056 & FDR <= 0.668)
summary(subset1)


# Specs: LogFC | 0.890-1.077 -- FDR | 0.056-0.668
subset2 <- geneData %>%
  filter(logFC >= 0.662 & logFC <= 0.788) %>%
  filter(FDR >= 0.156 & FDR <= 0.329)
summary(subset2)
