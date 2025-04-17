# *************************************
# *  A Manning Smith
# *  04/09/2025
# *  Tutorial 16.3
# *
# *  Supplemental code
# *
# *************************************


# Load Libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(reshape2)


geneData <- read_excel("Assignments/M16/tnbc_vs_normal_up_rnaseq.xlsx")
View(geneData)

summary(geneData)

subset_temp <- geneData %>%
  filter(logFC >= 1)
summary(subset_temp)

ggplot(subset_temp, aes(x = factor(1), y = logFC)) +
  geom_boxplot(outlier.shape = NA) +  # hides individual points beyond whiskers
  labs(x = NULL, y = "logFC") +
  scale_x_discrete(breaks = NULL) +
  scale_y_continuous(limits = c(1, 4), breaks = seq(0, 5, by = 1)) +
  labs(title = "LogFC >= 1 | Boxplot") +
  theme_minimal()


cor_matrix <- cor(subset_temp[, c("logFC", "FDR", "PValue")], use = "complete.obs")
cor_df <- melt(cor_matrix)

ggplot(cor_df, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1,1), space = "Lab", name = "Correlation") +
  geom_text(aes(label = round(value, 2)), color = "black", size = 4) +
  theme_minimal() +
  labs(title = "Correlation Matrix", x = "", y = "") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))


# Spec 1: Larger Sample
subset1 <- geneData %>%
  filter(logFC >= 1 & logFC <= 2.48) %>%
  filter(FDR <= 0.01)
summary(subset1)
write.table(subset1$GeneSymbol, "Assignments/M16/subset1_gene.txt", sep = "\t", row.names = FALSE, quote = FALSE)


# Spec 2: Smaller Sample
subset2 <- geneData %>%
  filter(logFC >= 2.48) %>%
  filter(FDR <= 0.01)
summary(subset3)
write.table(subset3$GeneSymbol, "Assignments/M16/subset2_gene.txt", sep = "\t", row.names = FALSE, quote = FALSE)