setwd("F:/5_ Jane_100322/03_Data_071924/Project_1549_SS_HS vs LS Brain")





#   Open Metadata File
Metadata <- read.table("metadata.txt",
                       head = TRUE,
                       sep = "\t")
str(Metadata)
head(Metadata)



#   Filter Metadata for Brain Cortex
Metadata_Cortex <- subset(Metadata, condition %in% c("BH_HS", "BH_LS"))
Metadata_Cortex


#   Rename the Metadata
rownames(Metadata_Cortex) <- c("SS_HS_2569",
                               "SS_HS_2570",
                               "SS_HS_2571",
                               "SS_HS_2573",
                               "SS_LS_2561",
                               "SS_LS_2562",
                               "SS_LS_2554",
                               "SS_LS_2555",
                               "SS_LS_2575")
Metadata_Cortex


#   Make the row names
Metadata_Cortex$X <- NULL
Metadata_Cortex










#   Open Count File
Counts <- read.table("Star_counts.txt",
                     head = TRUE,
                     sep = "\t")
head(Counts)
str(Counts)



#   Clean up the COunt File ahnd Save
rownames(Counts) <- Counts[, 1]
Counts <- Counts[7:24]
str(Counts)



#   Filter Counts for Cortex and TPPU
#   library(dplyr)  #   In order to use 'select()' function
Counts_Cortex <- dplyr::select(Counts, c("X.singular_genomics.P1549.STAR_bam.P1549_18_SS_HS_2569_S18_L001Aligned.sortedByCoord.out.bam",
                                         "X.singular_genomics.P1549.STAR_bam.P1549_19_SS_HS_2570_S19_L001Aligned.sortedByCoord.out.bam",
                                         "X.singular_genomics.P1549.STAR_bam.P1549_20_SS_HS_2571_S20_L001Aligned.sortedByCoord.out.bam",
                                         "X.singular_genomics.P1549.STAR_bam.P1549_21_SS_HS_2573_S21_L001Aligned.sortedByCoord.out.bam",
                                         "X.singular_genomics.P1549.STAR_bam.P1549_22_SS_LS_2561_S22_L001Aligned.sortedByCoord.out.bam",
                                         "X.singular_genomics.P1549.STAR_bam.P1549_23_SS_LS_2562_S23_L001Aligned.sortedByCoord.out.bam",
                                         "X.singular_genomics.P1549.STAR_bam.P1549_24_SS_LS_2554_S24_L001Aligned.sortedByCoord.out.bam",
                                         "X.singular_genomics.P1549.STAR_bam.P1549_25_SS_LS_2555_S25_L001Aligned.sortedByCoord.out.bam",
                                         "X.singular_genomics.P1549.STAR_bam.P1549_26_SS_LS_2575_S26_L001Aligned.sortedByCoord.out.bam"))
str(Counts_Cortex)



#   Rename the Counts_Cortex
colnames(Counts_Cortex) <- sub(".*_(SS_[HL]S_\\d+).*", "\\1", colnames(Counts_Cortex))
head(Counts_Cortex)
str(Counts_Cortex)

write.csv(Counts_Cortex, file = "HIPPOCAMPUS_HS vs LS/Raw Data_Cortex_STAR Counts.csv")






#   Merge the Counts_Cortex and Metadata_Cortex into ONE File
all(colnames(Counts_Cortex) %in% rownames(Metadata_Cortex))   #   Check the names match
all(colnames(Counts_Cortex) == rownames(Metadata_Cortex))   #   Check the names match



#   Normalize using DESeq2
library(DESeq2)
Dataset <- DESeqDataSetFromMatrix(countData = round(Counts_Cortex),
                                  colData = Metadata_Cortex,
                                  design = ~condition)
Dataset



#   Remove Low Count Reads
Clean <- rowSums(counts(Dataset)) >= 10   #   Form Row and Column Sums and Means (Form row and column sums and means for numeric arrays (or data frames))
Dataset <- Dataset[Clean, ]
Dataset



#   Setting reference for DEG Analysis
Dataset$Condition <- relevel(Dataset$condition, ref = "BH_LS")
DEG <- DESeq(Dataset)
Result <- results(DEG)
summary(Result)



#   Adjust threshold using alpha
Result <- results(DEG, alpha = 0.05)
summary(Result)



#   Annotation
library(org.Rn.eg.db)
Result.df <-as.data.frame(Result)
str(Result.df)
Result.df$Symbol <- mapIds(org.Rn.eg.db, #    Get the mapped ids (column) for a set of keys that are of a particular keytype. Usually returned as a named character vector.
                           rownames(Result.df),
                           keytype = "ENSEMBL", #   Discover which keytypes can be passed in to select or keys and the keytype argument
                           column = "SYMBOL")
str(Result.df)

write.csv(Result.df, "Hippocampus_HS vs LS/DEG Result_HS vs LS.csv")


###############################################################################################################

#   Volcano Plot with FoldChange    ###############################################################################



#   Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)
Result.df$diffexpressed <- "NO"



#   If log2FoldChange > 0.5 and pvalue 0.05, set as "UP"
Result.df$diffexpressed[Result.df$log2FoldChange > 0.5 & Result.df$pvalue < 0.05] <- "UP"



#   If log2FoldChange < -0.5 and pvalue < 0.05, set as"DOWN"
Result.df$diffexpressed[Result.df$log2FoldChange < -0.5 & Result.df$pvalue < 0.05] <- "DOWN"




#   All Top Genes with larger Fold Change
library(tidyverse)

Best_FoldChange <- Result.df %>%
  filter(pvalue < 0.05) %>%  # Filter significant p-values
  filter(log2FoldChange < -0.5 | log2FoldChange > 0.5) %>%  # Filter based on log2FoldChange threshold
  arrange(desc(abs(log2FoldChange))) %>%  # Sort by absolute log2FoldChange
  mutate(diffexpressed = ifelse(diffexpressed == "NO", NA, diffexpressed)) %>%  # Replace "NO" with NA in diffexpressed
  mutate(Symbol = ifelse(grepl("^LOC", Symbol), NA, Symbol)) %>%  # Replace Symbols starting with "LOC" with NA
  drop_na(diffexpressed) %>%  # Drop rows where diffexpressed is NA
  drop_na(Symbol)  # Drop rows where Symbol is NA


head(Best_FoldChange)
str(Best_FoldChange)


Best_FoldChange30 <- Result.df %>%
  filter(pvalue < 0.05) %>%
  filter(log2FoldChange < -0.5 | log2FoldChange > 0.5) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  mutate(diffexpressed = ifelse(diffexpressed == "NO", NA, diffexpressed)) %>%    # Replace "NO" with NA
  mutate(Symbol = ifelse(grepl("^LOC", Symbol), NA, Symbol)) %>%  # Replace Symbols starting with "LOC" with NA
  drop_na(diffexpressed) %>%      # Drop rows with NA in the diffexpressed column
  drop_na(Symbol) %>%
  head(30)

str(Best_FoldChange30)


#   Make Volcano column with all best genes
# Result.df$volcano <- ifelse(Result.df$Symbol %in% Best_FoldChange$Symbol, Result.df$Symbol, NA)


#   Make volcano column with 30 best genes
Result.df$volcano <- ifelse(Result.df$Symbol %in% Best_FoldChange30$Symbol, Result.df$Symbol, NA)





#Create a basic Volcano Plot
library(ggplot2)
ggplot(data = Result.df,
       aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point()




#   Final Volcano Plot
library(ggrepel)  #   In order to use 'geom_text_repel()" function
ggplot(data = Result.df,
       aes(x = log2FoldChange,
           y = -log10(pvalue),
           col = diffexpressed,
           label = volcano)) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "red", linetype = 'dashed') + #Add Threshold Lines
  geom_hline(yintercept = -log10(0.05), col = "red", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(name = "Expression",
                     values = c("#00AFBB", "#525251", "#f55142"), #   to set the colors of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  coord_cartesian(ylim = c(0, 15), xlim = c(-8, 8)) + #   Even if we have genes with p-values or log2FC outside the limits, they will still be plotted
  labs(title = "Hippocampus SS LS vs HS",
       x = expression("log"[2]*"FC"),
       y = expression("-log"[10]*"pvalue")) +
  theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5, vjust = 1.5),
        axis.title = element_text(size = 20, face = "bold"),
        axis.ticks = element_line(size = 1),
        axis.text = element_text(color = "black", size = 15),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.position = c(0.85, .85),
        axis.line = element_line(color = "black", size = 1.3),
        panel.background = element_blank()) +
  scale_x_continuous(breaks = seq(-8, 8, 2)) + #    Customize the breaks in the x axis
  geom_text_repel(max.overlaps = Inf)



###############################################################################################################

#   Volcano Plot with pvalue    ###############################################################################


Best_Pvalue <- Result.df %>%
  filter(pvalue < 0.05) %>%
  filter(log2FoldChange < -0.5 | log2FoldChange > 0.5) %>%
  arrange(pvalue) %>%
  mutate(diffexpressed = ifelse(diffexpressed == "NO", NA, diffexpressed)) %>%    # Replace "NO" with NA
  mutate(Symbol = ifelse(grepl("^LOC", Symbol), NA, Symbol)) %>%  # Replace Symbols starting with "LOC" with NA
  drop_na(diffexpressed) %>%      # Drop rows with NA in the diffexpressed column
  drop_na(Symbol)

Best_Pvalue30 <- Result.df %>%
  filter(pvalue < 0.05) %>%
  filter(log2FoldChange < -0.5 | log2FoldChange > 0.5) %>%
  arrange(pvalue) %>%
  mutate(diffexpressed = ifelse(diffexpressed == "NO", NA, diffexpressed)) %>%    # Replace "NO" with NA
  mutate(Symbol = ifelse(grepl("^LOC", Symbol), NA, Symbol)) %>%  # Replace Symbols starting with "LOC" with NA
  drop_na(diffexpressed) %>%      # Drop rows with NA in the diffexpressed column
  drop_na(Symbol) %>%
  head(30)


str(Best_Pvalue30)


Result.df$delabel <- ifelse(Result.df$Symbol %in% Best_Pvalue30$Symbol, Result.df$Symbol, NA)
Result.df$delabel <- ifelse(grepl("^LOC", Result.df$delabel), NA, Result.df$delabel)





#   Final Volcano Plot
library(ggrepel)  #   In order to use 'geom_text_repel()" function
ggplot(data = Result.df,
       aes(x = log2FoldChange,
           y = -log10(pvalue),
           col = diffexpressed,
           label = delabel)) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "red", linetype = 'dashed') + #Add Threshold Lines
  geom_hline(yintercept = -log10(0.05), col = "red", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(name = "Expression",
                     values = c("#00AFBB", "#525251", "#f55142"), #   to set the colors of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  coord_cartesian(ylim = c(0, 15), xlim = c(-8, 8)) + #   Even if we have genes with p-values or log2FC outside the limits, they will still be plotted
  labs(title = "Hippocampus SS LS vs HS",
       x = expression("log"[2]*"FC"),
       y = expression("-log"[10]*"pvalue")) +
  theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5, vjust = 1.5),
        axis.title = element_text(size = 20, face = "bold"),
        axis.ticks = element_line(size = 1),
        axis.text = element_text(color = "black", size = 15),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.position = c(0.85, .85),
        axis.line = element_line(color = "black", size = 1.3),
        panel.background = element_blank()) +
  scale_x_continuous(breaks = seq(-8, 8, 2)) + #    Customize the breaks in the x axis
  geom_text_repel(max.overlaps = Inf)



###############################################################################################################

#   Volcano Plot with BOTH Foldchange & pvalue    #############################################################


# Combine the conditions into a single column
Result.df$category <- ifelse(
  Result.df$Symbol %in% Best_Pvalue30$Symbol & Result.df$diffexpressed == "DOWN", "Pvalue_DOWN",
  ifelse(Result.df$Symbol %in% Best_Pvalue30$Symbol & Result.df$diffexpressed == "UP", "Pvalue_UP",
         ifelse(Result.df$Symbol %in% Best_FoldChange30$Symbol & Result.df$diffexpressed == "DOWN", "FoldChange_DOWN",
                ifelse(Result.df$Symbol %in% Best_FoldChange30$Symbol & Result.df$diffexpressed == "UP", "FoldChange_UP", Result.df$diffexpressed))))



# Add a single label column for labeling the graph
Result.df$label <- ifelse(
  Result.df$Symbol %in% Best_Pvalue30$Symbol & Result.df$diffexpressed == "DOWN" |
    Result.df$Symbol %in% Best_Pvalue30$Symbol & Result.df$diffexpressed == "UP" |
    Result.df$Symbol %in% Best_FoldChange30$Symbol & Result.df$diffexpressed == "DOWN" |
    Result.df$Symbol %in% Best_FoldChange30$Symbol & Result.df$diffexpressed == "UP",
  Result.df$Symbol, 
  NA)




#   Final Volcano Plot
library(ggrepel)  #   In order to use 'geom_text_repel()" function


ggplot(data = Result.df,
       aes(x = log2FoldChange,
           y = -log10(pvalue),
           col = category,
           label = label)) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "red", linetype = 'dashed') + #Add Threshold Lines
  geom_hline(yintercept = -log10(0.05), col = "red", linetype = 'dashed') +
  geom_point(data = subset(Result.df, category %in% c("Pvalue_DOWN", "Pvalue_UP", "FoldChange_DOWN", "FoldChange_UP", "Both")),
             size = 2) +  # Points included in the legend
  geom_point(data = subset(Result.df, category %in% c("UP", "DOWN", "NO")),
             size = 2, show.legend = FALSE) + # Points excluded from the legend
  scale_color_manual(name = "Expression",
                     values = c("Pvalue_DOWN" = "#00A9FF",
                                "Pvalue_UP" = "firebrick2", 
                                "FoldChange_DOWN" = "steelblue3", 
                                "FoldChange_UP" = "#FF68A1",
                                "Both" = "green3",
                                "UP" = "#f55142",
                                "DOWN" = "#00AFBB",
                                "NO" = "#525251"),
                     breaks = c("Pvalue_DOWN",
                                "Pvalue_UP",
                                "FoldChange_DOWN",
                                "FoldChange_UP",
                                "Both"),
                     labels = c("Pvalue DOWN",
                                "Pvalue UP",
                                "FoldChange DOWN",
                                "FoldChange UP",
                                "Both")) +
  coord_cartesian(ylim = c(0, 15), xlim = c(-8, 8)) + #   Even if we have genes with p-values or log2FC outside the limits, they will still be plotted
  labs(title = "Hippocampus SS LS vs HS",
       x = expression("log"[2]*"FC"),
       y = expression("-log"[10]*"pvalue")) +
  theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5, vjust = 1.5),
        axis.title = element_text(size = 20, face = "bold"),
        axis.ticks = element_line(size = 1),
        axis.text = element_text(color = "black", size = 15),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.position = c(0.85, .85),
        axis.line = element_line(color = "black", size = 1.3),
        panel.background = element_blank()) +
  scale_x_continuous(breaks = seq(-8, 8, 2)) + #    Customize the breaks in the x axis
  geom_text_repel(max.overlaps = Inf)


###############################################################################################################

#   Heatmap   #################################################################################################

library(ComplexHeatmap)


Best_Genes30 <- Result.df %>%
  arrange(padj) %>%
  drop_na(Symbol) %>%
  head(30)
Best_Genes30


mat <- counts(DEG, normalized = TRUE)[rownames(Best_Genes30), ]
head(mat)
mat.z <- t(apply(mat, 1, scale))    #   t = transpose
head(mat.z)
colnames(mat.z) <- c("SS_HS_2569",
                     "SS_HS_2570",
                     "SS_HS_2571",
                     "SS_HS_2573",
                     "SS_LS_2561",
                     "SS_LS_2562",
                     "SS_LS_2554",
                     "SS_LS_2555",
                     "SS_LS_2575")
Heatmap(mat.z,
        name = "HIPPOCAMPUS SS LS vs HS",
        cluster_rows = T,
        cluster_columns = T,
        column_labels = colnames(mat.z),
        row_labels = Best_Genes30$Symbol)



###############################################################################################################

#   Build PCA Plot    #########################################################################################


VSDF344vsAD <- vst(DEG, blind = FALSE)
plotPCA(VSDF344vsAD, intgroup = "condition") +
  geom_label(aes(label = colnames(mat.z)))



###############################################################################################################

#   DAVID Ontology    #########################################################################################

head(Best_FoldChange)
head(Best_Pvalue)

write.csv(Best_FoldChange, "HIPPOCAMPUS_HS vs LS/FOLDCHANGE/Best Genes with largest FoldChange_Hipocampus_SS_LS vs HS.csv")
write.csv(Best_Pvalue, "HIPPOCAMPUS_HS vs LS/PVALUE/Best Genes with lowest Pvalue_Hipocampus_SS_LS vs HS.csv")







#   FoldChange
url_foldchange <- read.table("https://davidbioinformatics.nih.gov/data/download/conv_04F65A412D081737572974006.txt",
                             header = TRUE,
                             sep = "\t",
                             quote = "\"")
colnames(url_foldchange) = c("Symbol", "ENSEMBL", "Species", "Description")
str(url_foldchange)
head(url_foldchange)



#   Convert Row Names of 'Best_FoldChange' into the First Column
ENSEMBL <- rownames(Best_FoldChange)
rownames(Best_FoldChange) <- NULL
ENSEMBL_col <- cbind(ENSEMBL, Best_FoldChange)



#   Merge the files by Gene Name
merged_data <- merge(ENSEMBL_col, url_foldchange, by = "ENSEMBL")
str(merged_data)


write.csv(merged_data, "HIPPOCAMPUS_HS vs LS/FOLDCHANGE/FINAL_FOLDCHANGE_Hippocampus_SS_LS vs HS.csv")








#   pvalue
url_pvalue <- read.table("https://davidbioinformatics.nih.gov/data/download/conv_04F65A412D081737573621385.txt",
                         header = TRUE,
                         sep = "\t",
                         quote = "\"")
colnames(url_pvalue) = c("Symbol", "ENSEMBL", "Species", "Description")
str(url_pvalue)



#   Convert Row Names of 'Best_Pvalue' into the First Column
ENSEMBL <- rownames(Best_Pvalue)
rownames(Best_Pvalue) <- NULL
ENSEMBL_col <- cbind(ENSEMBL, Best_Pvalue)



#   Merge the files by Gene Name
merged_data <- merge(ENSEMBL_col, url_pvalue, by = "ENSEMBL")
str(merged_data)


write.csv(merged_data, "HIPPOCAMPUS_HS vs LS/PVALUE/FINAL_PVALUE_Hippocampus_SS_LS vs HS.csv")



###############################################################################################################

#   Functional Analysis   #####################################################################################



#   BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(AnnotationDbi)

str(Result.df)
head(Result.df)








#   BP: Biological Processes (eg generation of neurons)   ############################################################

#   Upregulated Bar Plot
genes_to_test_up <- rownames(Best_Pvalue[Best_Pvalue$log2FoldChange > 0.5, ])
GO_results_up <- enrichGO(gene = genes_to_test_up, OrgDb = org.Rn.eg.db, keyType = "ENSEMBL", ont = "BP")
as.data.frame(GO_results_up)
head(GO_results_up)
barplot(GO_results_up,
        showCategory = 15) +
  ggtitle("Gene Ontology of up-regulated genes (Biological Processes)")



#   downregulated Bar Plot
genes_to_test_down <- rownames(Best_Pvalue[Best_Pvalue$log2FoldChange < -0.5, ])
GO_results_down <- enrichGO(gene = genes_to_test_down, OrgDb = org.Rn.eg.db, keyType = "ENSEMBL", ont = "BP")
head(GO_results_down)
as.data.frame(GO_results_down)

barplot(GO_results_down,
        showCategory = 15) +
  ggtitle("Gene Ontology of down-regulated genes (Biological Processes)")



#   MF: Molecular Function (eg actin binding)   ##################################################################

#   Upregulated Bar Plot
genes_to_test_up <- rownames(Best_Pvalue[Best_Pvalue$log2FoldChange > 0.5, ])
GO_results_up <- enrichGO(gene = genes_to_test_up, OrgDb = org.Rn.eg.db, keyType = "ENSEMBL", ont = "MF")
as.data.frame(GO_results_up)
head(GO_results_up)
barplot(GO_results_up,
        showCategory = 15) +
  ggtitle("Gene Ontology of up-regulated genes (Molecular Function)")



#   downregulated Bar Plot
genes_to_test_down <- rownames(Best_Pvalue[Best_Pvalue$log2FoldChange < -0.5, ])
GO_results_down <- enrichGO(gene = genes_to_test_down, OrgDb = org.Rn.eg.db, keyType = "ENSEMBL", ont = "MF")
head(GO_results_down)
as.data.frame(GO_results_down)

barplot(GO_results_down,
        showCategory = 15) +
  ggtitle("Gene Ontology of down-regulated genes (Molecular Function)")




#   CC: Cellular Component (eg cell jucntion)   ##################################################################   

#   Upregulated Bar Plot
genes_to_test_up <- rownames(Best_Pvalue[Best_Pvalue$log2FoldChange > 0.5, ])
GO_results_up <- enrichGO(gene = genes_to_test_up, OrgDb = org.Rn.eg.db, keyType = "ENSEMBL", ont = "CC")
as.data.frame(GO_results_up)
head(GO_results_up)
barplot(GO_results_up,
        showCategory = 15) +
  ggtitle("Gene Ontology of up-regulated genes (Cellular Component)")



#   downregulated Bar Plot
genes_to_test_down <- rownames(Best_Pvalue[Best_Pvalue$log2FoldChange < -0.5, ])
GO_results_down <- enrichGO(gene = genes_to_test_down, OrgDb = org.Rn.eg.db, keyType = "ENSEMBL", ont = "CC")
head(GO_results_down)
as.data.frame(GO_results_down)

barplot(GO_results_down,
        showCategory = 15) +
  ggtitle("Gene Ontology of down-regulated genes (Cellular Component)")
