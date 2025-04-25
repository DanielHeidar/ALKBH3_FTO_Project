setwd("C:/Users/danielheidar/OneDrive - Menntaský/RNAseq")
FTO_TPM_data <- read.table("FTO_DESeq2_TPM_values.tsv", header = TRUE, stringsAsFactors = FALSE, row.names = NULL)

FTO_TPM_data$siFTO_1_DAM_ratio <- FTO_TPM_data$siFTO_nuc1_DAM / FTO_TPM_data$siFTO_cyto1_DAM

FTO_TPM_data$siFTO_1_KAK_ratio <- FTO_TPM_data$siFTO_nuc1_KAK / FTO_TPM_data$siFTO_cyto1_KAK

FTO_TPM_data$siFTO_2_DAM_ratio <- FTO_TPM_data$siFTO_nuc2_DAM / FTO_TPM_data$siFTO_cyto2_DAM

FTO_TPM_data$siFTO_2_KAK_ratio <- FTO_TPM_data$siFTO_nuc2_KAK / FTO_TPM_data$siFTO_cyto2_KAK


FTO_TPM_data$mean_ratio <- rowMeans(FTO_TPM_data[, c("siFTO_1_DAM_ratio", "siFTO_1_KAK_ratio", "siFTO_2_DAM_ratio", "siFTO_2_KAK_ratio")])


search_value <- "ENSG00000111640"
result <- FTO_TPM_data[FTO_TPM_data$row.names == search_value, ]

# Print the row
print(result)
GADPH_mean <- result$mean_ratio




FTO_TPM_data$GADPH_ratio <- FTO_TPM_data$mean_ratio / GADPH_mean
  
search_value <- "ENSG00000111640"
result <- FTO_TPM_data[FTO_TPM_data$row.names == search_value, ]

# Print the row
print(result)
GADPH_mean <- result$mean_ratio

numeric_cols <- sapply(FTO_TPM_data, is.numeric)

last_col <- ncol(FTO_TPM_data)


FTO_TPM_data_clean <- FTO_TPM_data[!is.na(FTO_TPM_data[, last_col]) & !is.infinite(FTO_TPM_data[, last_col]), ]


write_csv(FTO_TPM_data_clean, "FTO_TPM_data_ratio.csv")


--------------------------------------------------------
  

FTO_TPM_data_clean$siCTRL_1_DAM_ratio <- FTO_TPM_data_clean$siCTRL_nuc1_DAM / FTO_TPM_data_clean$siCTRL_cyto1_DAM

FTO_TPM_data_clean$siCTRL_1_KAK_ratio <- FTO_TPM_data_clean$siCTRL_nuc1_KAK / FTO_TPM_data_clean$siCTRL_cyto1_KAK

FTO_TPM_data_clean$siCTRL_2_DAM_ratio <- FTO_TPM_data_clean$siCTRL_nuc2_DAM / FTO_TPM_data_clean$siCTRL_cyto2_DAM

FTO_TPM_data_clean$siCTRL_2_KAK_ratio <- FTO_TPM_data_clean$siCTRL_nuc2_KAK / FTO_TPM_data_clean$siCTRL_cyto2_KAK


FTO_TPM_data_clean$mean_ratio_ctrl <- rowMeans(FTO_TPM_data_clean[, c("siCTRL_1_DAM_ratio", "siCTRL_1_KAK_ratio", "siCTRL_2_DAM_ratio", "siCTRL_2_KAK_ratio")])


search_value <- "ENSG00000111640"
result <- FTO_TPM_data_clean[FTO_TPM_data_clean$row.names == search_value, ]

# Print the row
print(result)
GADPH_mean <- result$mean_ratio_ctrl




FTO_TPM_data_clean$GADPH_ratio_ctrl <- FTO_TPM_data_clean$mean_ratio_ctrl / GADPH_mean

FTO_TPM_data_clean$Final_ratio <- FTO_TPM_data_clean$GADPH_ratio / FTO_TPM_data_clean$GADPH_ratio_ctrl

write_csv(FTO_TPM_data_clean, "FTO_TPM_data_ratio.csv")



last_col <- ncol(FTO_TPM_data_clean)

FTO_TPM_data_clean <- FTO_TPM_data_clean[!is.na(FTO_TPM_data_clean[, last_col]) & !is.infinite(FTO_TPM_data_clean[, last_col]), ]


mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_names_RNAseq<- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
                       filters = "ensembl_gene_id",
                       values = FTO_TPM_data_clean$row.names,
                       mart = mart)
merged_data <- merge(
  data.frame(transcript_id = FTO_TPM_data_clean$row.names, FTO_TPM_data_clean),
  gene_names_RNAseq,
  by.x = "transcript_id",
  by.y = "ensembl_gene_id",
  all.x = TRUE
)

merged_data$external_gene_name[is.na(merged_data$external_gene_name)] <- "Unknown"

data 

write_csv(merged_data, "FTO_data_ratio_2.csv")

merged_data <- read.csv2("FTO_data_ratio_2.csv")

----------------------------------------------------------------------
  
merged_data_FTO <- merged_data

filtered_data_FTO <- merged_data_FTO[merged_data_FTO$Final_ratio > 1.5 & merged_data_FTO$Final_ratio < 11.2, ]

entrez_ids <- bitr(filtered_data_FTO$external_gene_name, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
head(entrez_ids)

library(org.Hs.eg.db)

mapped_genes <- bitr(entrez_ids$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
length(mapped_genes$ENTREZID)  # Ensure you have enough mapped genes

go_results <- enrichGO(
  gene          = entrez_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "ALL",    # Options: "BP" (Biological Process), "MF" (Molecular Function), "CC" (Cellular Component)
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

head(go_results)

test_genes <- c("TP53", "BRCA1", "EGFR", "MYC", "AKT1")
test_results <- enrichGO(
  gene = bitr(test_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "ALL",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)
head(go_results)
------------------------------------------------------

merged_data_FTO$f

ggplot(filtered_data_FTO, aes(x = Final_ratio)) +
  geom_histogram(binwidth = 0.5, fill = "skyblue", color = "black", alpha = 0.7) +
  labs(
    title = "Distribution of Numbers in Column",
    x = "Values",
    y = "Frequency"
  ) +
  theme_minimal()

ggplot(filtered_data_ALKBH3, aes(x = Final_ratio)) +
  geom_density(fill = "skyblue", alpha = 0.7) +
  labs(
    title = "Density Plot of Values",
    x = "Values",
    y = "Density"
  ) +
  theme_minimal()


ggplot(filtered_data_FTO, aes(x = Final_ratio)) +
  geom_histogram(aes(y = ..density..), binwidth = 0.5, fill = "skyblue", color = "black", alpha = 0.5) +
  geom_density(color = "darkblue", size = 1) +
  labs(
    title = "Histogram with Density Overlay",
    x = "Values",
    y = "Density"
  ) +
  theme_minimal()
----------------------------------------------------------------


gene_mapping_unique <- gene_names_RNAseq[!duplicated(gene_names_RNAseq$ensembl_transcript_id), ]

install.packages("tidyverse")
library(tidyverse)

# Load necessary libraries
library(dplyr)
library(ggplot2)

# Filter and process the data
dataf <- data.frame(merged_data$Difference.ratio.Experimental.control) %>%
  filter(
    !is.na(merged_data.Difference.ratio.Experimental.control),
    !is.nan(merged_data.Difference.ratio.Experimental.control),
    !is.infinite(merged_data.Difference.ratio.Experimental.control)
  ) %>%
  mutate(merged_data.Difference.ratio.Experimental.control = as.numeric(as.character(merged_data.Difference.ratio.Experimental.control))) %>%
  mutate(interval = case_when(
    merged_data.Difference.ratio.Experimental.control < 1.5 ~ "0 - 1.5",
    merged_data.Difference.ratio.Experimental.control >= 1.5 & merged_data.Difference.ratio.Experimental.control <= 2 ~ "1.5 - 2",
    merged_data.Difference.ratio.Experimental.control > 2 & merged_data.Difference.ratio.Experimental.control <= 5 ~ "2 - 5",
    merged_data.Difference.ratio.Experimental.control > 5 ~ "5+"
  ))

# Group data by interval and calculate counts
interval_counts <- dataf %>%
  group_by(interval) %>%
  summarise(count = n(), .groups = 'drop') %>%
  arrange(interval)  # Ensure intervals are ordered logically

# Create the barplot
ggplot(interval_counts, aes(x = interval, y = count, fill = interval)) +
  geom_bar(stat = "identity", color = "black", size = 0.3) +  # Add borders to bars
  geom_text(
    aes(label = count),
    vjust = -0.5,
    size = 4.5,
    fontface = "bold"
  ) +  # Add count labels above bars
  scale_fill_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3")) +  # Custom color palette
  labs(
    title = "Ratio Intervals in siFTO Genes",
    x = "Ratio Interval",
    y = "Number of Transcripts"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(face = "bold"),
    legend.position = "none",
    panel.grid.major.y = element_line(color = "gray80", linetype = "dashed")
  ) +
  ylim(0, max(interval_counts$count) + 5)  # Adjust y-axis limits for labels



str(dataf$FTO_TPM_data_clean.Final_ratio)
summary(dataf$FTO_TPM_data_clean.Final_ratio)
