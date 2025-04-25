setwd("C:/Users/danielheidar/OneDrive - Menntaský/RNAseq")
ALKBH3_TPM_data <- read.table("ALKBH3_DESeq2_TPM_values.tsv", header = TRUE, stringsAsFactors = FALSE, row.names = NULL)

ALKBH3_TPM_data$siALKBH3_1_DAM_ratio <- ALKBH3_TPM_data$siALKBH3_nuc1_DAM / ALKBH3_TPM_data$siALKBH3_cyto1_DAM

ALKBH3_TPM_data$siALKBH3_1_KAK_ratio <- ALKBH3_TPM_data$siALKBH3_nuc1_KAK / ALKBH3_TPM_data$siALKBH3_cyto1_KAK

ALKBH3_TPM_data$siALKBH3_2_DAM_ratio <- ALKBH3_TPM_data$siALKBH3_nuc2_DAM / ALKBH3_TPM_data$siALKBH3_cyto2_DAM

ALKBH3_TPM_data$siALKBH3_2_KAK_ratio <- ALKBH3_TPM_data$siALKBH3_nuc2_KAK / ALKBH3_TPM_data$siALKBH3_cyto2_KAK

ALKBH3_TPM_data$siALKBH3_3_DAM_ratio <- ALKBH3_TPM_data$siALKBH3_nuc3_DAM / ALKBH3_TPM_data$siALKBH3_cyto3_DAM


ALKBH3_TPM_data$mean_ratio <- rowMeans(ALKBH3_TPM_data[, c("siALKBH3_1_DAM_ratio", "siALKBH3_1_KAK_ratio", "siALKBH3_2_DAM_ratio", "siALKBH3_2_KAK_ratio", "siALKBH3_3_DAM_ratio")])


search_value <- "ENSG00000111640"
result <- ALKBH3_TPM_data[ALKBH3_TPM_data$row.names == search_value, ]

# Print the row
print(result)
GADPH_mean <- result$mean_ratio




ALKBH3_TPM_data$GADPH_ratio <- ALKBH3_TPM_data$mean_ratio / GADPH_mean


# Print the row
print(result)
GADPH_mean <- result$mean_ratio

numeric_cols <- sapply(ALKBH3_TPM_data, is.numeric)

last_col <- ncol(ALKBH3_TPM_data)


ALKBH3_TPM_data_clean <- ALKBH3_TPM_data[!is.na(ALKBH3_TPM_data[, last_col]) & !is.infinite(ALKBH3_TPM_data[, last_col]), ]


write_csv(FTO_TPM_data_clean, "FTO_TPM_data_ratio.csv")


--------------------------------------------------------
  
ALKBH3_TPM_data_clean$siCTRL_1_DAM_ratio <- ALKBH3_TPM_data_clean$siCTRL_nuc1_DAM / ALKBH3_TPM_data_clean$siCTRL_cyto1_DAM

ALKBH3_TPM_data_clean$siCTRL_1_KAK_ratio <- ALKBH3_TPM_data_clean$siCTRL_nuc1_KAK / ALKBH3_TPM_data_clean$siCTRL_cyto1_KAK

ALKBH3_TPM_data_clean$siCTRL_2_DAM_ratio <- ALKBH3_TPM_data_clean$siCTRL_nuc2_DAM / ALKBH3_TPM_data_clean$siCTRL_cyto2_DAM

ALKBH3_TPM_data_clean$siCTRL_2_KAK_ratio <- ALKBH3_TPM_data_clean$siCTRL_nuc2_KAK / ALKBH3_TPM_data_clean$siCTRL_cyto2_KAK

ALKBH3_TPM_data_clean$siCTRL_3_DAM_ratio <- ALKBH3_TPM_data_clean$siCTRL_nuc3_DAM / ALKBH3_TPM_data_clean$siCTRL_cyto3_DAM


ALKBH3_TPM_data_clean$mean_ratio_ctrl <- rowMeans(ALKBH3_TPM_data_clean[, c("siCTRL_1_DAM_ratio", "siCTRL_1_KAK_ratio", "siCTRL_2_DAM_ratio", "siCTRL_2_KAK_ratio", "siCTRL_3_DAM_ratio")])


search_value <- "ENSG00000111640"
result <- ALKBH3_TPM_data_clean[ALKBH3_TPM_data_clean$row.names == search_value, ]

# Print the row
print(result)
GADPH_mean <- result$mean_ratio_ctrl




ALKBH3_TPM_data_clean$GADPH_ratio_ctrl <- ALKBH3_TPM_data_clean$mean_ratio_ctrl / GADPH_mean

ALKBH3_TPM_data_clean$Final_ratio <- ALKBH3_TPM_data_clean$GADPH_ratio / ALKBH3_TPM_data_clean$GADPH_ratio_ctrl

write_csv(FTO_TPM_data_clean, "FTO_TPM_data_ratio.csv")



last_col <- ncol(ALKBH3_TPM_data_clean)

ALKBH3_TPM_data_clean <- ALKBH3_TPM_data_clean[!is.na(ALKBH3_TPM_data_clean[, last_col]) & !is.infinite(ALKBH3_TPM_data_clean[, last_col]), ]


mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_names_RNAseq<- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
                          filters = "ensembl_gene_id",
                          values = ALKBH3_TPM_data_clean$row.names,
                          mart = mart)

merged_data <- merge(
  data.frame(transcript_id = ALKBH3_TPM_data_clean$row.names, ALKBH3_TPM_data_clean),
  gene_names_RNAseq,
  by.x = "transcript_id",
  by.y = "ensembl_gene_id",
  all.x = TRUE
)

merged_data$external_gene_name[is.na(merged_data$external_gene_name)] <- "Unknown"

write_csv(merged_data, "ALKBH3_data_ratio_test.csv")


merged_data <- read.csv2("ALKBH3_data_ratio.csv")
-----
  
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
    title = "Ratio Intervals in siALKBH3 Genes",
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
