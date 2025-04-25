setwd("C:/Users/danielheidar/OneDrive - Menntaský/Nanopore")

significant_P <- read.csv2('Significant_P.csv', sep=";")

install.packages("AnnotationDbi")
BiocManager::install("org.Hs.eg.db")

# Load packages
library(AnnotationDbi)
library(org.Hs.eg.db)

gene_symbols <- select(org.Hs.eg.db, keys = significant_P$cleaned, keytype = "REFSEQ", columns = "SYMBOL")

# Check the results
head(gene_symbols)

significant_P$cleaned <- sub("\\..*", "", significant_P$ref_id)

-------------------
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

# Get gene symbols using RefSeq IDs
gene_symbols_biomart <- getBM(attributes = c("refseq_mrna", "hgnc_symbol"),
                              filters = "refseq_mrna",
                              values = significant_P$cleaned,
                              mart = ensembl)

# View the results
head(gene_symbols_biomart)

merged_data_nano <- merge(
  data.frame(transcript_id = significant_P$cleaned, significant_P),
  gene_symbols_biomart,
  by.x = "transcript_id",
  by.y = "refseq_mrna",
  all.x = TRUE
)

write_csv(merged_data_nano, "Significant_P_genes.csv")

common_ALK <- intersect(filtered_data_ALKBH3$external_gene_name, merged_data_nano$hgnc_symbol)

common_FTO <- intersect(filtered_data_FTO$external_gene_name, merged_data_nano$hgnc_symbol)