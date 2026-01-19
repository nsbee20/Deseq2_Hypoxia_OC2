
# Load required libraries for data manipulation and genomic annotation
library(dplyr)
library(tidyverse)
library(org.Hs.eg.db)
library(here)
library(AnnotationDbi)

# Load raw count data from TSV file
# Using as.data.frame to ensure compatibility with AnnotationDbi and base R merging
raw_data <- as.data.frame(read_tsv("raw_counts.tsv"))

# --- DATA FILTERING AND SAMPLE RENAMING ---

# Define mapping between GSM identifiers and descriptive sample names
# Samples include LNCaP and PC3 cell lines under Normoxia and Hypoxia conditions
gsm_map <- c(
  "GSM3145509" = "LNCaP_RNA-Seq_Empty_Vector_Normoxia_rep1",
  "GSM3145510" = "LNCaP_RNA-Seq_Empty_Vector_Normoxia_rep2",
  "GSM3145513" = "LNCaP_RNA-Seq_Empty_Vector_Hypoxia_rep1",
  "GSM3145514" = "LNCaP_RNA-Seq_Empty_Vector_Hypoxia_rep2",
  "GSM3145517" = "PC3_RNA-Seq_siCtrl_Normoxia_rep1",
  "GSM3145518" = "PC3_RNA-Seq_siCtrl_Normoxia_rep2",
  "GSM3145521" = "PC3_RNA-Seq_siCtrl_Hypoxia_rep1",
  "GSM3145522" = "PC3_RNA-Seq_siCtrl_Hypoxia_rep2"
)

# Identify target columns: the GeneID primary key and relevant GSM samples
gene_id_col <- colnames(raw_data)[1] 
gsm_cols <- intersect(colnames(raw_data), names(gsm_map))

# Subset dataset to retain only the identified columns
data_filtered <- raw_data[, c(gene_id_col, gsm_cols)]

# Rename GSM columns to descriptive names using the mapping vector
colnames(data_filtered) <- ifelse(colnames(data_filtered) %in% names(gsm_map), 
                                  gsm_map[colnames(data_filtered)], 
                                  colnames(data_filtered))

# --- GENE ANNOTATION ---

# Ensure Entrez IDs are formatted as characters for database compatibility
entrez_ids <- as.character(data_filtered$GeneID)

# Query org.Hs.eg.db for Symbol, Full Gene Name, and Gene Type
# Keytype is set to ENTREZID to match the input data
annotations <- AnnotationDbi::select(org.Hs.eg.db, 
                                     keys = entrez_ids, 
                                     columns = c("SYMBOL", "GENENAME", "GENETYPE"), 
                                     keytype = "ENTREZID")

# Merge annotations with expression data
# A left join (all.x = TRUE) is performed to preserve all original rows
data_annotated <- merge(data_filtered, annotations, by.x = "GeneID", by.y = "ENTREZID", all.x = TRUE)

# --- FILE EXPORT ---

# Export the final annotated data frame to a CSV file for downstream analysis
# row.names is set to FALSE to prevent writing the numeric index column
write.csv(data_annotated, 
          file = "annotated_raw_data.csv", 
          row.names = FALSE)
