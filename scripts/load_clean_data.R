
# Load required libraries for data manipulation and genomic annotation
library(dplyr)
library(tidyverse)
library(org.Hs.eg.db)
library(here)
library(AnnotationDbi)
library(conflicted)

# Load raw count data from TSV file
# Using as.data.frame to ensure compatibility with AnnotationDbi and base R merging
raw_data <- as.data.frame(read_tsv(here("data","raw_counts.tsv")))

# --- DATA FILTERING AND SAMPLE RENAMING ---

# 1. Create the mapping variable (First 2 columns only)
geo_metadata <- data.frame(
  GSM_Accession = c("GeneID",
    "GSM3145509", "GSM3145510", "GSM3145511", "GSM3145512", # LNCaP Normoxia
    "GSM3145513", "GSM3145514", "GSM3145515", "GSM3145516", # LNCaP Hypoxia
    "GSM3145517", "GSM3145518", "GSM3145519", "GSM3145520", # PC3 Normoxia
    "GSM3145521", "GSM3145522", "GSM3145523", "GSM3145524"  # PC3 Hypoxia
  ),
  Original_Title = c("GeneID",
    "LNCaP_RNA-Seq_Empty_Vector_Normoxia_rep1", "LNCaP_RNA-Seq_Empty_Vector_Normoxia_rep2",
    "LNCaP_RNA-Seq_Flag-OC2_Overexpression_Normoxia_rep1", "LNCaP_RNA-Seq_Flag-OC2_Overexpression_Normoxia_rep2",
    "LNCaP_RNA-Seq_Empty_Vector_Hypoxia_rep1", "LNCaP_RNA-Seq_Empty_Vector_Hypoxia_rep2",
    "LNCaP_RNA-Seq_Flag-OC2_Overexpression_Hypoxia_rep1", "LNCaP_RNA-Seq_Flag-OC2_Overexpression_Hypoxia_rep2",
    "PC3_RNA-Seq_siCtrl_Normoxia_rep1", "PC3_RNA-Seq_siCtrl_Normoxia_rep2",
    "PC3_RNA-Seq_siOC2_Normoxia_siRNA1", "PC3_RNA-Seq_siOC2_Normoxia_siRNA2",
    "PC3_RNA-Seq_siCtrl_Hypoxia_rep1", "PC3_RNA-Seq_siCtrl_Hypoxia_rep2",
    "PC3_RNA-Seq_siOC2_Hypoxia_siRNA1", "PC3_RNA-Seq_siOC2_Hypoxia_siRNA2"
  )
)

# 2. Clean titles: replace Empty_Vector -> EV and remove RNA-Seq
geo_metadata <- geo_metadata %>%
  mutate(Clean_Title = Original_Title %>%
           str_replace_all("Empty_Vector", "EV") %>%
           str_replace_all("_RNA-Seq", "")%>%
          str_replace_all("Flag-OC2_Overexpression","OC2_overexpress")) 

# 3. Assign the cleaned titles to your data matrix
# (Ensure your matrix columns are in the same order as these GSMs)
stopifnot(all(colnames(raw_data)==geo_metadata$GSM_Accession))
colnames(raw_data) <- geo_metadata$Clean_Title

head(raw_data)

# --- GENE ANNOTATION ---

# Ensure Entrez IDs are formatted as characters for database compatibility
entrez_ids <- as.character(raw_data$GeneID)

# Query org.Hs.eg.db for Symbol, Full Gene Name, and Gene Type
# Keytype is set to ENTREZID to match the input data
annotations <- AnnotationDbi::select(org.Hs.eg.db, 
                                     keys = entrez_ids, 
                                     columns = c("SYMBOL", "GENENAME", "GENETYPE"), 
                                     keytype = "ENTREZID")

# Merge annotations with expression data
# A left join (all.x = TRUE) is performed to preserve all original rows
annotated_data <- merge(raw_data, annotations, by.x = "GeneID", by.y = "ENTREZID", all.x = TRUE)

# --- FILE EXPORT ---

# Export the final annotated data frame to a CSV file for downstream analysis
# row.names is set to FALSE to prevent writing the numeric index column
write.csv(annotated_data, 
          file = here("data","annotated_raw_data.csv"), 
          row.names = FALSE)


# Cleaning_for_analysis_file ----------------------------------------------
# Force R to use dplyr for select and filter
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

#remove rows with no gene symbols
annotated_data <- annotated_data %>%
  filter(!is.na(SYMBOL))

# 1. Handle Duplicated Symbols
# make.unique() appends .1, .2, etc., to duplicate strings
# This ensures every row name is distinct while preserving the symbol
original_symbols <- annotated_data$SYMBOL
unique_symbols <- make.unique(ifelse(is.na(original_symbols), "Unknown", original_symbols), sep = "_dup")

# 2. Create the metadata_gene dataframe
# This preserves the link between the original IDs and the new unique symbols
metadata_gene <- annotated_data %>%
  select(GeneID, SYMBOL, GENENAME, GENETYPE) %>%
  mutate(Unique_Symbol = unique_symbols)

# 3. Prepare the count data matrix
# Remove the annotation columns and the original GeneID
count_data <- annotated_data %>%
  select(-GeneID, -SYMBOL, -GENENAME, -GENETYPE) %>%
  as.data.frame()

# Assign the unique symbols as row names
rownames(count_data) <- unique_symbols

#make sure no NA values in the df
stopifnot(!any(is.na(count_data)))

#check rownames are same as Unique_symbol column
stopifnot(identical(rownames(count_data), metadata_gene$Unique_Symbol))

# filtering based on low expression ---------------------------------------
keep_non_zero <- rowSums(count_data == 0) <= 7

# Apply the filter to both the expression matrix and the metadata
count_data_fil<- count_data[keep_non_zero, ]
metadata_gene_fil <- metadata_gene[keep_non_zero, ]
# Verify the results
message("Rows remaining: ", nrow(count_data_fil))
message("Rows removed (all zeros): ", sum(!keep_non_zero))

# Final alignment check
stopifnot(identical(rownames(count_data_fil), metadata_gene_fil$Unique_Symbol))

write.csv(count_data_fil, 
          file = here("data","analysis_ready_count_data.csv"), 
          row.names = T)

write.csv(metadata_gene_fil, 
          file = here("data","gene_meatadata_analysis_ready_count.csv"), 
          row.names = F)
