library(DESeq2)
library(ggplot2)
library(pheatmap)
library(patchwork)
library(ggplotify)
library(RColorBrewer)
library(clusterProfiler)
library(EnhancedVolcano)

QC_pca_heatmap <- function(dds,metadata,cell_line_name) {
  
  # 2. Transform data (rlog)
  rld <- rlog(dds, blind = TRUE)
  
  # 3. PCA Calculation
  pcaData <- plotPCA(rld, intgroup = c("Genetic_perturbation", "Treatment"), returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"), 2)
  
  # 4. Generate PCA Plot
  pca<-ggplot(pcaData, aes(PC1, PC2,
                           color = Genetic_perturbation,
                           shape = Treatment)) +
    geom_point(size = 4, alpha = 0.9) +
    xlab(paste0("PC1(", percentVar[1], "% variance)")) +
    ylab(paste0("PC2(", percentVar[2], "% variance)")) +
    labs(color = "Genetic perturbation",
         shape = "Treatment",title = paste("PCA:", cell_line_name)) +
    coord_fixed() + 
    theme_classic() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  # 5. Correlation of samples
  rld_mat <- assay(rld)
  rld_cor <- cor(rld_mat)
  
  # Clean names for visualization
  rownames(rld_cor) <- gsub("PC3_|LNCaP_", "", rownames(rld_cor))
  colnames(rld_cor) <- gsub("PC3_|LNCaP_", "", colnames(rld_cor))
  
  # Clean annotation names to match
  annotation_col <- metadata
  rownames(annotation_col) <- gsub("PC3_|LNCaP_", "", rownames(annotation_col))
  
  # Capture pheatmap as a ggplot-compatible object
  heat_p<-(as.ggplot(pheatmap(
    rld_cor,
    annotation_col = annotation_col,
    color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
    main = paste("Correlation Matrix:", cell_line_name),
    treeheight_col = 10,
    treeheight_row = 10,
    border_color = "grey",
    fontsize = 8,
    fontsize_row = 8,
    fontsize_col = 8,
    legend = TRUE,
    legend_width = 0.001,silent = TRUE
  ))+ theme_void() + 
    theme(
      # Ensure the title still shows up and is positioned correctly
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 10)),
      # Force white background for GitHub rendering
      plot.background = element_rect(fill = "white", color = NA)
    ))
  
  return(list(pca = pca, heatmap = heat_p))
}


# Function to generate consistent volcano plots
plot_volcano <- function(res, title_str,...) {
  # Convert to data frame to handle logical indexing
  df <- as.data.frame(res)
  top_labels <- head(rownames(res[order(res$padj), ]), 30)
  
  # Create a custom color vector
  keyvals <- ifelse(
    df$padj < 0.05 & df$log2FoldChange > 1.0, 'red3',
    ifelse(df$padj < 0.05 & df$log2FoldChange < -1.0, 'royalblue',
           'grey'))
  
  #Convert any remaining NAs to 'grey'
  keyvals[is.na(keyvals)] <- 'grey'
  
  # Ensure the vector has names for the legend
  names(keyvals)[keyvals == 'red3'] <- 'Upregulated'
  names(keyvals)[keyvals == 'grey'] <- 'NS / Low Fold'
  names(keyvals)[keyvals == 'royalblue'] <- 'Downregulated'
  
  EnhancedVolcano(res,
                  lab = rownames(res),
                  selectLab = top_labels,
                  x = 'log2FoldChange',
                  y = 'padj', 
                  title = title_str,
                  subtitle = NULL,      # Removes the "EnhancedVolcano" subtitle
                  caption = NULL,       # Removes the "total variables" footer
                  ylab = bquote(~-Log[10] ~ italic(padj)), 
                  pCutoff = 0.05,
                  FCcutoff = 1.0, 
                  colCustom = keyvals, # Use the custom color vector
                  labSize = 5,
                  colAlpha = 0.6,
                  drawConnectors = TRUE,
                  widthConnectors = 0.5,
                  max.overlaps = 20,          # Limits label clutter
                  lengthConnectors = unit(0.01, "npc"), # Length of the lines
                  ...
  )
}

#function for plotting kegg pathways enriched in the DEGs
plot_kegg_pathways<- function(gene_list, title_str, p_cutoff = 0.05, max_pathways = 20) {
  
  # 1. Run the KEGG enrichment inside the function
  kegg_res <- clusterProfiler::enrichKEGG(
    gene = gene_list, 
    organism = 'hsa', 
    pvalueCutoff = p_cutoff
  )
  
  # 2. Convert to data frame
  df <- as.data.frame(kegg_res)
  
  # 3. Check if any pathways were found
  if (nrow(df) == 0) {
    return(message("No enriched pathways found for ", title_str, " at p < ", p_cutoff))
  }
  
  # 4. Filter to top N and calculate -log10(padj)
  df <- df[1:min(max_pathways, nrow(df)), ] 
  df$neg_log_padj <- -log10(df$p.adjust)
  
  # 5. Reorder factor levels for the y-axis (most significant at top)
  df$Description <- factor(df$Description, 
                           levels = df$Description[order(df$neg_log_padj)])
  
  # 6. Build the Plot
  ggplot2::ggplot(df, aes(x = neg_log_padj, y = Description)) +
    geom_segment(aes(x = 0, xend = neg_log_padj, y = Description, yend = Description), 
                 color = "lightgrey") +
    # Color is outside aes() to avoid the legend
    geom_point(aes(size = Count), color = "#f8766d") + 
    theme_minimal(base_size = 16) +
    labs(
      title = title_str,
      x = expression(-log[10](italic(padj))), 
      y = "Pathway",
      size = "Gene Count"
    ) +
    theme(
      plot.title = element_text(hjust = 0, face = "bold", size = 12),
      #panel.grid.minor = element_blank(),
      #panel.grid.major.y = element_blank(),
      legend.position = "right"
    ) +
    # Dashed line at the significance threshold
    geom_vline(xintercept = -log10(p_cutoff), linetype = "dashed", color = "black")
}