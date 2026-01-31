library(DESeq2)
library(ggplot2)
library(pheatmap)
library(patchwork)
library(ggplotify)
library(RColorBrewer)

QC_pca_heatmap <- function(dds,metadata) {
  
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
    # Use paste0 to combine the text with your calculated variance
    xlab(paste0("PC1(", percentVar[1], "% variance)")) +
    ylab(paste0("PC2(", percentVar[2], "% variance)")) +
    labs(color = "Genetic perturbation",
         shape = "Treatment",
         title = "PCA between samples") +
    coord_fixed() + # Important: ensures the axes are proportional to the variance
    theme_classic() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  # 5. Correlation Heatmap Preparation
  rld_mat <- assay(rld)
  rld_cor <- cor(rld_mat)
  
  # Clean names for visualization
  rownames(rld_cor) <- gsub("PC3_|LNCaP_", "", rownames(rld_cor))
  colnames(rld_cor) <- gsub("PC3_|LNCaP_", "", colnames(rld_cor))
  
  # Clean annotation names to match
  annotation_col <- metadata
  rownames(annotation_col) <- gsub("PC3_|LNCaP_", "", rownames(annotation_col))
  
  # Capture pheatmap as a ggplot-compatible object
  heat_p<-as.ggplot(pheatmap(
    rld_cor,
    annotation_col = annotation_col,
    color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
    main = "Correlation between samples",
    treeheight_col = 10,
    treeheight_row = 10,
    border_color = "grey",
    fontsize = 8,
    fontsize_row = 8,
    fontsize_col = 8,
    legend = TRUE,
    legend_width = 0.001,silent = TRUE
  ))
  
  return(list(pca = pca, heatmap = heat_p))
}