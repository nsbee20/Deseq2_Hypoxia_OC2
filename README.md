# **Impact of ONECUT2 (OC2) on the Hypoxic Transcriptome in Prostate Cancer**

## **Project Overview**
This project investigates the transcriptional response of prostate cancer cell lines (LNCaP and PC3) to hypoxic conditions (low oxygen) and explores the regulatory role of the transcription factor ONECUT2 (OC2). Using a DESeq2-based workflow, we characterize how OC2 knockdown modulates the expression of hypoxia-induced genes, ultimately identifying a "Rescue Signature" where OC2 inhibition reverses hypoxic signaling.

## **Data Source & Biological Context**
The transcriptomic data used in this analysis is derived from a study published in Nature Communications (2019): https://pubmed.ncbi.nlm.nih.gov/30655535/

## **Key Biological Question and Findings**
1) Which genes and pathways are regulated in prostrate cancer cell lines under hypoxic conditions?

   Identified genes regulated under hypoxic conditions in both cell lines using Deseq2 package. 
     Confirmed robust activation of HIF-1 signaling and Glycolysis pathways across both cell lines.

2) How does OC2 regulate transcriptomic profile under hypoxic conditions?
   
   Identified 100s of DEGs in PC3 cells under hypoxia and normal conditions when OC2 is knocked down. 
   
3) What are the genes that OC2 regulates under hypoxia conditions? Can we identify a "rescue gene signature":
   
   Discovered a specific gene set where OC2 knockdown "rescues" (reverses) the transcriptomic changes driven by hypoxia. This analysis can be utlized to further understand the mechanism of action of OC2 in tumor progression under hypoxia conditions.  

## **Workflow & Methodology**
1. Quality Control & Library Composition
Before differential expression is performed, raw data was evaluated for technical artifacts:
  a) Library Depth: Assessed sequencing depth across replicates to ensure statistical power.
  b) Composition Analysis: Quantified Mitochondrial and Ribosomal content and nuclear/other RNA content. 
  c) Most of this analysis is available in Initial_QC qmd and md files. 

2. Differential Expression Analysis using DESeq2: 
  a) Utilized a Generalized Linear Model (GLM) to estimate fold changes to understand the affect of treatment and genetic perturbation for both cell types.
  b) Visualization, QC, stepwise procedure for Deseq2: Displayed step wise procedure for identifying DEGs using Deseq2 package. Confirmed and validated findings at each step using various       QC steps and visualization tools. 
  
4. Pathway Enrichment: To move from gene lists to biological mechanisms, we performed enrichment analysis: Validated that regulated genes were enriched for Glycolysis and HIF-1 target genes (well known mechanism under hypoxia condition).

5. OC2 Rescue Signature under hypoxix conditions: Using this analysis we identified 54 genes that are regulated in the opposite direction by OC2 Knockdown under hypoxic conditions. 

Step 2-4 analysis are available in DGEs_Deseq2 qmd and md files. 

## **Repository Structure**
/data: Contains analysis-ready count matrices and metadata.

/scripts: Custom R functions for visualization (functions.R).

/analysis: Contains md and qmd files of analysis. 
