# **Impact of ONECUT2 (OC2) on the Hypoxic Transcriptome in Prostate Cancer**

## **Project Overview**
This project investigates the transcriptional response of prostate cancer cell lines (LNCaP and PC3) to hypoxic conditions (low oxygen) and explores the regulatory role of the transcription factor ONECUT2 (OC2). Using a DESeq2-based workflow, we characterize how OC2 knockdown modulates the expression of hypoxia-induced genes, ultimately identifying a "Rescue Signature" where OC2 inhibition reverses hypoxic signaling.

## **Data Source & Biological Context**
The transcriptomic data used in this analysis is derived from a study published in Nature Communications (2019): https://pubmed.ncbi.nlm.nih.gov/30655535/

## **Key Biological Findings**
1) Hypoxia Response: Confirmed robust activation of HIF-1 signaling and Glycolysis pathways across both cell lines.
2) OC2 Regulation: Identified 729 DEGs in PC3 cells significantly regulated by OC2 knockdown under hypoxia.
3) Rescue Signature: Discovered a specific gene set where OC2 knockdown "rescues" (reverses) the effects of hypoxia, suggesting OC2 is a critical mediator of metabolic adaptation in prostate cancer.

## **Workflow & Methodology**
1. Quality Control & Library Composition
Before differential expression, we evaluated the raw data for technical artifacts:

Library Depth: Assessed sequencing depth across replicates to ensure statistical power.

Composition Analysis: Quantified Mitochondrial and Ribosomal content. High mitochondrial percentages can indicate low-quality samples, but our data showed consistent nuclear-dominated libraries.

2. Differential Expression Analysis (DESeq2)
We utilized a Generalized Linear Model (GLM) to estimate fold changes:

Design: ~ Treatment + Genetic_perturbation + Treatment:Genetic_perturbation.

Shrinkage: Applied apeglm shrinkage to handle low-count genes and provide conservative effect size estimates.

Visualization: Used MA plots and Volcano plots to validate model behavior and identify significant genes (padj < 0.05, |log2FC| > 1).

3. Pathway Enrichment
To move from gene lists to biological mechanisms, we performed enrichment analysis:

KEGG/GO: Validated that upregulated genes were enriched for anaerobic metabolism and HIF-1 target genes.

4. OC2 Rescue Signature
The highlight of this analysis is the identification of genes where:

Hypoxia causes Up/Down-regulation.

OC2 Knockdown causes the opposite effect (Down/Up-regulation).
This was visualized using ComplexHeatmap to show the reversal of Z-scores across conditions.

## **Repository Structure**
/data: Contains analysis-ready count matrices and metadata.

/scripts: Custom R functions for visualization (functions.R).

/
