# NASH Progression Atlas

A multi-cohort liver transcriptomics meta-analysis framework to identify
robust disease progression signatures, integrate serum miRNA regulation,
and characterize cell-type dynamics across NAFLD → NASH → fibrosis stages.

## Overview

This repository implements an end-to-end computational pipeline for:

• Harmonizing liver transcriptomic cohorts  
• Performing per-cohort differential expression analysis  
• Conducting random-effects meta-analysis across cohorts  
• Identifying a conserved NASH progression programme  
• Integrating circulating miRNA regulatory signals  
• Performing pathway enrichment and module analysis  
• Estimating liver cell-type enrichment via ssGSEA  
• Linking molecular signatures to fibrosis and disease severity  

The analysis is designed for reproducibility, manuscript-grade figure generation,
and extension to additional cohorts.

## Biological Goals

1. Identify reproducible transcriptomic changes across independent liver cohorts
2. Derive a conserved progression signature associated with fibrosis severity
3. Quantify a progression score across patients
4. Integrate serum miRNA panels with liver target networks
5. Characterize cell-type shifts during disease progression
6. Generate interpretable pathway-level and regulatory insights

## Pipeline Structure

The scripts are designed to be run sequentially.

### 00 — Phenotype Harmonization & QC
- Standardizes disease labels across cohorts
- Verifies available contrasts
- Ensesses cohort sample counts

**Output:** Harmonized phenotype object


### 01 — Expression Processing & Gene Mapping
- Builds expression matrices
- Maps ENSEMBL IDs to gene symbols
- Creates unified expression store

**Output:** Combined SYMBOL-level expression store

### 02 — Per-Cohort Differential Expression (DGE)
- Uses limma-based linear modeling
- Computes logFC, p-values, adjusted p-values
- Generates volcano plots

**Output:** Per-cohort DGE tables + figures


### 03 — Random-Effects Meta-Analysis
- Integrates cohort-level effect sizes
- Uses inverse-variance weighting
- Identifies robust cross-cohort DE genes

**Output:** Meta-analysis DEG tables + meta-volcano plots


### 04 — Pathway Enrichment Analysis
- Hallmark / KEGG / GO enrichment
- Contrast-specific pathway discovery
- Cross-contrast pathway comparison

**Output:** Enrichment tables + pathway visualizations


### 05 — Progression Programme & Core Genes
- Identifies genes conserved across severity contrasts
- Constructs progression score (PC1)
- Visualizes score vs fibrosis / NAS

**Output:** Core gene set + progression score + figure panels


### 06 — Serum miRNA Panel Integration
- Curates serum miRNA candidates
- Identifies inverse miRNA–target pairs
- Maps targets to progression modules

**Output:** miRNA-target network tables + figures


### 07 — Liver miRNA Target Network Analysis
- Integrates liver meta-DEGs with miRNA regulation
- Performs enrichment and module-level analysis
- Generates regulatory network figures

**Output:** Integrated miRNA–gene network


### 08 — Cell-Type Deconvolution
- Performs ssGSEA/GSVA-based scoring
- Estimates immune/stromal cell shifts
- Correlates cell-type scores with fibrosis

**Output:** Cell-type enrichment matrices + association plots

## How to Run

Run scripts sequentially:

```bash
Rscript 00_pheno_harmonization_qc.R
Rscript 01_expression_processing_symbol_mapping.R
Rscript 02_cohort_dge_analysis.R
Rscript 03_meta_analysis_random_effects.R
Rscript 04_pathway_enrichment_analysis.R
Rscript 05_progression_programme_core_genes.R
Rscript 06_serum_miRNA_integration.R
Rscript 07_miRNA_liver_network_analysis.R
Rscript 08_celltype_deconvolution_analysis.R
