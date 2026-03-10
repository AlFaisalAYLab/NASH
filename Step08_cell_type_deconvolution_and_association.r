# ═══════════════════════════════════════════════════════════════════════════════
# Script_08_SingleCell_Deconvolution.R
# Cell-type deconvolution of bulk liver cohorts using published scRNA-seq
# signatures, linked to fibrosis stage and progression score
# Output: Figure 8A–8D + comprehensive tables (PDF + CSV)
# Requires:
#   - expr_store_SYMBOL_all8.rds from Script_01
#   - pheno_harmonized.rds from Script_00
#   - module_scores.rds from Script_07 (optional)
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(GSVA)
  library(patchwork)
  library(ggpubr)
  library(ggdist)
  library(corrplot)
  library(gridExtra)
  library(grid)
  library(scales)
})

select <- dplyr::select  # prevent Biobase masking

proj_root <- "D:/NASH_metaanalysis"
rds_dir   <- file.path(proj_root, "Data", "RDS")
fig_dir   <- file.path(proj_root, "Figures", "Deconvolution")
tbl_dir   <- file.path(proj_root, "Tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tbl_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n═══════════════════════════════════════════════════════════════════════════════\n")
cat("SCRIPT 08: SINGLE-CELL SIGNATURE DECONVOLUTION\n")
cat("═══════════════════════════════════════════════════════════════════════════════\n\n")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 1: Load bulk data + phenotype
# ═══════════════════════════════════════════════════════════════════════════════

cat("STEP 1: Loading bulk expression + phenotype...\n")
expr_store <- readRDS(file.path(rds_dir, "expr_store_SYMBOL_all8.rds"))
pheno <- readRDS(file.path(rds_dir, "pheno_harmonized.rds"))

cat("  Cohorts:", length(expr_store), "\n")
cat("  Samples:", nrow(pheno), "\n")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 2: Curated liver cell-type signatures
# From MacParland et al. 2018 (Nature Commun), Ramachandran et al. 2019 (Nature),
# and Aizarani et al. 2019 (Nature)
# Top 30 marker genes per cell type (scRNA-seq DEG FDR < 0.01, logFC > 1)
# ═══════════════════════════════════════════════════════════════════════════════

cat("\nSTEP 2: Building liver cell-type signature gene sets...\n")

cell_signatures <- list(
  Hepatocytes = c(
    "ALB", "APOA1", "APOB", "APOC3", "APOE", "TF", "TTR", "SERPINA1",
    "HP", "FGA", "FGB", "FGG", "CYP3A4", "CYP2E1", "CYP1A2", "ADH1B",
    "ADH4", "ALDOB", "TAT", "ASS1", "SDS", "HAL", "HGD", "HPD",
    "ASGR1", "PCK1", "G6PC", "SLC2A2", "AKR1C4", "HNF4A"
  ),
  HSC_activated = c(  # Hepatic stellate cells (activated/myofibroblast)
    "ACTA2", "COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A3",
    "TAGLN", "MYL9", "TPM2", "PDGFRB", "LOXL2", "LOX", "ELN",
    "TIMP1", "TIMP3", "BGN", "DCN", "LUM", "POSTN", "SPARC",
    "THBS2", "FN1", "VIM", "DES", "RGS5", "MMP2", "CTGF",
    "TGFB1", "TGFBI", "IGFBP3"
  ),
  HSC_quiescent = c(  # Quiescent stellate cells
    "LRAT", "RBP1", "RELN", "NGFR", "PPARg", "VIPR1", "HAND2",
    "NR2F2", "GFAP", "CXCL12", "HGF", "COLEC11", "IGFBP5", "IGFBP7",
    "CYGB", "RSPO3", "GJA1", "OLFML3", "PTCH1", "SLIT2",
    "ECSCR", "ADAMTS1", "PLAC9", "RAMP1", "VIPR2", "LHFP",
    "SYNPO2", "MEG3", "RBMS3", "PCDH7"
  ),
  Macrophages_resident = c(  # Kupffer cells
    "CD68", "CD163", "MARCO", "TIMD4", "VSIG4", "CD5L", "MERTK",
    "MRC1", "CLEC10A", "FOLR2", "SLC40A1", "HMOX1", "NR1H3",
    "SLCO2B1", "LILRB5", "SIGLEC1", "STAB1", "F13A1", "SLC11A1",
    "C1QA", "C1QB", "C1QC", "AIF1", "CSF1R", "CD14", "FCGR3A",
    "FCER1G", "LST1", "HLA-DRA", "HLA-DRB1"
  ),
  Macrophages_infiltrating = c(  # Monocyte-derived (TREM2+ scar-associated)
    "TREM2", "SPP1", "GPNMB", "CD9", "FABP5", "LGALS3", "CTSB",
    "CTSD", "CTSL", "APOE", "LPL", "LIPA", "NUPR1", "PLD3",
    "CCL2", "CCL3", "CCL4", "CXCL10", "IL1B", "TNF",
    "NFKBIA", "SOCS3", "THBS1", "VCAN", "MMP9", "S100A9",
    "S100A8", "FCN1", "EREG", "G0S2"
  ),
  Cholangiocytes = c(
    "KRT19", "KRT7", "EPCAM", "SOX9", "CFTR", "SLC4A2", "SCTR",
    "TFF1", "TFF2", "TFF3", "MUC1", "MUC6", "CLDN4", "CLDN10",
    "SPP1", "TACSTD2", "FXYD2", "DEFB1", "PIGR", "PKD2",
    "PROM1", "ANXA4", "KRT18", "KRT8", "CDH1", "MMP7",
    "HNF1B", "ONECUT1", "GGT1", "AQP1"
  ),
  Endothelial_LSEC = c(  # Liver sinusoidal endothelial cells
    "PECAM1", "CDH5", "VWF", "KDR", "FLT1", "ERG", "ENG",
    "CLEC4G", "CLEC4M", "CLEC1B", "STAB2", "FCGR2B", "OIT3",
    "LYVE1", "DNASE1L3", "FCN3", "FCN2", "MRC1", "CD36",
    "ICAM1", "PLVAP", "CALCRL", "RAMP2", "ACKR1", "MMRN1",
    "EMCN", "CD34", "CRHBP", "GJA4", "PODXL"
  ),
  T_NK_cells = c(
    "CD3D", "CD3E", "CD3G", "CD2", "CD7", "IL7R", "LEF1",
    "CD8A", "CD8B", "GZMA", "GZMB", "GZMK", "GZMH", "PRF1",
    "NKG7", "KLRD1", "KLRB1", "KLRK1", "KLRF1", "GNLY",
    "NCAM1", "IFNG", "CCL5", "CCR7", "SELL", "TCF7",
    "CD4", "FOXP3", "CTLA4", "IL2RA"
  ),
  B_Plasma_cells = c(
    "CD19", "CD79A", "CD79B", "MS4A1", "PAX5", "BANK1", "BLK",
    "CD22", "FCRL5", "IGHD", "IGHM", "IGHG1", "IGKC", "IGLC1",
    "JCHAIN", "MZB1", "XBP1", "SDC1", "PRDM1", "IRF4",
    "TNFRSF17", "CD38", "DERL3", "FKBP11", "SSR4", "SEC11C",
    "SPAG4", "TXNDC5", "BHLHA15", "CD27"
  )
)

# Check how many signature genes exist in our expression data
shared_genes <- rownames(expr_store[[1]])
sig_coverage <- sapply(cell_signatures, function(genes) {
  avail <- sum(genes %in% shared_genes)
  sprintf("%d/%d (%.0f%%)", avail, length(genes), 100 * avail / length(genes))
})
cat("\nSignature coverage in expression data:\n")
for (ct in names(sig_coverage)) cat("  ", ct, ":", sig_coverage[ct], "\n")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 3: Run GSVA (ssGSEA mode) for cell-type scoring
# ═══════════════════════════════════════════════════════════════════════════════

cat("\nSTEP 3: Running GSVA/ssGSEA cell-type scoring...\n")

deconv_list <- list()
for (gse in names(expr_store)) {
  mat <- expr_store[[gse]]

  gsva_res <- tryCatch({
    gsva(
      gsvaParam(mat, cell_signatures, kcdf = "Gaussian"),
      verbose = FALSE
    )
  }, error = function(e) {
    # Fallback: try older gsva() API
    tryCatch(
      gsva(mat, cell_signatures, method = "ssgsea", kcdf = "Gaussian", verbose = FALSE),
      error = function(e2) {
        cat("  ⚠️ GSVA failed for", gse, ":", e2$message, "\n")
        NULL
      }
    )
  })

  if (is.null(gsva_res)) next

  # Transpose: samples × cell types
  df <- as.data.frame(t(gsva_res)) %>%
    mutate(GSE = gse, GSM = rownames(t(gsva_res)))

  deconv_list[[gse]] <- df
  cat("  ", gse, ":", ncol(gsva_res), "samples scored\n")
}

deconv_all <- bind_rows(deconv_list)

# Merge with phenotype
deconv_all <- deconv_all %>%
  left_join(dplyr::select(pheno, GSE, GSM, 
                           disease = disease_harmonized,
                           fibrosis = fibrosis_stage,
                           fibrosis_grp = fibrosis_group,
                           NAS = nas_score),
            by = c("GSE", "GSM"))

cat("  Total samples with deconvolution:", nrow(deconv_all), "\n")
saveRDS(deconv_all, file.path(rds_dir, "deconv_celltype_scores.rds"))

# Cell type columns
ct_cols <- names(cell_signatures)

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 4: Link to progression score (from Script_05, if available)
# ═══════════════════════════════════════════════════════════════════════════════

cat("\nSTEP 4: Linking deconvolution to progression score...\n")

prog_file <- file.path(rds_dir, "progression_score.rds")
if (file.exists(prog_file)) {
  prog_scores <- readRDS(prog_file)
  deconv_all <- deconv_all %>%
    left_join(prog_scores, by = c("GSE", "GSM"))
  cat("  Progression scores merged\n")
} else {
  cat("  ⚠️ progression_score.rds not found — skipping correlation with ProgScore\n")
  deconv_all$ProgScore_PC1 <- NA
}

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 8A: Complex heatmap — cell-type scores × fibrosis stage
# ═══════════════════════════════════════════════════════════════════════════════

cat("\nGenerating Figure 8A: Cell-type heatmap across fibrosis stages...\n")

# Compute mean score per cell type × fibrosis stage
hm_data <- deconv_all %>%
  filter(!is.na(fibrosis)) %>%
  mutate(fibrosis = factor(fibrosis)) %>%
  group_by(fibrosis) %>%
  summarise(across(all_of(ct_cols), mean, na.rm = TRUE), .groups = "drop") %>%
  tibble::column_to_rownames("fibrosis") %>%
  as.matrix()

# Z-score per cell type (column)
hm_z <- scale(hm_data)

# Fibrosis stage annotation


fib_colors <- c("0" = "#fee0b6", "1" = "#f1a340", "2" = "#f46d43",
                "3" = "#d73027", "4" = "#a50026")
fib_stages <- rownames(hm_z)

left_ha <- rowAnnotation(
  Stage = factor(fib_stages, levels = names(fib_colors)),
  col = list(Stage = fib_colors),
  annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
  show_legend = TRUE
)


# Cell type group annotation
ct_groups <- c(
  Hepatocytes = "Parenchymal",
  HSC_activated = "Mesenchymal", HSC_quiescent = "Mesenchymal",
  Macrophages_resident = "Immune", Macrophages_infiltrating = "Immune",
  Cholangiocytes = "Parenchymal",
  Endothelial_LSEC = "Vascular",
  T_NK_cells = "Immune", B_Plasma_cells = "Immune"
)
ct_group_colors <- c(Parenchymal = "#66c2a5", Mesenchymal = "#fc8d62",
                      Immune = "#8da0cb", Vascular = "#e78ac3")

top_ha <- HeatmapAnnotation(
  Compartment = ct_groups[colnames(hm_z)],
  col = list(Compartment = ct_group_colors),
  annotation_name_gp = gpar(fontsize = 9, fontface = "bold")
)

col_fun_8a <- colorRamp2(c(-2, -1, 0, 1, 2), c("#2166ac", "#67a9cf", "white", "#ef8a62", "#b2182b"))

ht8a <- Heatmap(
  hm_z,
  name = "Z-score",
  col = col_fun_8a,
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  top_annotation = top_ha,
  left_annotation = left_ha,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_gp = gpar(fontsize = 9),
  column_names_rot = 45,
  rect_gp = gpar(col = "white", lwd = 1),
  border = TRUE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", hm_z[i, j]), x, y, gp = gpar(fontsize = 7))
  },
  row_title = "Fibrosis stage",
  row_title_gp = gpar(fontsize = 11, fontface = "bold"),
  column_title = "Cell-type enrichment scores (Z-scaled)",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(
    title = "Z-score",
    at = c(-2, -1, 0, 1, 2),
    legend_height = unit(3, "cm"),
    title_gp = gpar(fontsize = 9, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  )
)

pdf(file.path(fig_dir, "Figure_8A_CellType_Fibrosis_Heatmap.pdf"), width = 11, height = 6)
draw(ht8a, heatmap_legend_side = "right", annotation_legend_side = "right",
     padding = unit(c(10, 10, 10, 20), "mm"), merge_legend = TRUE)
dev.off()
cat("  Figure 8A saved\n")

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 8B: Raincloud plots — HSC_activated + Macrophages_infiltrating by fibrosis
# ═══════════════════════════════════════════════════════════════════════════════

cat("\nGenerating Figure 8B: Cell-type rainclouds by fibrosis...\n")

key_celltypes <- c("HSC_activated", "Macrophages_infiltrating", 
                    "Hepatocytes", "Macrophages_resident")
fib_pal <- c("0" = "#fee0b6", "1" = "#f1a340", "2" = "#f46d43",
             "3" = "#d73027", "4" = "#a50026")

rain_data <- deconv_all %>%
  filter(!is.na(fibrosis)) %>%
  mutate(fibrosis = factor(fibrosis)) %>%
  dplyr::select(GSE, GSM, fibrosis, all_of(key_celltypes)) %>%
  pivot_longer(cols = all_of(key_celltypes), names_to = "cell_type", values_to = "score") %>%
  mutate(cell_type = gsub("_", " ", cell_type))

plot_list_8b <- list()
for (ct in gsub("_", " ", key_celltypes)) {
  df_ct <- rain_data %>% filter(cell_type == ct)

  p <- ggplot(df_ct, aes(x = fibrosis, y = score, fill = fibrosis)) +
    stat_halfeye(adjust = 0.5, width = 0.6, .width = 0, justification = -0.2,
                 point_colour = NA, alpha = 0.7) +
    geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.5) +
    geom_jitter(width = 0.05, alpha = 0.2, size = 0.5, shape = 16, color = "grey30") +
    stat_compare_means(
      comparisons = list(c("0", "4"), c("1", "3"), c("0", "2")),
      method = "wilcox.test", label = "p.signif",
      tip.length = 0.01, step.increase = 0.06, size = 3
    ) +
    scale_fill_manual(values = fib_pal) +
    labs(x = "Fibrosis stage", y = "Enrichment score", title = ct) +
    theme_classic(base_size = 10) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 11),
      axis.title = element_text(face = "bold")
    )
  plot_list_8b[[ct]] <- p
}

p8b <- wrap_plots(plot_list_8b, ncol = 2) +
  plot_annotation(
    title = "Cell-type enrichment across fibrosis stages",
    subtitle = "ssGSEA scores | Wilcoxon rank-sum test",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "grey40")
    )
  )

ggsave(file.path(fig_dir, "Figure_8B_CellType_Raincloud_Fibrosis.pdf"), p8b, 
       width = 10, height = 10)
cat("  Figure 8B saved\n")

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 8C: Correlation matrix — cell types vs disease features
# ═══════════════════════════════════════════════════════════════════════════════

cat("\nGenerating Figure 8C: Correlation matrix...\n")

cor_data <- deconv_all %>%
  filter(!is.na(fibrosis)) %>%
  mutate(
    fibrosis_num = as.numeric(as.character(fibrosis)),
    NAS_num = as.numeric(as.character(NAS))
  )

# Spearman correlation: cell types vs fibrosis + NAS
cor_mat <- matrix(NA, nrow = length(ct_cols), ncol = 2,
                  dimnames = list(ct_cols, c("Fibrosis", "NAS")))
pval_mat <- cor_mat

for (ct in ct_cols) {
  if (sum(!is.na(cor_data[[ct]]) & !is.na(cor_data$fibrosis_num)) > 10) {
    res_f <- cor.test(cor_data[[ct]], cor_data$fibrosis_num, method = "spearman", exact = FALSE)
    cor_mat[ct, "Fibrosis"] <- res_f$estimate
    pval_mat[ct, "Fibrosis"] <- res_f$p.value
  }
  if (sum(!is.na(cor_data[[ct]]) & !is.na(cor_data$NAS_num)) > 10) {
    res_n <- cor.test(cor_data[[ct]], cor_data$NAS_num, method = "spearman", exact = FALSE)
    cor_mat[ct, "NAS"] <- res_n$estimate
    pval_mat[ct, "NAS"] <- res_n$p.value
  }
}

# Stars for significance
star_mat <- matrix("", nrow = nrow(pval_mat), ncol = ncol(pval_mat))
star_mat[pval_mat < 0.05] <- "*"
star_mat[pval_mat < 0.01] <- "**"
star_mat[pval_mat < 0.001] <- "***"

# ComplexHeatmap version
col_fun_cor <- colorRamp2(c(-0.5, -0.25, 0, 0.25, 0.5), 
                           c("#2166ac", "#92c5de", "white", "#f4a582", "#b2182b"))

ht8c <- Heatmap(
  cor_mat,
  name = "Spearman rho",
  col = col_fun_cor,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 11, fontface = "bold"),
  column_names_rot = 0,
  column_names_centered = TRUE,
  rect_gp = gpar(col = "white", lwd = 2),
  border = TRUE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (!is.na(cor_mat[i, j])) {
      grid.text(sprintf("%.2f", cor_mat[i, j]), x, y - unit(1, "mm"),
                gp = gpar(fontsize = 9, fontface = "bold"))
      grid.text(star_mat[i, j], x, y + unit(3, "mm"),
                gp = gpar(fontsize = 11, col = "black"))
    }
  },
  column_title = "Cell-type correlation with disease severity",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(
    title = "Spearman\nrho",
    at = c(-0.5, -0.25, 0, 0.25, 0.5),
    legend_height = unit(3, "cm"),
    title_gp = gpar(fontsize = 9, fontface = "bold")
  ),
  left_annotation = rowAnnotation(
    Compartment = ct_groups[ct_cols],
    col = list(Compartment = ct_group_colors),
    annotation_name_gp = gpar(fontsize = 9)
  ),
  width = unit(5, "cm")
)

pdf(file.path(fig_dir, "Figure_8C_CellType_Correlation_Matrix.pdf"), width = 8, height = 7)
draw(ht8c, heatmap_legend_side = "right", annotation_legend_side = "right",
     padding = unit(c(10, 10, 10, 20), "mm"), merge_legend = TRUE)
dev.off()
cat("  Figure 8C saved\n")

# ──────────────────────────────────────────────────────────────────────────────
# FIGURE 8D: Per-cohort trend lines — HSC_activated + Macrophages_infiltrating
# ──────────────────────────────────────────────────────────────────────────────
# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 8D: Per-cohort trend lines — HSC_activated + Macrophages_infiltrating
# ══════════════════════════════════════════════════════════════════════════════

cat("\nGenerating Figure 8D: Per-cohort trend lines...\n")

cohort_pal <- c(
  GSE135251 = "#E41A1C", GSE162694 = "#377EB8", GSE130970 = "#4DAF4A",
  GSE49541  = "#984EA3", GSE89632  = "#FF7F00", GSE48452  = "#A65628",
  GSE126848 = "#F781BF", GSE167523 = "#1B9E77"
)

trend_data <- deconv_all %>%
  filter(!is.na(fibrosis)) %>%
  mutate(fib_num = as.numeric(gsub("[^0-9]", "", as.character(fibrosis))))

cat("  fib_num range:", range(trend_data$fib_num, na.rm = TRUE), "\n")
cat("  Non-NA fib_num:", sum(!is.na(trend_data$fib_num)), "/", nrow(trend_data), "\n")

make_trend <- function(ct, ct_label) {
  df <- trend_data[!is.na(trend_data[[ct]]) & !is.na(trend_data$fib_num), , drop = FALSE]
  df$score <- df[[ct]]
  
  sp <- df %>%
    group_by(GSE) %>%
    summarise(
      n_pairs = sum(!is.na(fib_num) & !is.na(score)),
      rho = if (n_pairs > 5)
        cor(fib_num, score, method = "spearman", use = "pair")
      else NA_real_,
      p = if (n_pairs > 5)
        cor.test(fib_num, score, method = "spearman", exact = FALSE)$p.value
      else NA_real_,
      .groups = "drop"
    ) %>%
    filter(!is.na(rho)) %>%
    mutate(
      stars = ifelse(p < 0.001, "***",
                     ifelse(p < 0.01,  "**",
                            ifelse(p < 0.05,  "*", "ns"))),
      label = paste0(GSE, "  rho=", round(rho, 2), " ", stars)
    )
  
  if (nrow(sp) == 0) {
    warning("No cohorts with valid pairs for ", ct_label)
    return(ggplot() + theme_void() + ggtitle(paste("No data for", ct_label)))
  }
  
  clabels <- setNames(sp$label, sp$GSE)
  
  ggplot(df[df$GSE %in% sp$GSE, ],
         aes(x = fib_num, y = score, color = GSE, group = GSE)) +
    stat_summary(fun = mean, geom = "line", linewidth = 1.2, alpha = 0.9) +
    stat_summary(fun = mean, geom = "point", size = 2.5) +
    stat_summary(fun.data = mean_se, geom = "errorbar",
                 width = 0.15, linewidth = 0.6, alpha = 0.7) +
    scale_x_continuous(breaks = 0:4, labels = paste0("F", 0:4)) +
    scale_color_manual(values = cohort_pal, labels = clabels) +
    labs(x = "Fibrosis stage",
         y = paste0(ct_label, " score (mean +/- SE)"),
         color = NULL) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 7, family = "mono"),
      axis.title = element_text(face = "bold")
    )
}

p8d_hsc <- make_trend("HSC_activated", "Activated HSC") +
  ggtitle("Activated HSC enrichment") +
  theme(plot.title = element_text(face = "bold", size = 12))

p8d_mac <- make_trend("Macrophages_infiltrating", "Infiltrating macrophages") +
  ggtitle("Infiltrating macrophage enrichment") +
  theme(plot.title = element_text(face = "bold", size = 12))

p8d_hep <- make_trend("Hepatocytes", "Hepatocytes") +
  ggtitle("Hepatocyte enrichment") +
  theme(plot.title = element_text(face = "bold", size = 12))

p8d_kup <- make_trend("Macrophages_resident", "Kupffer cells") +
  ggtitle("Kupffer cell enrichment") +
  theme(plot.title = element_text(face = "bold", size = 12))

p8d <- (p8d_hsc | p8d_mac) / (p8d_hep | p8d_kup) +
  plot_annotation(
    title = "Cell-type enrichment trends across fibrosis stages",
    subtitle = "Per-cohort mean +/- SE | Spearman rho per cohort",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "grey40")
    )
  )

ggsave(file.path(fig_dir, "Figure_8D_CellType_Fibrosis_Trends.pdf"),
       p8d, width = 14, height = 10)
cat("  Figure 8D saved\n")


# ═══════════════════════════════════════════════════════════════════════════════
# TABLE 1 (PDF): Cell-type deconvolution summary statistics
# ═══════════════════════════════════════════════════════════════════════════════

cat("\nGenerating comprehensive PDF tables...\n")

# Table: mean score per fibrosis stage × cell type
summary_table <- deconv_all %>%
  filter(!is.na(fibrosis)) %>%
  mutate(fibrosis = factor(fibrosis)) %>%
  group_by(fibrosis) %>%
  summarise(
    N = n(),
    across(all_of(ct_cols), ~ sprintf("%.3f (%.3f)", mean(.x, na.rm = TRUE), sd(.x, na.rm = TRUE))),
    .groups = "drop"
  )

# Pretty column names
pretty_names <- c("Fibrosis\nStage", "N",
                  "Hepato-\ncytes", "HSC\n(activated)", "HSC\n(quiescent)",
                  "Kupffer\ncells", "Infiltrating\nmacrophages", "Cholangio-\ncytes",
                  "LSEC", "T/NK\ncells", "B/Plasma\ncells")

pdf(file.path(tbl_dir, "Table_CellType_Deconvolution_Summary.pdf"), width = 16, height = 6)
grid.newpage()

title_grob <- textGrob("Cell-Type Enrichment Scores by Fibrosis Stage",
                        gp = gpar(fontsize = 14, fontface = "bold"))
subtitle_grob <- textGrob("Values shown as mean (SD) | ssGSEA enrichment scores",
                           gp = gpar(fontsize = 10, col = "grey40"))

tg <- tableGrob(
  summary_table,
  rows = NULL,
  theme = ttheme_default(
    core = list(
      fg_params = list(cex = 0.65, fontfamily = "mono"),
      bg_params = list(fill = c("#f7f7f7", "white"))
    ),
    colhead = list(
      fg_params = list(cex = 0.7, fontface = "bold"),
      bg_params = list(fill = "#d9e2f3")
    )
  )
)

grid.arrange(title_grob, subtitle_grob, tg, ncol = 1, 
             heights = c(0.06, 0.04, 0.90))
dev.off()

write.csv(summary_table, file.path(tbl_dir, "Table_CellType_Deconvolution_Summary.csv"), 
          row.names = FALSE)
cat("  Table 1 saved (PDF + CSV)\n")

# ═══════════════════════════════════════════════════════════════════════════════
# TABLE 2 (PDF): Correlation results (cell type × disease features)
# ═══════════════════════════════════════════════════════════════════════════════

cor_table <- data.frame(
  Cell_Type = rownames(cor_mat),
  Compartment = ct_groups[rownames(cor_mat)],
  Fibrosis_rho = sprintf("%.3f", cor_mat[, "Fibrosis"]),
  Fibrosis_p = sprintf("%.1e", pval_mat[, "Fibrosis"]),
  Fibrosis_sig = star_mat[, 1],
  NAS_rho = sprintf("%.3f", cor_mat[, "NAS"]),
  NAS_p = sprintf("%.1e", pval_mat[, "NAS"]),
  NAS_sig = star_mat[, 2],
  stringsAsFactors = FALSE
)

pdf(file.path(tbl_dir, "Table_CellType_Correlation_Results.pdf"), width = 12, height = 5)
grid.newpage()

title_grob2 <- textGrob("Cell-Type vs Disease Severity: Spearman Correlations",
                          gp = gpar(fontsize = 14, fontface = "bold"))
subtitle_grob2 <- textGrob("* p<0.05  ** p<0.01  *** p<0.001 | All 8 cohorts pooled",
                             gp = gpar(fontsize = 10, col = "grey40"))

tg2 <- tableGrob(
  cor_table,
  rows = NULL,
  theme = ttheme_default(
    core = list(
      fg_params = list(cex = 0.7, fontfamily = "mono"),
      bg_params = list(fill = c("#f7f7f7", "white"))
    ),
    colhead = list(
      fg_params = list(cex = 0.75, fontface = "bold"),
      bg_params = list(fill = "#fde0dd")
    )
  )
)

grid.arrange(title_grob2, subtitle_grob2, tg2, ncol = 1,
             heights = c(0.08, 0.05, 0.87))
dev.off()

write.csv(cor_table, file.path(tbl_dir, "Table_CellType_Correlation_Results.csv"),
          row.names = FALSE)
cat("  Table 2 saved (PDF + CSV)\n")

# ═══════════════════════════════════════════════════════════════════════════════
# TABLE 3 (PDF): Signature gene sets used
# ═══════════════════════════════════════════════════════════════════════════════

sig_table <- data.frame(
  Cell_Type = rep(names(cell_signatures), sapply(cell_signatures, length)),
  Gene = unlist(cell_signatures),
  stringsAsFactors = FALSE
) %>%
  group_by(Cell_Type) %>%
  summarise(
    N_genes = n(),
    N_in_data = sum(Gene %in% shared_genes),
    Coverage = sprintf("%.0f%%", 100 * N_in_data / N_genes),
    Top_markers = paste(head(Gene[Gene %in% shared_genes], 8), collapse = ", "),
    .groups = "drop"
  )

pdf(file.path(tbl_dir, "Table_CellType_Signature_Genes.pdf"), width = 14, height = 5)
grid.newpage()

title_grob3 <- textGrob("Cell-Type Signature Gene Sets",
                          gp = gpar(fontsize = 14, fontface = "bold"))
subtitle_grob3 <- textGrob("Curated from MacParland 2018, Ramachandran 2019, Aizarani 2019",
                             gp = gpar(fontsize = 10, col = "grey40"))

tg3 <- tableGrob(
  sig_table,
  rows = NULL,
  theme = ttheme_default(
    core = list(
      fg_params = list(cex = 0.7),
      bg_params = list(fill = c("#f7f7f7", "white"))
    ),
    colhead = list(
      fg_params = list(cex = 0.75, fontface = "bold"),
      bg_params = list(fill = "#d4edda")
    )
  )
)

grid.arrange(title_grob3, subtitle_grob3, tg3, ncol = 1,
             heights = c(0.08, 0.05, 0.87))
dev.off()

write.csv(sig_table, file.path(tbl_dir, "Table_CellType_Signature_Genes.csv"),
          row.names = FALSE)
cat("  Table 3 saved (PDF + CSV)\n")

# ═══════════════════════════════════════════════════════════════════════════════
# SUMMARY
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n═══════════════════════════════════════════════════════════════════════════════\n")
cat("SCRIPT 08 COMPLETE — SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════════════════════\n\n")
cat("Cell-type signatures: ", length(cell_signatures), " types\n")
cat("Samples scored:       ", nrow(deconv_all), "\n")
cat("\nFigures saved to:", fig_dir, "\n")
cat("  8A — Cell-type × fibrosis heatmap (annotated, Z-scored, with values)\n")
cat("  8B — Raincloud plots (4 key cell types × fibrosis)\n")
cat("  8C — Spearman correlation matrix (cell type × fibrosis + NAS)\n")
cat("  8D — Per-cohort trend lines (HSC + macrophages × fibrosis)\n")
cat("\nTables saved to:", tbl_dir, "\n")
cat("  Table_CellType_Deconvolution_Summary (PDF + CSV)\n")
cat("  Table_CellType_Correlation_Results (PDF + CSV)\n")
cat("  Table_CellType_Signature_Genes (PDF + CSV)\n")
cat("\n🎉 ALL 8 SCRIPTS COMPLETE — Manuscript figures and tables ready!\n")
cat("\nFinal manuscript figure lineup:\n")
cat("  Fig 1: Study design (BioRender/Illustrator)\n")
cat("  Fig 2: DGE volcanos (Script_02)\n")
cat("  Fig 3: Meta-analysis volcanos + overlap (Script_03)\n")
cat("  Fig 4: Pathway enrichment + alluvial (Script_04)\n")
cat("  Fig 5: Progression programme (Script_05)\n")
cat("  Fig 6: Serum miRNA panel + integration (Script_06)\n")
cat("  Fig 7: miRNA-liver target network (Script_07)\n")
cat("  Fig 8: Cell-type deconvolution (Script_08)\n")

