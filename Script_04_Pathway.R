# =============================================================================
# Script 04 – Meta-pathway analysis (Hallmark only) + rich figures
# =============================================================================
# Inputs (already created):
#   Tables/Meta_Fibrosis_AdvancedVsMild_results.csv
#   Tables/Meta_NAS_HighVsLow_results.csv
#   Tables/Meta_NASH_vs_Control_results.csv
#
# Outputs:
#   Tables/Pathways_<contrast>.csv
#   Tables/Pathways_allContrasts.csv
#   Figures/MetaPathways/MetaDEG_Overlap_Heatmap_Venn.pdf
#   Figures/MetaPathways/PathwayBar_<contrast>.pdf
#   Figures/MetaPathways/PathwaySharing_Heatmap.pdf
#   Figures/MetaPathways/Alluvial_Contrast_to_Pathways.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(ComplexHeatmap)
  library(circlize)
  library(ggvenn)
  library(ggalluvial)
  library(msigdbr)
  library(clusterProfiler)
})

# ─────────────────────────────────────────────────────────────────────────────
# 0. Paths
# ─────────────────────────────────────────────────────────────────────────────
proj_root <- "D:/NASH_metaanalysis"
tab_dir   <- file.path(proj_root, "Tables")
fig_dir   <- file.path(proj_root, "Figures", "MetaPathways")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# 1. Load meta results (three contrasts)
# ─────────────────────────────────────────────────────────────────────────────
meta_files <- c(
  Fibrosis_AdvancedVsMild = "Meta_Fibrosis_AdvancedVsMild_results.csv",
  NAS_HighVsLow           = "Meta_NAS_HighVsLow_results.csv",
  NASH_vs_Control         = "Meta_NASH_vs_Control_results.csv"
)

meta_list <- lapply(names(meta_files), function(lbl) {
  f <- file.path(tab_dir, meta_files[[lbl]])
  if (!file.exists(f)) stop("Missing meta file: ", f)
  df <- read_csv(f, show_col_types = FALSE)
  
  # Standardize gene column
  if ("SYMBOL" %in% colnames(df)) {
    df <- df %>% rename(gene = SYMBOL)
  } else if (!"gene" %in% colnames(df)) {
    stop("Meta file for ", lbl, " lacks gene/SYMBOL column.")
  }
  
  # Standardize FDR and logFC column names if needed
  if ("metafdr" %in% colnames(df)) {
    df <- df %>% rename(meta_FDR = metafdr)
  }
  if ("metalogFC" %in% colnames(df)) {
    df <- df %>% rename(meta_logFC = metalogFC)
  }
  
  df$contrast <- lbl
  df
})
names(meta_list) <- names(meta_files)

# ─────────────────────────────────────────────────────────────────────────────
# 2. Define meta-DEG sets (for overlaps and enrichment)
# ─────────────────────────────────────────────────────────────────────────────
fc_cut  <- 0.2    # meta log2FC threshold
fdr_cut <- 0.05   # FDR threshold

meta_deg <- lapply(meta_list, function(df) {
  df %>%
    mutate(
      dir = case_when(
        meta_FDR < fdr_cut & meta_logFC >  fc_cut ~ "Up",
        meta_FDR < fdr_cut & meta_logFC < -fc_cut ~ "Down",
        TRUE ~ "NS"
      )
    )
})

deg_sets_up <- lapply(meta_deg, function(df) df %>% filter(dir == "Up") %>% pull(gene) %>% unique())
deg_sets_dn <- lapply(meta_deg, function(df) df %>% filter(dir == "Down") %>% pull(gene) %>% unique())

# ─────────────────────────────────────────────────────────────────────────────
# 3. Overlap heatmap + Venn (Figure-style, not volcanos)
# ─────────────────────────────────────────────────────────────────────────────
deg_counts <- tibble(
  contrast = names(meta_deg),
  Up   = sapply(deg_sets_up, length),
  Down = sapply(deg_sets_dn, length)
)

mat_counts <- as.matrix(deg_counts[, c("Up", "Down")])
rownames(mat_counts) <- deg_counts$contrast

ht <- Heatmap(
  mat_counts,
  name = "DEG count",
  col = colorRamp2(c(min(mat_counts), max(mat_counts)),
                   c("#fee0d2", "#de2d26")),
  cluster_rows = FALSE,
  cluster_columns = FALSE
)

pdf(file.path(fig_dir, "MetaDEG_Overlap_Heatmap_Venn.pdf"),
    width = 7, height = 7)
draw(ht, heatmap_legend_side = "right")

if (length(meta_deg) == 3) {
  up_list  <- deg_sets_up
  dn_list  <- deg_sets_dn
  names(up_list) <- names(meta_deg)
  names(dn_list) <- names(meta_deg)
  
  grid::grid.newpage()
  print(
    ggvenn(up_list,
           fill_color = c("#fc9272", "#9ecae1", "#a1d99b"),
           show_percentage = FALSE) +
      ggtitle("Meta-DEGs Upregulated (FDR<0.05, |log2FC|>0.2)")
  )
  
  grid::grid.newpage()
  print(
    ggvenn(dn_list,
           fill_color = c("#fc9272", "#9ecae1", "#a1d99b"),
           show_percentage = FALSE) +
      ggtitle("Meta-DEGs Downregulated (FDR<0.05, |log2FC|>0.2)")
  )
}
dev.off()

# ─────────────────────────────────────────────────────────────────────────────
# 4. Hallmark gene sets
# ─────────────────────────────────────────────────────────────────────────────
msig_h <- msigdbr(species = "Homo sapiens",
                  collection = "H") %>%
  dplyr::select(gs_name, gene_symbol)

# Background genes per contrast = all genes meta-tested
bg_genes <- lapply(meta_list, function(df) unique(df$gene))

run_enrich_for_contrast <- function(contrast_name, deg_genes, background_genes,
                                    msig_df) {
  
  if (length(deg_genes) < 10) {
    message("Too few DEGs for ", contrast_name, " – skipping enrichment.")
    return(NULL)
  }
  
  enr <- enricher(gene          = deg_genes,
                  TERM2GENE     = msig_df,
                  universe      = background_genes,
                  pAdjustMethod = "BH")
  
  if (is.null(enr) || nrow(as.data.frame(enr)) == 0) {
    message("No significant Hallmark pathways for ", contrast_name)
    return(NULL)
  }
  
  as_tibble(as.data.frame(enr)) %>%
    mutate(
      contrast   = contrast_name,
      collection = "Hallmark"
    )
}

# ─────────────────────────────────────────────────────────────────────────────
# 5. Run enrichment (Up meta-DEGs) and save tables
# ─────────────────────────────────────────────────────────────────────────────
pathway_results <- list()

for (cn in names(meta_deg)) {
  message("Enrichment for ", cn)
  up_genes <- deg_sets_up[[cn]]
  bg       <- bg_genes[[cn]]
  
  res_h <- run_enrich_for_contrast(cn, up_genes, bg, msig_h)
  
  if (!is.null(res_h) && nrow(res_h) > 0) {
    pathway_results[[cn]] <- res_h
    write_csv(res_h, file.path(tab_dir, paste0("Pathways_", cn, ".csv")))
  }
}

if (length(pathway_results) == 0) {
  stop("No pathway results produced. Check DEG thresholds or meta files.")
}

pathway_combined <- bind_rows(pathway_results)
write_csv(pathway_combined, file.path(tab_dir, "Pathways_allContrasts.csv"))

# ─────────────────────────────────────────────────────────────────────────────
# 6. Pathway barplots per contrast (Hallmark only)
# ─────────────────────────────────────────────────────────────────────────────
for (cn in names(pathway_results)) {
  df <- pathway_results[[cn]]
  
  top_h <- df %>%
    arrange(p.adjust) %>%
    slice_head(n = 20)
  
  for_plot <- top_h %>%
    mutate(
      term = str_replace(Description, "HALLMARK_", ""),
      term = str_replace_all(term, "_", " "),
      term = factor(term, levels = rev(unique(term))),
      minus_log10_FDR = -log10(p.adjust)
    )
  
  p <- ggplot(for_plot,
              aes(x = minus_log10_FDR, y = term)) +
    geom_col(fill = "#fc9272") +
    labs(
      title = paste0("Top Hallmark pathways – ", cn, " (Up meta-DEGs)"),
      x = "-log10(FDR)",
      y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0))
  
  ggsave(file.path(fig_dir, paste0("PathwayBar_", cn, ".pdf")),
         p, width = 7, height = 6)
}

# ─────────────────────────────────────────────────────────────────────────────
# 7. Pathway sharing heatmap (Hallmark)
# ─────────────────────────────────────────────────────────────────────────────
top_h_all <- pathway_combined %>%
  group_by(contrast) %>%
  arrange(p.adjust, .by_group = TRUE) %>%
  slice_head(n = 20) %>%
  ungroup()

all_terms <- unique(top_h_all$Description)
mat_term <- matrix(0,
                   nrow = length(all_terms),
                   ncol = length(meta_deg),
                   dimnames = list(all_terms, names(meta_deg)))

for (cn in names(meta_deg)) {
  df <- top_h_all %>% filter(contrast == cn)
  mat_term[df$Description, cn] <- -log10(df$p.adjust)
}

ht2 <- Heatmap(
  mat_term,
  name = "-log10(FDR)",
  col = colorRamp2(c(0, max(mat_term)), c("white", "#08519c")),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_gp = grid::gpar(fontsize = 6)
)

pdf(file.path(fig_dir, "PathwaySharing_Heatmap.pdf"),
    width = 7, height = 8)
draw(ht2, heatmap_legend_side = "right")
dev.off()

# ─────────────────────────────────────────────────────────────────────────────
# 8. Alluvial plot: contrast → Hallmark category → specific term
# ─────────────────────────────────────────────────────────────────────────────
alluvial_df <- top_h_all %>%
  mutate(
    category = str_split_fixed(Description, "_", 2)[, 1],
    term_clean = str_replace(Description, "HALLMARK_", ""),
    term_clean = str_replace_all(term_clean, "_", " "),
    minus_log10_FDR = -log10(p.adjust)
  )

if (nrow(alluvial_df) > 0) {
  p_alluvial <- ggplot(alluvial_df,
                       aes(y = minus_log10_FDR,
                           axis1 = contrast,
                           axis2 = category,
                           axis3 = term_clean)) +
    geom_alluvium(aes(fill = contrast), width = 0.2, alpha = 0.7) +
    geom_stratum(width = 0.2, fill = "grey90", color = "grey40") +
    geom_text(stat = "stratum",
              aes(label = after_stat(stratum)),
              size = 3) +
    scale_x_discrete(limits = c("Contrast", "Category", "Pathway"),
                     expand = c(0.05, 0.05)) +
    labs(
      title = "Meta-analysis: shared and contrast-specific Hallmark pathways",
      y = "Flow weight (~ -log10(FDR))",
      x = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", hjust = 0))
  
  ggsave(file.path(fig_dir, "Alluvial_Contrast_to_Pathways.pdf"),
         p_alluvial, width = 10, height = 6)
}

cat("\n✅ Script 04 complete.\n")
cat("  Tables: Pathways_<contrast>.csv, Pathways_allContrasts.csv\n")
cat("  Figures: MetaDEG_Overlap_Heatmap_Venn.pdf, PathwayBar_*.pdf,\n")
cat("           PathwaySharing_Heatmap.pdf, Alluvial_Contrast_to_Pathways.pdf\n")



# ─────────────────────────────────────────────────────────────────────────────
# 3. Overlap heatmap + Venn (nicer colours)
# ─────────────────────────────────────────────────────────────────────────────
deg_counts <- tibble(
  contrast = names(meta_deg),
  Up   = sapply(deg_sets_up, length),
  Down = sapply(deg_sets_dn, length)
)

mat_counts <- as.matrix(deg_counts[, c("Up", "Down")])
rownames(mat_counts) <- deg_counts$contrast

# Use a blue-white-red scale centered on the median count
min_val <- min(mat_counts)
max_val <- max(mat_counts)
mid_val <- median(mat_counts)

ht <- Heatmap(
  mat_counts,
  name = "DEG count",
  col  = colorRamp2(c(min_val, mid_val, max_val),
                    c("#2166ac", "white", "#b2182b")),
  cluster_rows    = FALSE,
  cluster_columns = FALSE,
  row_names_side  = "left",
  column_title    = NULL
)

pdf(file.path(fig_dir, "MetaDEG_Overlap_Heatmap_Venn.pdf"),
    width = 7, height = 7)
draw(ht, heatmap_legend_side = "right")

# Venns stay the same…
if (length(meta_deg) == 3) {
  up_list  <- deg_sets_up
  dn_list  <- deg_sets_dn
  names(up_list) <- names(meta_deg)
  names(dn_list) <- names(meta_deg)
  
  grid::grid.newpage()
  print(
    ggvenn(up_list,
           fill_color = c("#fb9a99", "#a6cee3", "#b2df8a"),
           show_percentage = FALSE) +
      ggtitle("Meta-DEGs Upregulated (FDR<0.05, |log2FC|>0.2)")
  )
  
  grid::grid.newpage()
  print(
    ggvenn(dn_list,
           fill_color = c("#fb9a99", "#a6cee3", "#b2df8a"),
           show_percentage = FALSE) +
      ggtitle("Meta-DEGs Downregulated (FDR<0.05, |log2FC|>0.2)")
  )
}
dev.off()

