# =============================================================================
# Script 03 – Meta-analysis of DGE
# =============================================================================
# Start with: Fibrosis_AdvancedVsMild, NAS_HighVsLow, NASH_vs_Control
# =============================================================================

library(dplyr)
library(metafor)
library(ggplot2)
library(ggrepel)
library(readr)

proj_root <- "D:/NASH_metaanalysis"
rds_dir   <- file.path(proj_root, "Data", "RDS")
fig_dir   <- file.path(proj_root, "Figures", "Meta")
tab_dir   <- file.path(proj_root, "Tables")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

# Which contrasts to meta-analyze
meta_contrasts <- c("Fibrosis_AdvancedVsMild",
                    "NAS_HighVsLow",
                    "NASH_vs_Control")

# Helper: build meta input table for one contrast
build_meta_input <- function(cname) {
  de_list <- readRDS(file.path(rds_dir, paste0("DE_", cname, ".rds")))
  cohorts <- names(de_list)
  
  # Intersect genes across all cohorts
  gene_lists <- lapply(de_list, function(x) x$gene)
  common_genes <- Reduce(intersect, gene_lists)
  
  cat("\n", cname, ": cohorts =", length(cohorts),
      "| common genes =", length(common_genes), "\n")
  
  # Long table: one row per gene × cohort
  meta_tab <- bind_rows(lapply(cohorts, function(gse) {
    tt <- de_list[[gse]]
    tt <- tt[match(common_genes, tt$gene), ]
    stopifnot(all(tt$gene == common_genes))
    
    # limma provides logFC and standard error (from t and logFC)
    se <- abs(tt$logFC / tt$t)  # se = |logFC / t|
    
    tibble(
      gene   = tt$gene,
      GSE    = gse,
      logFC  = tt$logFC,
      se     = se,
      pval   = tt$P.Value
    )
  }))
  
  meta_tab
}

# Helper: run random-effects meta-analysis per gene
run_meta <- function(meta_tab, cname) {
  genes <- unique(meta_tab$gene)
  
  res_list <- lapply(genes, function(g) {
    df <- meta_tab %>% filter(gene == g)
    
    # metafor random-effects on logFC with known SE
    tryCatch({
      m <- rma(yi = logFC, sei = se, method = "REML", data = df)
      
      tibble(
        gene      = g,
        k         = nrow(df),
        meta_logFC = as.numeric(m$b),
        meta_se    = m$se,
        meta_z     = m$zval,
        meta_pval  = m$pval,
        Q          = m$QE,
        Q_pval     = m$QEp,
        I2         = m$I2
      )
    }, error = function(e) {
      NULL
    })
  })
  
  res <- bind_rows(res_list)
  res <- res %>%
    mutate(
      meta_FDR = p.adjust(meta_pval, method = "BH"),
      contrast = cname
    )
  
  res
}

# Helper: volcano for meta-analysis
plot_meta_volcano <- function(meta_res, cname,
                              fc_cut = 0.2, p_cut = 0.05) {
  df <- meta_res %>%
    mutate(
      sig = case_when(
        meta_FDR < p_cut & meta_logFC >  fc_cut ~ "Up",
        meta_FDR < p_cut & meta_logFC < -fc_cut ~ "Down",
        TRUE ~ "NS"
      ),
      sig = factor(sig, levels = c("Up", "Down", "NS"))
    )
  
  n_up   <- sum(df$sig == "Up")
  n_down <- sum(df$sig == "Down")
  
  top_genes <- df %>%
    filter(sig != "NS") %>%
    arrange(meta_FDR) %>%
    head(20)
  
  ggplot(df, aes(x = meta_logFC, y = -log10(meta_FDR), color = sig)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_manual(
      values = c("Up" = "#E41A1C", "Down" = "#377EB8", "NS" = "grey70"),
      labels = c(
        paste0("Up (", n_up, ")"),
        paste0("Down (", n_down, ")"),
        "NS"
      )
    ) +
    geom_hline(yintercept = -log10(p_cut), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = "dashed", color = "grey40") +
    geom_text_repel(
      data = top_genes,
      aes(label = gene),
      size = 2.5, max.overlaps = 20
    ) +
    labs(
      title = paste0("Meta-analysis: ", cname),
      x = "Meta log2FC (random-effects)",
      y = "-log10(FDR)",
      color = "Direction"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
}

# ─────────────────────────────────────────────────────────────────────────────
# MAIN LOOP
# ─────────────────────────────────────────────────────────────────────────────
all_meta <- list()

for (cname in meta_contrasts) {
  cat("\n============================================\n")
  cat("META-ANALYSIS:", cname, "\n")
  cat("============================================\n")
  
  meta_tab <- build_meta_input(cname)
  meta_res <- run_meta(meta_tab, cname)
  
  cat("  Meta-tested genes:", nrow(meta_res), "\n")
  cat("  Significant at FDR<0.05:", sum(meta_res$meta_FDR < 0.05, na.rm = TRUE), "\n")
  
  all_meta[[cname]] <- meta_res
  
  # Save tables
  write_csv(meta_res,
            file.path(tab_dir, paste0("Meta_", cname, "_results.csv")))
  saveRDS(meta_res,
          file.path(rds_dir, paste0("Meta_", cname, "_results.rds")))
  
  # Volcano
  p <- plot_meta_volcano(meta_res, cname)
  ggsave(file.path(fig_dir, paste0("MetaVolcano_", cname, ".pdf")),
         p, width = 7, height = 6)
}

# Combined meta table
meta_combined <- bind_rows(all_meta)
write_csv(meta_combined, file.path(tab_dir, "Meta_allContrasts_results.csv"))
saveRDS(all_meta, file.path(rds_dir, "Meta_allContrasts_results.rds"))

cat("\n✅ Script 03 complete.\n")
cat("  - Per-contrast meta tables in Tables/Meta_*_results.csv\n")
cat("  - Meta volcano plots in Figures/Meta\n")
