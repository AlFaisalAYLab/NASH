# =============================================================================
# Script 02 – Differential Gene Expression (limma) per Cohort per Contrast
# =============================================================================
# Contrasts:
#   1. Fibrosis_AdvancedVsMild   (fibrosis_group)
#   2. NAS_HighVsLow             (nas_group)
#   3. NASH_vs_Control           (disease_harmonized)
#   4. NASH_vs_NAFL              (disease_harmonized)
#   5. Ballooning_YesVsNo        (ballooning_group)
#   6. Inflammation_YesVsNo      (inflammation_group)
# =============================================================================

library(dplyr)
library(readr)
library(stringr)
library(limma)
library(ggplot2)
library(ggrepel)

# ─────────────────────────────────────────────────────────────────────────────
# 0. PATHS
# ─────────────────────────────────────────────────────────────────────────────
proj_root  <- "D:/NASH_metaanalysis"
new_rds    <- file.path(proj_root, "Data", "RDS")
fig_dir    <- file.path(proj_root, "Figures", "DGE")
tab_dir    <- file.path(proj_root, "Tables")

if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
if (!dir.exists(tab_dir)) dir.create(tab_dir, recursive = TRUE)

# ─────────────────────────────────────────────────────────────────────────────
# 1. LOAD DATA
# ─────────────────────────────────────────────────────────────────────────────
pheno <- readRDS(file.path(new_rds, "pheno_harmonized.rds"))
expr_store <- readRDS(file.path(new_rds, "expr_store_SYMBOL_all8.rds"))

cat("Pheno:", nrow(pheno), "rows\n")
cat("Expr store:", length(expr_store), "cohorts\n")

# ─────────────────────────────────────────────────────────────────────────────
# 2. DEFINE CONTRASTS
# ─────────────────────────────────────────────────────────────────────────────
# Each contrast: name, pheno_column, group1 (numerator), group2 (denominator)
contrast_defs <- list(
  list(name = "Fibrosis_AdvancedVsMild",
       col  = "fibrosis_group",
       g1   = "Advanced",
       g2   = "Mild"),
  
  list(name = "NAS_HighVsLow",
       col  = "nas_group",
       g1   = "High_NAS",
       g2   = "Low_NAS"),
  
  list(name = "NASH_vs_Control",
       col  = "disease_harmonized",
       g1   = "NASH",
       g2   = "Control"),
  
  list(name = "NASH_vs_NAFL",
       col  = "disease_harmonized",
       g1   = "NASH",
       g2   = "NAFL"),
  
  list(name = "Ballooning_YesVsNo",
       col  = "ballooning_group",
       g1   = "Ballooning",
       g2   = "No_Ballooning"),
  
  list(name = "Inflammation_YesVsNo",
       col  = "inflammation_group",
       g1   = "Inflammation",
       g2   = "No_Inflammation")
)

# ─────────────────────────────────────────────────────────────────────────────
# 3. HELPER: Run limma for one cohort × one contrast
# ─────────────────────────────────────────────────────────────────────────────
run_limma <- function(mat, pheno_sub, contrast_def, gse,
                      min_expr = 1, min_samples_pct = 0.2) {
  
  col_name <- contrast_def$col
  g1 <- contrast_def$g1
  g2 <- contrast_def$g2
  
  # Subset pheno to the two groups
  ph <- pheno_sub %>%
    filter(.data[[col_name]] %in% c(g1, g2)) %>%
    mutate(group = factor(.data[[col_name]], levels = c(g2, g1)))
  
  if (nrow(ph) < 6 || sum(ph$group == g1) < 3 || sum(ph$group == g2) < 3) {
    return(NULL)
  }
  
  # Subset expression to matching samples
  common_samples <- intersect(colnames(mat), ph$GSM)
  if (length(common_samples) < 6) return(NULL)
  
  ph  <- ph %>% filter(GSM %in% common_samples) %>% arrange(match(GSM, common_samples))
  mat <- mat[, common_samples, drop = FALSE]
  
  # Filter low-expression genes
  min_n <- ceiling(min_samples_pct * ncol(mat))
  keep  <- rowSums(mat > min_expr) >= min_n
  mat   <- mat[keep, , drop = FALSE]
  
  if (nrow(mat) < 100) {
    cat("      Only", nrow(mat), "genes after filtering — skipping\n")
    return(NULL)
  }
  
  # Design matrix
  design <- model.matrix(~ group, data = ph)
  
  # Fit limma
  fit  <- lmFit(mat, design)
  fit  <- eBayes(fit)
  
  # Extract results for the group coefficient (g1 vs g2)
  coef_name <- colnames(design)[2]  # "groupAdvanced", "groupHigh_NAS", etc.
  tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
  
  # Add metadata
  tt$gene   <- rownames(tt)
  tt$GSE    <- gse
  tt$contrast <- contrast_def$name
  tt$n_g1   <- sum(ph$group == g1)
  tt$n_g2   <- sum(ph$group == g2)
  tt$n_total <- nrow(ph)
  tt$n_genes_tested <- nrow(mat)
  
  rownames(tt) <- NULL
  tt
}

# ─────────────────────────────────────────────────────────────────────────────
# 4. HELPER: Volcano plot
# ─────────────────────────────────────────────────────────────────────────────
make_volcano <- function(tt, gse, contrast_name,
                         fc_cut = 0.585, pval_cut = 0.05) {
  
  tt <- tt %>%
    mutate(
      sig = case_when(
        adj.P.Val < pval_cut & logFC >  fc_cut ~ "Up",
        adj.P.Val < pval_cut & logFC < -fc_cut ~ "Down",
        TRUE ~ "NS"
      ),
      sig = factor(sig, levels = c("Up", "Down", "NS"))
    )
  
  n_up   <- sum(tt$sig == "Up")
  n_down <- sum(tt$sig == "Down")
  
  # Top genes to label
  top_genes <- tt %>%
    filter(sig != "NS") %>%
    arrange(adj.P.Val) %>%
    head(20)
  
  p <- ggplot(tt, aes(x = logFC, y = -log10(adj.P.Val), color = sig)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_manual(
      values = c("Up" = "#E41A1C", "Down" = "#377EB8", "NS" = "grey70"),
      labels = c(
        paste0("Up (", n_up, ")"),
        paste0("Down (", n_down, ")"),
        "NS"
      )
    ) +
    geom_hline(yintercept = -log10(pval_cut), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = "dashed", color = "grey40") +
    geom_text_repel(
      data = top_genes,
      aes(label = gene),
      size = 2.5, max.overlaps = 15,
      color = "black", segment.color = "grey50"
    ) +
    labs(
      title    = paste0(gse, ": ", contrast_name),
      subtitle = paste0("n=", tt$n_total[1],
                        "  (", tt$n_g1[1], " vs ", tt$n_g2[1], ")"),
      x = "log2 Fold Change",
      y = "-log10(adj. P-value)",
      color = "Direction"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 12, face = "bold"))
  
  p
}

# ─────────────────────────────────────────────────────────────────────────────
# 5. MAIN LOOP – RUN ALL CONTRASTS × COHORTS
# ─────────────────────────────────────────────────────────────────────────────
all_DE_results <- list()
de_summary <- data.frame()

for (cdef in contrast_defs) {
  cat("\n", strrep("=", 60), "\n")
  cat(" CONTRAST:", cdef$name, "(", cdef$g1, "vs", cdef$g2, ")\n")
  cat(strrep("=", 60), "\n")
  
  for (gse in names(expr_store)) {
    mat       <- expr_store[[gse]]
    pheno_sub <- pheno %>% filter(GSE == gse)
    
    cat("  ", gse, ": ")
    
    tt <- tryCatch(
      run_limma(mat, pheno_sub, cdef, gse),
      error = function(e) {
        cat("ERROR:", conditionMessage(e), "\n")
        NULL
      }
    )
    
    if (is.null(tt)) {
      cat("skipped (insufficient samples)\n")
      next
    }
    
    n_sig_005 <- sum(tt$adj.P.Val < 0.05, na.rm = TRUE)
    n_sig_fc  <- sum(tt$adj.P.Val < 0.05 & abs(tt$logFC) > 0.585, na.rm = TRUE)
    
    cat(sprintf("%d genes tested | %d DEGs (FDR<0.05) | %d DEGs (FDR<0.05 & |logFC|>0.585)\n",
                nrow(tt), n_sig_005, n_sig_fc))
    
    # Store results
    key <- paste0(cdef$name, "__", gse)
    all_DE_results[[key]] <- tt
    
    # Summary row
    de_summary <- bind_rows(de_summary, tibble(
      contrast     = cdef$name,
      GSE          = gse,
      n_g1         = tt$n_g1[1],
      n_g2         = tt$n_g2[1],
      n_total      = tt$n_total[1],
      genes_tested = nrow(tt),
      DEG_FDR05    = n_sig_005,
      DEG_FDR05_FC = n_sig_fc
    ))
    
    # Volcano plot
    p <- make_volcano(tt, gse, cdef$name)
    fname <- file.path(fig_dir, paste0("Volcano_", cdef$name, "_", gse, ".pdf"))
    ggsave(fname, p, width = 7, height = 6)
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# 6. SAVE ALL DE RESULTS
# ─────────────────────────────────────────────────────────────────────────────

# Combined list
saveRDS(all_DE_results, file.path(new_rds, "DE_allContrasts_allCohorts.rds"))

# Combined table
DE_combined <- bind_rows(all_DE_results)
saveRDS(DE_combined, file.path(new_rds, "DE_combined_table.rds"))
write_csv(DE_combined, file.path(tab_dir, "DE_combined_table.csv"))

# Summary table
write_csv(de_summary, file.path(tab_dir, "DE_summary.csv"))

cat("\n========== DE SUMMARY ==========\n")
print(de_summary, n = Inf)

# ─────────────────────────────────────────────────────────────────────────────
# 7. PER-CONTRAST GENE LISTS (for meta-analysis input)
# ─────────────────────────────────────────────────────────────────────────────
cat("\n========== SAVING PER-CONTRAST DE LISTS ==========\n")

for (cname in unique(de_summary$contrast)) {
  keys <- names(all_DE_results)[grepl(paste0("^", cname, "__"), names(all_DE_results))]
  contrast_results <- all_DE_results[keys]
  names(contrast_results) <- gsub(paste0("^", cname, "__"), "", names(contrast_results))
  
  out_file <- file.path(new_rds, paste0("DE_", cname, ".rds"))
  saveRDS(contrast_results, out_file)
  cat("  Saved:", basename(out_file),
      "—", length(contrast_results), "cohorts\n")
}

# ─────────────────────────────────────────────────────────────────────────────
# 8. QUICK GENE OVERLAP CHECK ACROSS COHORTS PER CONTRAST
# ─────────────────────────────────────────────────────────────────────────────
cat("\n========== GENES TESTED PER CONTRAST ==========\n")
for (cname in unique(de_summary$contrast)) {
  keys <- names(all_DE_results)[grepl(paste0("^", cname, "__"), names(all_DE_results))]
  gene_lists <- lapply(all_DE_results[keys], function(x) x$gene)
  common <- Reduce(intersect, gene_lists)
  cohort_names <- gsub(paste0("^", cname, "__"), "", keys)
  cat(sprintf("  %-30s  cohorts=%d  common_genes=%d\n",
              cname, length(gene_lists), length(common)))
}

cat("\n✅ Script 02 complete. Outputs:\n")
cat("  RDS: DE_allContrasts_allCohorts.rds, DE_combined_table.rds, DE_<contrast>.rds\n")
cat("  CSV: DE_combined_table.csv, DE_summary.csv\n")
cat("  PDF: Volcano plots in", fig_dir, "\n")
cat("  Next: Script 03 (Meta-analysis)\n")

# Print summary
print(as.data.frame(de_summary))

# Finish the rest of Script 02
cat("\n========== SAVING PER-CONTRAST DE LISTS ==========\n")

for (cname in unique(de_summary$contrast)) {
  keys <- names(all_DE_results)[grepl(paste0("^", cname, "__"), names(all_DE_results))]
  contrast_results <- all_DE_results[keys]
  names(contrast_results) <- gsub(paste0("^", cname, "__"), "", names(contrast_results))
  
  out_file <- file.path(new_rds, paste0("DE_", cname, ".rds"))
  saveRDS(contrast_results, out_file)
  cat("  Saved:", basename(out_file), "—", length(contrast_results), "cohorts\n")
}

cat("\n========== GENES TESTED PER CONTRAST ==========\n")
for (cname in unique(de_summary$contrast)) {
  keys <- names(all_DE_results)[grepl(paste0("^", cname, "__"), names(all_DE_results))]
  gene_lists <- lapply(all_DE_results[keys], function(x) x$gene)
  common <- Reduce(intersect, gene_lists)
  cat(sprintf("  %-30s  cohorts=%d  common_genes=%d\n",
              cname, length(gene_lists), length(common)))
}

cat("\n✅ Script 02 complete.\n")
cat("  Next: Script 03 (Meta-analysis)\n")

