# =============================================================================
# Script 01 – Build expr_store_SYMBOL_all8.rds
# =============================================================================
# 1) Rebuild Ensembl→SYMBOL annotation (fresh from org.Hs.eg.db)
# 2) Rebuild RNA-seq EXPR_SYMBOL files with clean annotation
# 3) Load all 8 pre-built EXPR_SYMBOL files
# 4) Align columns to pheno GSM IDs (via eset pData where needed)
# 5) Log2-transform raw counts
# 6) Save + QC
# =============================================================================
# ═══════════════════════════════════════════════════════════════════════════
# FIX GSE126848: Ensembl → SYMBOL with positional GSM assignment
# ═══════════════════════════════════════════════════════════════════════════
cat("=== Rebuilding GSE126848 with positional GSM mapping ===\n")

ens_to_sym_fresh <- readRDS(file.path(old_rds, "Ensembl_to_HGNC_Annotation.rds"))

# Use raw expr (has ENSG rownames, 57 cols)
expr_126_raw <- readRDS(file.path(geo_cache, "GSE126848_expr.rds"))

# Step 1: Assign GSM IDs positionally (57 = 57)
pd_126 <- pData(eset_126)
colnames(expr_126_raw) <- rownames(pd_126)  # GSM IDs
cat("Columns now:", head(colnames(expr_126_raw), 3), "\n")

# Step 2: Ensembl → SYMBOL
ens_ids <- gsub("\\.\\d+$", "", rownames(expr_126_raw))

# Check overlap
cat("Ensembl overlap:", sum(ens_ids %in% ens_to_sym_fresh$ENSEMBL), 
    "/", length(ens_ids), "\n")

df <- data.frame(ENSEMBL = ens_ids, expr_126_raw, check.names = FALSE)
merged <- merge(ens_to_sym_fresh, df, by = "ENSEMBL")
merged$ENSEMBL <- NULL
cat("After merge:", nrow(merged), "rows\n")

sym_df <- aggregate(. ~ SYMBOL, data = merged, FUN = mean, na.rm = TRUE)
rownames(sym_df) <- sym_df$SYMBOL
sym_df$SYMBOL <- NULL
mat_126_sym <- as.matrix(sym_df)

cat("GSE126848 final:", dim(mat_126_sym), "\n")
cat("Colnames:", head(colnames(mat_126_sym), 3), "\n")
cat("Genes:", head(rownames(mat_126_sym), 5), "\n")
cat("Range:", range(mat_126_sym, na.rm = TRUE), "\n")

saveRDS(mat_126_sym, file.path(old_rds, "GSE126848_EXPR_SYMBOL.rds"))
cat("✅ GSE126848 fixed\n")

library(dplyr)
library(readr)
library(stringr)
library(tibble)
library(Biobase)
library(AnnotationDbi)
library(org.Hs.eg.db)

# ─────────────────────────────────────────────────────────────────────────────
# 0. PATHS
# ─────────────────────────────────────────────────────────────────────────────
proj_root  <- "D:/NASH_metaanalysis"
geo_cache  <- file.path(proj_root, "GEO_cache")
old_rds    <- file.path(proj_root, "old_analysis", "NASH_Outputs", "RDS_Objects")
new_rds    <- file.path(proj_root, "Data", "RDS")
new_idmaps <- file.path(proj_root, "Data", "ID_maps")

if (!dir.exists(new_rds))    dir.create(new_rds, recursive = TRUE)
if (!dir.exists(new_idmaps)) dir.create(new_idmaps, recursive = TRUE)

# ─────────────────────────────────────────────────────────────────────────────
# 1. LOAD PHENO
# ─────────────────────────────────────────────────────────────────────────────
pheno <- readRDS(file.path(new_rds, "pheno_harmonized.rds"))
cat("Pheno loaded:", nrow(pheno), "rows\n")

# ─────────────────────────────────────────────────────────────────────────────
# 2. REBUILD ENSEMBL → SYMBOL ANNOTATION (fresh, no Excel corruption)
# ─────────────────────────────────────────────────────────────────────────────
cat("\n========== REBUILDING ENSEMBL ANNOTATION ==========\n")

ens_to_sym_fresh <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = AnnotationDbi::keys(org.Hs.eg.db, keytype = "ENSEMBL"),
  columns = c("ENSEMBL", "SYMBOL"),
  keytype = "ENSEMBL"
) %>%
  filter(!is.na(SYMBOL), SYMBOL != "") %>%
  distinct(ENSEMBL, .keep_all = TRUE)

cat("Fresh Ensembl annotation:", nrow(ens_to_sym_fresh), "mappings\n")
saveRDS(ens_to_sym_fresh, file.path(old_rds, "Ensembl_to_HGNC_Annotation.rds"))

# ─────────────────────────────────────────────────────────────────────────────
# 3. REBUILD RNA-SEQ EXPR_SYMBOL FILES (Ensembl cohorts)
# ─────────────────────────────────────────────────────────────────────────────
cat("\n========== REBUILDING ENSEMBL COHORT SYMBOL FILES ==========\n")

ens_cohorts <- c("GSE135251", "GSE162694", "GSE126848")

for (gse in ens_cohorts) {
  expr_raw <- readRDS(file.path(geo_cache, paste0(gse, "_expr.rds")))
  ens_ids  <- gsub("\\.\\d+$", "", rownames(expr_raw))
  df       <- data.frame(ENSEMBL = ens_ids, expr_raw, check.names = FALSE)
  merged   <- merge(ens_to_sym_fresh, df, by = "ENSEMBL")
  merged$ENSEMBL <- NULL
  
  sym_df <- aggregate(. ~ SYMBOL, data = merged, FUN = mean, na.rm = TRUE)
  rownames(sym_df) <- sym_df$SYMBOL
  sym_df$SYMBOL <- NULL
  mat <- as.matrix(sym_df)
  
  cat(sprintf("  %s: %d genes × %d samples  range: %.1f–%.1f\n",
              gse, nrow(mat), ncol(mat),
              min(mat, na.rm = TRUE), max(mat, na.rm = TRUE)))
  saveRDS(mat, file.path(old_rds, paste0(gse, "_EXPR_SYMBOL.rds")))
}

# ─────────────────────────────────────────────────────────────────────────────
# 4. LOAD ALL 8 PRE-BUILT EXPR_SYMBOL FILES
# ─────────────────────────────────────────────────────────────────────────────
cat("\n========== LOADING EXPR_SYMBOL FILES ==========\n")

gse_ids <- c("GSE135251","GSE162694","GSE130970","GSE126848",
             "GSE167523","GSE49541","GSE89632","GSE48452")

expr_store_SYMBOL <- list()

for (gse in gse_ids) {
  f <- file.path(old_rds, paste0(gse, "_EXPR_SYMBOL.rds"))
  if (file.exists(f)) {
    expr_store_SYMBOL[[gse]] <- readRDS(f)
    cat(sprintf("  Loaded %-12s  %d genes × %d samples\n",
                gse, nrow(expr_store_SYMBOL[[gse]]),
                ncol(expr_store_SYMBOL[[gse]])))
  } else {
    cat("  MISSING:", f, "\n")
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# 5. BUILD COLUMN → GSM MAPS FROM ESET pData
# ─────────────────────────────────────────────────────────────────────────────
cat("\n========== LOADING ESET pData FOR COLUMN MAPPING ==========\n")

eset_maps <- list()
for (gse in gse_ids) {
  eset_file <- file.path(geo_cache, paste0(gse, "_eset.rds"))
  if (file.exists(eset_file)) {
    eset <- readRDS(eset_file)
    if (is(eset, "ExpressionSet")) {
      pd <- pData(eset)
      eset_maps[[gse]] <- pd
      cat(sprintf("  %-12s  eset pData: %d rows\n", gse, nrow(pd)))
    }
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# 6. ALIGN COLUMNS TO PHENO GSM IDs
# ─────────────────────────────────────────────────────────────────────────────
align_columns <- function(mat, pheno_sub, gse, eset_pd = NULL) {
  ecols  <- colnames(mat)
  gsms   <- pheno_sub$GSM
  titles <- pheno_sub$sample_title
  
  gsm_hit   <- length(intersect(ecols, gsms))
  title_hit <- length(intersect(ecols, titles))
  
  # Case 1: columns are already GSM IDs
  if (gsm_hit >= 0.5 * length(gsms)) {
    keep <- intersect(gsms, ecols)
    cat("    GSM match → keeping", length(keep), "of", length(gsms), "\n")
    return(mat[, keep, drop = FALSE])
  }
  
  # Case 2: columns match sample_titles exactly
  if (title_hit >= 0.5 * length(titles)) {
    keep <- intersect(titles, ecols)
    mat_sub <- mat[, keep, drop = FALSE]
    name_map <- setNames(pheno_sub$GSM, pheno_sub$sample_title)
    colnames(mat_sub) <- name_map[colnames(mat_sub)]
    cat("    Title→GSM remap → keeping", ncol(mat_sub), "\n")
    return(mat_sub)
  }
  
  # Case 3: columns match last token of sample_title
  title_tokens <- str_extract(titles, "\\S+$")
  token_hit <- length(intersect(ecols, title_tokens))
  if (token_hit >= 0.5 * length(titles)) {
    keep <- intersect(title_tokens, ecols)
    mat_sub <- mat[, keep, drop = FALSE]
    name_map <- setNames(pheno_sub$GSM, title_tokens)
    colnames(mat_sub) <- name_map[colnames(mat_sub)]
    cat("    Title-token→GSM remap → keeping", ncol(mat_sub), "\n")
    return(mat_sub)
  }
  
  # Case 4: use eset pData to find column → GSM mapping
  if (!is.null(eset_pd)) {
    pd_gsms <- rownames(eset_pd)
    for (col_name in colnames(eset_pd)) {
      pd_vals <- as.character(eset_pd[[col_name]])
      overlap <- length(intersect(ecols, pd_vals))
      if (overlap >= 0.5 * ncol(mat)) {
        map_df <- data.frame(
          expr_col = pd_vals,
          GSM      = pd_gsms,
          stringsAsFactors = FALSE
        ) %>% filter(expr_col %in% ecols, GSM %in% gsms)
        
        keep_cols <- map_df$expr_col
        mat_sub <- mat[, keep_cols, drop = FALSE]
        colnames(mat_sub) <- map_df$GSM
        cat("    eset pData col '", col_name, "' → keeping", ncol(mat_sub), "\n")
        return(mat_sub)
      }
    }
    
    # Try cleaned pData values (last word)
    if ("geo_accession" %in% colnames(eset_pd)) {
      geo_acc <- as.character(eset_pd$geo_accession)
      for (col_name in colnames(eset_pd)) {
        pd_vals_clean <- str_extract(as.character(eset_pd[[col_name]]), "\\S+$")
        overlap <- length(intersect(ecols, pd_vals_clean))
        if (overlap >= 0.5 * ncol(mat)) {
          map_df <- data.frame(
            expr_col = pd_vals_clean,
            GSM      = geo_acc,
            stringsAsFactors = FALSE
          ) %>% filter(expr_col %in% ecols, GSM %in% gsms)
          
          keep_cols <- map_df$expr_col
          mat_sub <- mat[, keep_cols, drop = FALSE]
          colnames(mat_sub) <- map_df$GSM
          cat("    eset pData cleaned col '", col_name, "' → keeping",
              ncol(mat_sub), "\n")
          return(mat_sub)
        }
      }
    }
  }
  
  # Case 5: same count → assume order matches
  if (ncol(mat) == nrow(pheno_sub)) {
    cat("    WARNING: no name overlap — assuming column order matches pheno\n")
    colnames(mat) <- pheno_sub$GSM
    return(mat)
  }
  
  cat("    ERROR: cannot align columns (expr=", ncol(mat),
      ", pheno=", nrow(pheno_sub), ")\n")
  NULL
}

cat("\n========== ALIGNING COLUMNS ==========\n")
for (gse in names(expr_store_SYMBOL)) {
  cat("  ", gse, ":\n")
  ph_sub  <- pheno %>% filter(GSE == gse)
  eset_pd <- eset_maps[[gse]]
  mat_aligned <- align_columns(expr_store_SYMBOL[[gse]], ph_sub, gse, eset_pd)
  
  if (!is.null(mat_aligned)) {
    expr_store_SYMBOL[[gse]] <- mat_aligned
  } else {
    cat("    FAILED — removing from store\n")
    expr_store_SYMBOL[[gse]] <- NULL
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# 7. LOG2 TRANSFORM RAW COUNTS WHERE NEEDED
# ─────────────────────────────────────────────────────────────────────────────
cat("\n========== LOG2 TRANSFORM CHECK ==========\n")
for (gse in names(expr_store_SYMBOL)) {
  mat <- expr_store_SYMBOL[[gse]]
  na_ct <- sum(is.na(mat))
  
  if (na_ct == length(mat)) {
    cat(sprintf("  %-12s  ALL NA — skipping (needs rebuild)\n", gse))
    next
  }
  
  med <- median(mat, na.rm = TRUE)
  mx  <- max(mat, na.rm = TRUE)
  
  if (mx > 1000 || med > 100) {
    cat(sprintf("  %-12s  median=%.1f max=%.1f → APPLYING log2(x+1)\n",
                gse, med, mx))
    expr_store_SYMBOL[[gse]] <- log2(mat + 1)
    new_med <- median(expr_store_SYMBOL[[gse]], na.rm = TRUE)
    new_mx  <- max(expr_store_SYMBOL[[gse]], na.rm = TRUE)
    cat(sprintf("    → After: median=%.2f max=%.2f\n", new_med, new_mx))
  } else {
    cat(sprintf("  %-12s  median=%.2f max=%.2f → already log-scale ✓\n",
                gse, med, mx))
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# 8. SAVE
# ─────────────────────────────────────────────────────────────────────────────
out_path <- file.path(new_rds, "expr_store_SYMBOL_all8.rds")
saveRDS(expr_store_SYMBOL, out_path)
cat("\n✅ Saved:", out_path, "\n")

# ─────────────────────────────────────────────────────────────────────────────
# 9. FINAL QC
# ─────────────────────────────────────────────────────────────────────────────
cat("\n========== EXPR STORE SUMMARY ==========\n")
for (gse in names(expr_store_SYMBOL)) {
  mat     <- expr_store_SYMBOL[[gse]]
  ph_sub  <- pheno %>% filter(GSE == gse)
  matched <- sum(colnames(mat) %in% ph_sub$GSM)
  cat(sprintf("  %-12s  genes=%-6d  samples=%-4d  pheno_match=%d/%d\n",
              gse, nrow(mat), ncol(mat), matched, nrow(ph_sub)))
}

cat("\n========== VALUE RANGE CHECK ==========\n")
for (gse in names(expr_store_SYMBOL)) {
  mat   <- expr_store_SYMBOL[[gse]]
  rng   <- range(mat, na.rm = TRUE)
  na_ct <- sum(is.na(mat))
  cat(sprintf("  %-12s  min=%.2f  median=%.2f  max=%.2f  NA=%d\n",
              gse, rng[1], median(mat, na.rm = TRUE), rng[2], na_ct))
}

cat("\n========== GENE NAME SAMPLES ==========\n")
for (gse in names(expr_store_SYMBOL)) {
  rn <- rownames(expr_store_SYMBOL[[gse]])
  cat(sprintf("  %-12s  n=%d  first5: %s\n",
              gse, length(rn), paste(head(rn, 5), collapse = ", ")))
}

cat("\n========== GENE OVERLAP ==========\n")
gene_lists <- lapply(expr_store_SYMBOL, rownames)
common_genes <- Reduce(intersect, gene_lists)
cat("  Common genes across all", length(expr_store_SYMBOL),
    "cohorts:", length(common_genes), "\n")

gse_names <- names(gene_lists)
cat("\n  Pairwise gene overlap:\n")
for (i in seq_along(gse_names)) {
  for (j in seq_along(gse_names)) {
    if (j > i) {
      ov <- length(intersect(gene_lists[[i]], gene_lists[[j]]))
      cat(sprintf("    %s ∩ %s = %d\n", gse_names[i], gse_names[j], ov))
    }
  }
}

cat("\n✅ Script 01 complete. Proceed to Script 02 (DGE).\n")
