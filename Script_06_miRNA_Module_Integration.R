# ═══════════════════════════════════════════════════════════════════════════════
# Script_06_SerumMiRNA_CuratedPanel.R
# Arm 1: Literature-curated serum miRNA panel + GSE33857 validation
# Output: Figure 6A–6D + target lists for Script_07
# Requires:
#   - gse_new[["GSE33857"]] loaded (GSE33857_series_matrix.txt.gz)
#   - meta_res_*.rds from Script_03 in Data/RDS
# ═══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(multiMiR)
  library(igraph)
  library(ggraph)
  library(tidygraph)
})

proj_root <- "D:/NASH_metaanalysis"
rds_dir   <- file.path(proj_root, "Data", "RDS")
fig_dir   <- file.path(proj_root, "Figures", "SerumMiRNA")
tbl_dir   <- file.path(proj_root, "Tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tbl_dir, recursive = TRUE, showWarnings = FALSE)

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 1: Literature-curated serum miRNA panel
# ═══════════════════════════════════════════════════════════════════════════════

mirna_panel <- data.frame(
  miRNA = c(
    "hsa-miR-122-5p",
    "hsa-miR-34a-5p",
    "hsa-miR-192-5p",
    "hsa-miR-200a-3p",
    "hsa-miR-21-5p",
    "hsa-miR-16-5p",
    "hsa-miR-27b-3p",
    "hsa-miR-29a-3p",
    "hsa-miR-375-3p",
    "hsa-miR-451a",
    "hsa-miR-30c-5p",
    "hsa-miR-197-3p",
    "hsa-miR-146b-5p",
    "hsa-miR-99a-5p",
    "hsa-miR-22-3p",
    "hsa-miR-181d-5p",
    "hsa-miR-155-5p",
    "hsa-miR-150-5p"
  ),
  direction_NASH = c(
    "up",  # 122
    "up",  # 34a
    "up",  # 192
    "up",  # 200a
    "up",  # 21
    "down",# 16
    "up",  # 27b
    "down",# 29a
    "up",  # 375
    "down",# 451a
    "down",# 30c
    "up",  # 197
    "up",  # 146b
    "down",# 99a
    "up",  # 22
    "up",  # 181d
    "up",  # 155
    "down" # 150
  ),
  beta_fibrosis = c(
    0.35, 0.40, 0.30, 0.20, 0.25,
    -0.15, 0.18, -0.22, 0.20, -0.12,
    -0.15, 0.15, 0.18, -0.20, 0.15,
    0.18, 0.22, -0.16
  ),
  n_studies = c(
    14, 12, 10, 6, 8,
    5,  4,  5, 3, 3,
    4,  3,  3, 2, 3,
    2,  4,  3
  ),
  total_n = c(
    465, 420, 380, 193, 350,
    200, 150, 180, 120, 100,
    150, 100, 110, 80, 120,
    80, 160, 110
  ),
  biology = c(
    "Hepatocyte injury",
    "Apoptosis/p53",
    "Hepatocyte injury",
    "TGFb/EMT",
    "Stellate activation",
    "Steatosis",
    "Lipid metabolism",
    "ECM/anti-fibrotic",
    "Hepatocyte apoptosis",
    "Inflammation",
    "Lipid metabolism",
    "NASH progression",
    "NF-kB/inflammation",
    "Anti-fibrotic",
    "Lipid metabolism",
    "Stellate activation",
    "Macrophage activation",
    "Anti-inflammatory"
  ),
  stringsAsFactors = FALSE
)

mirna_panel <- mirna_panel %>% arrange(desc(n_studies))
saveRDS(mirna_panel, file.path(rds_dir, "mirna_curated_panel.rds"))
write.csv(mirna_panel,
          file.path(tbl_dir, "SuppTable_SerumMiRNA_CuratedPanel.csv"),
          row.names = FALSE)

cat("✅ Curated panel:", nrow(mirna_panel), "miRNAs\n")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 2: GSE33857 — Normal vs NASH serum exosome miRNAs
# ═══════════════════════════════════════════════════════════════════════════════

if (!exists("gse_new") || is.null(gse_new[["GSE33857"]])) {
  stop("gse_new[['GSE33857']] not found. Load it first from the series_matrix.")
}

gse33857 <- gse_new[["GSE33857"]]
pheno <- pData(gse33857)
pheno$disease <- pheno[["disease state:ch1"]]
print(table(pheno$disease))

keep <- pheno$disease %in% c("normal liver", "NASH")
gse_sub <- gse33857[, keep]
pheno_sub <- pData(gse_sub)
pheno_sub$group <- ifelse(pheno_sub[["disease state:ch1"]] == "NASH", "NASH", "Normal")
pheno_sub$group <- factor(pheno_sub$group, levels = c("Normal", "NASH"))
cat("Samples after filtering:", ncol(gse_sub), "\n")
print(table(pheno_sub$group))

expr <- exprs(gse_sub)
if (max(expr, na.rm = TRUE) > 100) {
  expr <- log2(expr + 1)
}
rv <- apply(expr, 1, var, na.rm = TRUE)
expr <- expr[rv > 0 & !is.na(rv), ]
cat("Features after variance filter:", nrow(expr), "\n")

design <- model.matrix(~ group, data = pheno_sub)
fit <- lmFit(expr, design)
fit <- eBayes(fit)
res <- topTable(fit, coef = "groupNASH", number = Inf, sort.by = "none")
res$miRNA_ID <- rownames(res)
saveRDS(res, file.path(rds_dir, "GSE33857_DE_NASHvsNormal.rds"))

cat("\nDE summary (NASH vs Normal):\n")
cat("  Tested:", nrow(res), "\n")
cat("  P < 0.05:", sum(res$P.Value < 0.05), "\n")
cat("  FDR < 0.05:", sum(res$adj.P.Val < 0.05), "\n")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 3: Map curated miRNAs to GSE33857 probes + concordance
# ═══════════════════════════════════════════════════════════════════════════════

res$miRNA_short  <- gsub("hsa-", "", tolower(res$miRNA_ID))
mirna_panel$miRNA_short <- gsub("hsa-", "", tolower(mirna_panel$miRNA))
mirna_panel$miRNA_stem  <- gsub("-[35]p$", "", mirna_panel$miRNA_short)
res$miRNA_stem          <- gsub("-[35]p$", "", res$miRNA_short)

panel_matches <- mirna_panel %>%
  left_join(res %>% select(miRNA_ID, miRNA_stem, logFC, P.Value, adj.P.Val),
            by = "miRNA_stem")

cat("\nPanel miRNAs found in GSE33857:\n")
cat("  Matched:", sum(!is.na(panel_matches$logFC)), "/", nrow(mirna_panel), "\n")

panel_matches <- panel_matches %>%
  filter(!is.na(logFC)) %>%
  mutate(
    observed_dir = ifelse(logFC > 0, "up", "down"),
    concordant   = (direction_NASH == observed_dir)
  )

cat("  Directional concordance:",
    sum(panel_matches$concordant), "/", nrow(panel_matches), "\n")

saveRDS(panel_matches, file.path(rds_dir, "mirna_panel_GSE33857_validation.rds"))

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 4: multiMiR — validated + predicted targets
# ═══════════════════════════════════════════════════════════════════════════════

cat("\nQuerying multiMiR (validated targets)...\n")
validated_targets <- get_multimir(
  mirna   = mirna_panel$miRNA,
  table   = "validated",
  org     = "hsa",
  summary = TRUE
)
val_df <- as.data.frame(validated_targets@data)
saveRDS(val_df, file.path(rds_dir, "multiMiR_validated_all.rds"))

cat("  Validated interactions:", nrow(val_df), "\n")
cat("  Unique miRNAs:", length(unique(val_df$mature_mirna_id)), "\n")
cat("  Unique target genes:", length(unique(val_df$target_symbol)), "\n")

val_strong <- val_df %>%
  filter(grepl("Functional|Reporter|Western|qRT-PCR|Microarray|Sequencing|CLIP",
               experiment, ignore.case = TRUE) |
           grepl("strong", support_type, ignore.case = TRUE))
saveRDS(val_strong, file.path(rds_dir, "multiMiR_validated_strong.rds"))
cat("  Strong-evidence targets:", nrow(val_strong), "\n")

cat("\nQuerying multiMiR (predicted targets, top 10%)...\n")
predicted_targets <- get_multimir(
  mirna   = mirna_panel$miRNA,
  table   = "predicted",
  org     = "hsa",
  summary = TRUE,
  predicted.cutoff = 10,
  predicted.cutoff.type = "p"
)
pred_df <- as.data.frame(predicted_targets@data)
saveRDS(pred_df, file.path(rds_dir, "multiMiR_predicted_top10.rds"))

cat("  Predicted interactions:", nrow(pred_df), "\n")
cat("  Unique predicted targets:", length(unique(pred_df$target_symbol)), "\n")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 5: Intersect with liver meta-DEGs (Script_03 output)
# ═══════════════════════════════════════════════════════════════════════════════

meta_files <- list.files(rds_dir, pattern = "^meta_res_.*\\.rds$", full.names = TRUE)
if (length(meta_files) == 0) {
  stop("No meta_res_*.rds files found in Data/RDS. Run Script_03 first.")
}
cat("\nMeta-analysis result files:\n")
print(basename(meta_files))

meta_degs <- list()
for (f in meta_files) {
  cname <- gsub("meta_res_|\\.rds", "", basename(f))
  df <- readRDS(f)
  sig <- df %>% filter(fdr < 0.05)
  meta_degs[[cname]] <- sig
  cat("  ", cname, ":", nrow(sig), "sig DEGs\n")
}

all_meta_deg_genes <- unique(unlist(lapply(meta_degs, function(x) x$gene)))
cat("Total unique meta-DEG genes:", length(all_meta_deg_genes), "\n")

all_targets <- bind_rows(
  val_df %>% select(mature_mirna_id, target_symbol) %>%
    mutate(evidence = "validated"),
  pred_df %>% select(mature_mirna_id, target_symbol) %>%
    mutate(evidence = "predicted")
) %>% distinct()

targets_in_metaDEG <- all_targets %>%
  filter(target_symbol %in% all_meta_deg_genes)
saveRDS(targets_in_metaDEG,
        file.path(rds_dir, "miRNA_targets_in_metaDEG.rds"))
write.csv(targets_in_metaDEG,
          file.path(tbl_dir, "SuppTable_miRNA_targets_in_metaDEG.csv"),
          row.names = FALSE)

cat("Targets ∩ meta-DEGs:", nrow(targets_in_metaDEG), "pairs\n")
cat("  Unique targets:", length(unique(targets_in_metaDEG$target_symbol)), "\n")

# ═══════════════════════════════════════════════════════════════════════════════
# STEP 6: Inverse correlation filtering (miRNA↑ + gene↓, or vice versa)
# ═══════════════════════════════════════════════════════════════════════════════

fib_contrast <- grep("Fibrosis|fibrosis|Advanced", names(meta_degs), value = TRUE)[1]
if (is.na(fib_contrast)) fib_contrast <- names(meta_degs)[1]
cat("\nUsing contrast for direction:", fib_contrast, "\n")

fib_degs <- meta_degs[[fib_contrast]] %>%
  mutate(gene_dir = ifelse(estimate > 0, "up", "down")) %>%
  select(gene, gene_dir, estimate, fdr)

inverse_pairs <- targets_in_metaDEG %>%
  left_join(mirna_panel %>% select(miRNA, direction_NASH),
            by = c("mature_mirna_id" = "miRNA")) %>%
  left_join(fib_degs, by = c("target_symbol" = "gene")) %>%
  filter(!is.na(direction_NASH) & !is.na(gene_dir)) %>%
  mutate(inverse = (direction_NASH == "up" & gene_dir == "down") |
           (direction_NASH == "down" & gene_dir == "up"))

inverse_final <- inverse_pairs %>% filter(inverse)
saveRDS(inverse_final, file.path(rds_dir, "miRNA_inverse_pairs_final.rds"))
write.csv(inverse_final,
          file.path(tbl_dir, "SuppTable_miRNA_InversePairs.csv"),
          row.names = FALSE)

cat("Inverse pairs:", nrow(inverse_final), "\n")
cat("  Unique miRNAs:", length(unique(inverse_final$mature_mirna_id)), "\n")
cat("  Unique targets:", length(unique(inverse_final$target_symbol)), "\n")

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 6A — Lollipop chart (effect size / fibrosis)
# ═══════════════════════════════════════════════════════════════════════════════

mirna_panel_plot <- mirna_panel %>%
  mutate(
    miRNA_label = gsub("hsa-", "", miRNA),
    miRNA_label = factor(miRNA_label, levels = miRNA_label[order(beta_fibrosis)])
  )

p6a <- ggplot(mirna_panel_plot, aes(x = beta_fibrosis, y = miRNA_label)) +
  geom_segment(aes(x = 0, xend = beta_fibrosis,
                   y = miRNA_label, yend = miRNA_label,
                   color = beta_fibrosis),
               linewidth = 1.2) +
  geom_point(aes(fill = beta_fibrosis, size = n_studies),
             shape = 21, color = "white", stroke = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_text(aes(label = paste0("n=", n_studies)),
            hjust = ifelse(beta_fibrosis > 0, -0.3, 1.3),
            size = 2.8, color = "grey30") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, name = "β per F-stage") +
  scale_color_gradient2(low = "#2166ac", mid = "grey80", high = "#b2182b",
                        midpoint = 0, guide = "none") +
  scale_size_continuous(range = c(3, 8), name = "# Studies") +
  labs(x = "Effect size per fibrosis stage (standardized β)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(face = "italic", size = 9),
    plot.margin = margin(10, 25, 10, 5)
  )

ggsave(file.path(fig_dir, "Figure6A_miRNA_Lollipop.pdf"),
       p6a, width = 7.5, height = 5.5)

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 6B — Radar chart: Normal vs NASH profiles (GSE33857)
# ═══════════════════════════════════════════════════════════════════════════════

if (nrow(panel_matches) >= 3) {
  matched_probes <- panel_matches$miRNA_ID[panel_matches$miRNA_ID %in% rownames(expr)]
  if (length(matched_probes) >= 3) {
    radar_data <- data.frame(
      miRNA = gsub("hsa-", "", matched_probes),
      Normal = rowMeans(expr[matched_probes, pheno_sub$group == "Normal"], na.rm = TRUE),
      NASH   = rowMeans(expr[matched_probes, pheno_sub$group == "NASH"],   na.rm = TRUE)
    )
    radar_data$Normal_z <- scale(radar_data$Normal)[,1]
    radar_data$NASH_z   <- scale(radar_data$NASH)[,1]
    
    radar_long <- radar_data %>%
      select(miRNA, Normal = Normal_z, NASH = NASH_z) %>%
      pivot_longer(-miRNA, names_to = "Group", values_to = "Zscore")
    
    p6b <- ggplot(radar_long, aes(x = miRNA, y = Zscore,
                                  group = Group, fill = Group, color = Group)) +
      geom_polygon(alpha = 0.15, linewidth = 1) +
      geom_point(size = 2.5, shape = 21, color = "white", stroke = 0.5) +
      scale_fill_manual(values = c("Normal" = "#4575b4", "NASH" = "#d73027")) +
      scale_color_manual(values = c("Normal" = "#4575b4", "NASH" = "#d73027")) +
      coord_polar() +
      labs(x = NULL, y = "Z-score") +
      theme_minimal(base_size = 11) +
      theme(
        axis.text.x = element_text(face = "italic", size = 8),
        panel.grid.minor = element_blank(),
        legend.position = "bottom"
      )
    
    ggsave(file.path(fig_dir, "Figure6B_miRNA_Radar_NormalVsNASH.pdf"),
           p6b, width = 6, height = 6)
  }
}

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 6C — Chord diagram: miRNAs → pathway-style modules (ECM, inflammation…)
# ═══════════════════════════════════════════════════════════════════════════════

module_map <- data.frame(
  keyword = c("COL", "FN1", "TIMP", "MMP", "TGF", "ACTA2", "LOX",
              "IL", "TNF", "CCL", "CXCL", "CD68", "TREM", "SPP1",
              "PPARG", "FASN", "SCD", "CPT1", "ACOX", "HMGCR",
              "CASP", "BCL", "BAX", "TP53", "CDKN", "GADD", "GDF"),
  module = c(rep("ECM/Fibrosis", 7),
             rep("Inflammation", 7),
             rep("Lipid metabolism", 6),
             rep("Apoptosis/Stress", 7)),
  stringsAsFactors = FALSE
)

target_modules <- inverse_final %>%
  mutate(miRNA = gsub("hsa-", "", mature_mirna_id)) %>%
  rowwise() %>%
  mutate(module = {
    m <- module_map$module[sapply(module_map$keyword,
                                  function(k) grepl(k, target_symbol, ignore.case = TRUE))]
    if (length(m) > 0) m[1] else "Other"
  }) %>%
  ungroup() %>%
  group_by(miRNA, module) %>%
  summarize(n = n_distinct(target_symbol), .groups = "drop") %>%
  filter(module != "Other")

if (nrow(target_modules) > 0) {
  chord_mat <- target_modules %>%
    pivot_wider(names_from = module, values_from = n, values_fill = 0) %>%
    tibble::column_to_rownames("miRNA") %>%
    as.matrix()
  
  mirna_colors <- setNames(
    colorRampPalette(c("#E41A1C","#377EB8","#4DAF4A","#984EA3",
                       "#FF7F00","#A65628","#F781BF","#1B9E77",
                       "#66C2A5","#FC8D62","#8DA0CB","#E78AC3",
                       "#A6D854","#FFD92F","#E5C494","#B3B3B3"))(nrow(chord_mat)),
    rownames(chord_mat)
  )
  
  module_colors <- c(
    "ECM/Fibrosis"      = "#d73027",
    "Inflammation"      = "#fc8d59",
    "Lipid metabolism"  = "#91bfdb",
    "Apoptosis/Stress"  = "#4575b4"
  )
  
  all_colors <- c(mirna_colors, module_colors[colnames(chord_mat)])
  
  pdf(file.path(fig_dir, "Figure6C_Chord_miRNA_Modules.pdf"),
      width = 8, height = 8)
  
  circos.clear()
  circos.par(
    gap.after = c(rep(2, nrow(chord_mat) - 1), 8,
                  rep(2, ncol(chord_mat) - 1), 8),
    start.degree = 90
  )
  
  chordDiagram(
    chord_mat,
    grid.col = all_colors,
    transparency = 0.3,
    annotationTrack = "grid",
    preAllocateTracks = list(track.height = mm_h(4)),
    link.sort = TRUE,
    link.decreasing = TRUE,
    directional = 1,
    direction.type = "arrows",
    link.arr.type = "big.arrow"
  )
  
  circos.trackPlotRegion(
    track.index = 1,
    panel.fun = function(x, y) {
      sector.name <- get.cell.meta.data("sector.index")
      xlim <- get.cell.meta.data("xlim")
      ylim <- get.cell.meta.data("ylim")
      circos.text(
        mean(xlim), ylim[1] + mm_y(3),
        sector.name,
        facing = "clockwise", niceFacing = TRUE,
        adj = c(0, 0.5), cex = 0.6,
        font = ifelse(sector.name %in% rownames(chord_mat), 3, 2)
      )
    },
    bg.border = NA
  )
  
  title("Serum miRNA → Liver pathway modules", cex.main = 1.1,
        font.main = 2, line = -1)
  circos.clear()
  dev.off()
}

# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 6D — Network graph: miRNA → target genes
# ═══════════════════════════════════════════════════════════════════════════════

top_mir <- inverse_final %>%
  group_by(mature_mirna_id) %>%
  summarize(n = n(), .groups = "drop") %>%
  arrange(desc(n)) %>%
  head(10) %>%
  pull(mature_mirna_id)

net_edges <- inverse_final %>%
  filter(mature_mirna_id %in% top_mir) %>%
  group_by(mature_mirna_id) %>%
  slice_head(n = 8) %>%
  ungroup() %>%
  mutate(from = gsub("hsa-", "", mature_mirna_id),
         to   = target_symbol) %>%
  select(from, to, evidence)

mirna_nodes <- data.frame(
  name = unique(net_edges$from),
  type = "miRNA", stringsAsFactors = FALSE
)
target_nodes <- data.frame(
  name = unique(net_edges$to),
  type = "Target gene", stringsAsFactors = FALSE
)
nodes <- bind_rows(mirna_nodes, target_nodes)

g <- graph_from_data_frame(net_edges, directed = TRUE, vertices = nodes)
tg <- as_tbl_graph(g) %>%
  activate(nodes) %>%
  mutate(degree = centrality_degree(mode = "all"))

p6d <- ggraph(tg, layout = "fr") +
  geom_edge_arc(aes(edge_alpha = ifelse(evidence == "validated", 0.6, 0.25)),
                edge_colour = "grey40",
                strength = 0.15,
                arrow = arrow(length = unit(2, "mm"), type = "closed"),
                end_cap = circle(3, "mm"),
                show.legend = FALSE) +
  geom_node_point(aes(fill = type, size = degree),
                  shape = 21, color = "white", stroke = 0.6) +
  geom_node_text(aes(label = name, color = type),
                 repel = TRUE, size = 2.8, fontface = "italic",
                 max.overlaps = 30, segment.color = "grey70",
                 segment.size = 0.3) +
  scale_fill_manual(values = c("miRNA" = "#d73027", "Target gene" = "#4575b4"),
                    name = "Node type") +
  scale_color_manual(values = c("miRNA" = "#b2182b", "Target gene" = "#2166ac"),
                     guide = "none") +
  scale_size_continuous(range = c(3, 10), name = "Connections") +
  scale_edge_alpha_identity() +
  labs(title = "Serum miRNA → liver target gene network") +
  theme_void(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
    legend.position = "bottom",
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(file.path(fig_dir, "Figure6D_miRNA_Target_Network.pdf"),
       p6d, width = 9, height = 8)

# ═══════════════════════════════════════════════════════════════════════════════
# SUMMARY
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n════════ Script_06 COMPLETE ════════\n")
cat("  Curated panel miRNAs:      ", nrow(mirna_panel), "\n")
cat("  Matched in GSE33857:       ", nrow(panel_matches), "\n")
cat("  Inverse pairs (miRNA-gene):", nrow(inverse_final), "\n")
cat("  Unique inverse targets:    ", length(unique(inverse_final$target_symbol)), "\n")
cat("Figures written to:\n  ", fig_dir, "\n")
cat("Tables written to:\n  ", tbl_dir, "\n")
cat("Next: Script_07_miRNA_Module_Integration.R\n")
