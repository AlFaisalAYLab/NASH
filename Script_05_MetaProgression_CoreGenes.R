# ═══════════════════════════════════════════════════════════════════════════════
# STEP: Generate Figure 5 panels with fixed colours
# ═══════════════════════════════════════════════════════════════════════════════

library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ggdist)
library(ggpubfigs)
library(ggpubr)
library(dplyr)

# ── 1. Three heatmap colour options (5A) ─────────────────────────────────────

rng <- range(meta_mat_ordered, na.rm = TRUE)
lim <- min(max(abs(rng)), 1.2)

col_opt1 <- colorRamp2(c(-lim, -lim/2, 0, lim/2, lim),
                       c("#2166ac","#92c5de","white","#f4a582","#b2182b"))

col_opt2 <- colorRamp2(c(-lim, -lim/2, 0, lim/2, lim),
                       c("#542788","#b2abd2","#f7f7f7","#fdb863","#e66101"))

col_opt3 <- colorRamp2(c(-lim, 0, lim),
                       c("#1a9876","white","#d73027"))

opt_names <- c("Opt1_Navy_White_Red","Opt2_Purple_White_Orange","Opt3_Teal_White_Coral")
opt_cols  <- list(col_opt1, col_opt2, col_opt3)

for (i in 1:3) {
  ht <- Heatmap(
    meta_mat_ordered, name = "Meta log2FC", col = opt_cols[[i]],
    row_split = row_split, row_gap = unit(2,"mm"),
    row_title_rot = 0, row_title_gp = gpar(fontsize=9, fontface="bold"),
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = FALSE, show_column_names = TRUE,
    column_names_rot = 0, column_names_centered = TRUE,
    column_names_gp = gpar(fontsize=9, fontface="bold"),
    top_annotation = top_ha, left_annotation = left_ha, right_annotation = right_ha,
    border = TRUE, rect_gp = gpar(col="white", lwd=0.3),
    heatmap_legend_param = list(
      title = "Meta log2FC",
      at = round(c(-lim, -lim/2, 0, lim/2, lim), 2),
      legend_height = unit(3,"cm"),
      title_gp = gpar(fontsize=9, fontface="bold"),
      labels_gp = gpar(fontsize=8)
    )
  )
  pdf(file.path(fig_dir, paste0("Figure5A_Heatmap_",opt_names[i],".pdf")), width=7, height=8)
  draw(ht, heatmap_legend_side="right", annotation_legend_side="right",
       padding=unit(c(5,5,5,15),"mm"), merge_legend=TRUE)
  dev.off()
  cat("✅ Figure5A_Heatmap_", opt_names[i], ".pdf\n")
}


# ── 2. Raincloud: Prog score vs Fibrosis (5B) ────────────────────────────────

fib_pal <- c("0"="#fee0b6","1"="#f1a340","2"="#f46d43","3"="#d73027","4"="#a50026")

p_rain_fib <- prog_annot %>%
  filter(!is.na(FibStage)) %>%
  ggplot(aes(x = FibStage, y = ProgScore_PC1, fill = FibStage)) +
  stat_halfeye(adjust=0.5, width=0.6, .width=0, justification=-0.2,
               point_colour=NA, alpha=0.7) +
  geom_boxplot(width=0.12, outlier.shape=NA, alpha=0.5) +
  geom_jitter(width=0.05, alpha=0.25, size=0.6, shape=16, color="grey30") +
  stat_compare_means(
    comparisons = list(c("0","4"), c("1","4"), c("0","2")),
    method="wilcox.test", label="p.signif", tip.length=0.01,
    step.increase=0.06, size=3.5) +
  scale_fill_manual(values = fib_pal) +
  labs(x="Fibrosis stage", y="Progression score (PC1)") +
  theme_classic(base_size=11) +
  theme(legend.position="none", axis.title=element_text(face="bold"),
        plot.margin=margin(10,15,10,10))

ggsave(file.path(fig_dir, "Figure5B_ProgScore_vs_Fibrosis_Raincloud.pdf"),
       p_rain_fib, width=5.5, height=5)
cat("✅ Figure5B\n")


# ── 3. Raincloud: Prog score vs NAS (5C) ─────────────────────────────────────

nas_levels <- sort(unique(as.numeric(as.character(prog_annot$NAS[!is.na(prog_annot$NAS)]))))
nas_pal <- setNames(
  colorRampPalette(c("#c6dbef","#4292c6","#08519c","#08306b"))(length(nas_levels)),
  as.character(nas_levels))

p_rain_nas <- prog_annot %>%
  filter(!is.na(NAS)) %>%
  ggplot(aes(x = NAS, y = ProgScore_PC1, fill = NAS)) +
  stat_halfeye(adjust=0.5, width=0.6, .width=0, justification=-0.2,
               point_colour=NA, alpha=0.7) +
  geom_boxplot(width=0.12, outlier.shape=NA, alpha=0.5) +
  geom_jitter(width=0.05, alpha=0.25, size=0.6, shape=16, color="grey30") +
  stat_compare_means(
    comparisons = list(c(as.character(min(nas_levels)), as.character(max(nas_levels)))),
    method="wilcox.test", label="p.signif", tip.length=0.01, size=3.5) +
  scale_fill_manual(values = nas_pal) +
  labs(x="NAS score", y="Progression score (PC1)") +
  theme_classic(base_size=11) +
  theme(legend.position="none", axis.title=element_text(face="bold"),
        plot.margin=margin(10,15,10,10))

ggsave(file.path(fig_dir, "Figure5C_ProgScore_vs_NAS_Raincloud.pdf"),
       p_rain_nas, width=6, height=5)
cat("✅ Figure5C\n")


# ── 4. Trend line with high-contrast colours (5D) ────────────────────────────

cohort_pal <- c(
  "GSE135251"="#E41A1C", "GSE162694"="#377EB8", "GSE130970"="#4DAF4A",
  "GSE49541" ="#984EA3", "GSE89632" ="#FF7F00", "GSE48452" ="#A65628",
  "GSE126848"="#F781BF", "GSE167523"="#1B9E77"
)

df_trend <- prog_annot %>%
  filter(!is.na(FibStage)) %>%
  mutate(FibStage_num = as.numeric(as.character(FibStage)))

spearman_labs <- df_trend %>%
  group_by(GSE) %>%
  summarize(
    rho = cor(FibStage_num, ProgScore_PC1, method="spearman", use="pair"),
    p   = cor.test(FibStage_num, ProgScore_PC1, method="spearman")$p.value,
    .groups="drop") %>%
  mutate(
    stars = ifelse(p<0.001,"***", ifelse(p<0.01,"**", ifelse(p<0.05,"*","ns"))),
    label = paste0(GSE,"  ρ=", round(rho,2)," ", stars))

cohort_labels <- setNames(spearman_labs$label, spearman_labs$GSE)

p_trend <- ggplot(df_trend,
                  aes(x=FibStage_num, y=ProgScore_PC1, color=GSE, group=GSE)) +
  stat_summary(fun=mean, geom="line", linewidth=1.2, alpha=0.9) +
  stat_summary(fun=mean, geom="point", size=2.5) +
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.15, linewidth=0.6, alpha=0.7) +
  scale_x_continuous(breaks=0:4, labels=paste0("F",0:4)) +
  scale_color_manual(values=cohort_pal, labels=cohort_labels) +
  labs(x="Fibrosis stage", y="Progression score (mean ± SE)", color=NULL) +
  theme_classic(base_size=11) +
  theme(legend.position="right",
        legend.text=element_text(size=7, family="mono"),
        axis.title=element_text(face="bold"))

ggsave(file.path(fig_dir, "Figure5D_ProgScore_Fibrosis_TrendLine_perCohort.pdf"),
       p_trend, width=7.5, height=4.5)
cat("✅ Figure5D\n")

cat("\n🎉 All Figure 5 panels generated! Check the 3 heatmap options and tell me which one.\n")
