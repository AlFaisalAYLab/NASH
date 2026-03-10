# =============================================================================
# Script 00 – Project Setup + Full Phenotype Harmonization (All 8 Cohorts)
# =============================================================================
# Cohorts: GSE48452, GSE49541, GSE89632, GSE126848, GSE130970,
#          GSE135251, GSE162694, GSE167523
# =============================================================================

library(dplyr)
library(readr)
library(stringr)
library(purrr)

# ─────────────────────────────────────────────────────────────────────────────
# 0. PATHS
# ─────────────────────────────────────────────────────────────────────────────
proj_root  <- "D:/NASH_metaanalysis/"
new_cache  <- file.path(proj_root, "Data", "GEO_cache")
new_rds    <- file.path(proj_root, "Data", "RDS")
new_idmaps <- file.path(proj_root, "Data", "ID_maps")
new_qc     <- file.path(proj_root, "Data", "QC")
new_figs   <- file.path(proj_root, "Figures")
new_tabs   <- file.path(proj_root, "Tables")
new_logs   <- file.path(proj_root, "Logs")
old_cache  <- "D:/NASH_metaanalysis/GEO_cache"

# ─────────────────────────────────────────────────────────────────────────────
# 1. CREATE FOLDER STRUCTURE
# ─────────────────────────────────────────────────────────────────────────────
dirs_needed <- c(new_cache, new_rds, new_idmaps, new_qc,
                 new_figs, new_tabs, new_logs,
                 file.path(new_figs, "Meta"),
                 file.path(new_figs, "DGE"),
                 file.path(new_figs, "WGCNA"))

for (d in dirs_needed) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
    cat("Created:", d, "\n")
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# 2. COPY / SYMLINK CACHE FROM OLD LOCATION
# ─────────────────────────────────────────────────────────────────────────────
old_files <- list.files(old_cache, full.names = TRUE)
for (f in old_files) {
  dest <- file.path(new_cache, basename(f))
  if (!file.exists(dest)) file.copy(f, dest)
}
cat("Cache files available:", length(list.files(new_cache)), "\n")

# ─────────────────────────────────────────────────────────────────────────────
# 3. HELPERS
# ─────────────────────────────────────────────────────────────────────────────

# Read series matrix lines (handles both .txt and .txt.gz)
read_sm_lines <- function(gse_id, cache_dir = new_cache) {
  f <- list.files(cache_dir,
                  pattern = paste0(gse_id, ".*series_matrix.*\\.txt"),
                  full.names = TRUE)
  f_use <- f[!grepl("\\.gz$", f)]
  if (length(f_use) == 0) f_use <- f[1]
  if (length(f_use) == 0 || is.na(f_use)) {
    warning("No series matrix found for ", gse_id); return(character(0))
  }
  read_lines(f_use)
}

# Extract one characteristics field; pattern matched INSIDE quoted values
extract_char_field <- function(lines, pattern) {
  cl  <- lines[grepl("^!Sample_characteristics_ch1", lines)]
  hit <- cl[grepl(pattern, cl, ignore.case = TRUE)][1]
  if (is.na(hit)) return(NULL)
  parts <- str_extract_all(hit, '"([^"]+)"')[[1]]
  str_remove_all(parts, '"')
}

# Extract !Sample_geo_accession → GSM IDs
extract_gsm <- function(lines) {
  l <- lines[grepl("^!Sample_geo_accession", lines)]
  if (length(l) == 0) return(character(0))
  parts <- str_extract_all(l[1], '"([^"]+)"')[[1]]
  str_remove_all(parts, '"')
}

# Extract !Sample_title
extract_titles <- function(lines) {
  l <- lines[grepl("^!Sample_title", lines)]
  if (length(l) == 0) return(character(0))
  parts <- str_extract_all(l[1], '"([^"]+)"')[[1]]
  str_remove_all(parts, '"')
}

# ─────────────────────────────────────────────────────────────────────────────
# 4. BUILD BASE PHENO (GSM + sample_title + GSE)
# ─────────────────────────────────────────────────────────────────────────────
gse_ids <- c("GSE48452","GSE49541","GSE89632","GSE126848",
             "GSE130970","GSE135251","GSE162694","GSE167523")

pheno_list <- lapply(gse_ids, function(gse) {
  lines  <- read_sm_lines(gse)
  gsms   <- extract_gsm(lines)
  titles <- extract_titles(lines)
  n      <- max(length(gsms), length(titles))
  if (n == 0) { warning("Empty for ", gse); return(NULL) }
  tibble(
    GSE          = gse,
    GSM          = if (length(gsms)   == n) gsms   else rep(NA_character_, n),
    sample_title = if (length(titles) == n) titles else rep(NA_character_, n)
  )
})

pheno <- bind_rows(pheno_list) %>%
  mutate(
    disease_label        = NA_character_,
    disease_harmonized   = NA_character_,
    fibrosis_stage       = NA_character_,
    fibrosis_group       = NA_character_,
    nas_score            = NA_integer_,
    nas_group            = NA_character_,
    ballooning           = NA_integer_,
    ballooning_group     = NA_character_,
    lobular_inflammation = NA_integer_,
    inflammation_group   = NA_character_,
    steatosis            = NA_integer_,
    sex                  = NA_character_,
    age                  = NA_real_,
    bmi                  = NA_real_
  )

cat("Base pheno built:", nrow(pheno), "rows\n")

# ─────────────────────────────────────────────────────────────────────────────
# 5. GSE49541  –  fibrosis stage from "Stage:" field
# ─────────────────────────────────────────────────────────────────────────────
cat("\n── GSE49541 ──\n")
lines49 <- read_sm_lines("GSE49541")
idx49   <- which(pheno$GSE == "GSE49541")

stage49 <- extract_char_field(lines49, "Stage:")
if (!is.null(stage49)) {
  stage_clean <- str_remove(stage49, "^Stage:\\s*")
  pheno$fibrosis_stage[idx49] <- case_when(
    str_detect(stage_clean, "advanced") ~ "F3",
    str_detect(stage_clean, "mild")     ~ "F1",
    TRUE                               ~ NA_character_
  )
  pheno$disease_label[idx49] <- "NAFLD"
  cat("  Fibrosis:\n"); print(table(pheno$fibrosis_stage[idx49], useNA = "ifany"))
}

# ─────────────────────────────────────────────────────────────────────────────
# 6. GSE48452  –  group / fibrosis / NAS / inflammation / fat / sex / age
# ─────────────────────────────────────────────────────────────────────────────
cat("\n── GSE48452 ──\n")
lines48 <- read_sm_lines("GSE48452")
idx48   <- which(pheno$GSE == "GSE48452")

grp48 <- extract_char_field(lines48, "group:")
if (!is.null(grp48)) {
  pheno$disease_label[idx48] <- str_remove(grp48, "^group:\\s*")
  cat("  Disease:\n"); print(table(pheno$disease_label[idx48], useNA = "ifany"))
}

fib48 <- extract_char_field(lines48, "fibrosis:")
if (!is.null(fib48)) {
  pheno$fibrosis_stage[idx48] <- paste0("F", str_remove(fib48, "^fibrosis:\\s*"))
  cat("  Fibrosis:\n"); print(table(pheno$fibrosis_stage[idx48], useNA = "ifany"))
}

nas48 <- extract_char_field(lines48, "nas:")
if (!is.null(nas48))
  pheno$nas_score[idx48] <- as.integer(str_extract(nas48, "[0-9]+$"))

infl48 <- extract_char_field(lines48, "inflammation:")
if (!is.null(infl48))
  pheno$lobular_inflammation[idx48] <- as.integer(str_extract(infl48, "[0-9]+$"))

fat48 <- extract_char_field(lines48, "^fat:")
if (!is.null(fat48))
  pheno$steatosis[idx48] <- as.integer(str_extract(fat48, "[0-9]+$"))

sex48 <- extract_char_field(lines48, "Sex:")
if (!is.null(sex48))
  pheno$sex[idx48] <- toupper(substr(str_remove(sex48, "^Sex:\\s*"), 1, 1))

age48 <- extract_char_field(lines48, "age:")
if (!is.null(age48))
  pheno$age[idx48] <- as.integer(str_extract(age48, "[0-9]+$"))

bmi48 <- extract_char_field(lines48, "bmi:")
if (!is.null(bmi48))
  pheno$bmi[idx48] <- as.numeric(str_extract(bmi48, "[0-9.]+$"))

# ─────────────────────────────────────────────────────────────────────────────
# 7. GSE89632  –  diagnosis / fibrosis / ballooning / inflammation / NAS /
#                 steatosis / gender / age
# ─────────────────────────────────────────────────────────────────────────────
cat("\n── GSE89632 ──\n")
lines89 <- read_sm_lines("GSE89632")
idx89   <- which(pheno$GSE == "GSE89632")

diag89 <- extract_char_field(lines89, "diagnosis:")
if (!is.null(diag89)) {
  pheno$disease_label[idx89] <- str_remove(diag89, "^diagnosis:\\s*")
  cat("  Diagnosis:\n"); print(table(pheno$disease_label[idx89]))
}

fib89 <- extract_char_field(lines89, "fibrosis \\(stage\\)")
if (!is.null(fib89)) {
  pheno$fibrosis_stage[idx89] <- paste0("F", str_extract(fib89, "[0-9]+$"))
  cat("  Fibrosis:\n"); print(table(pheno$fibrosis_stage[idx89]))
}

lob89 <- extract_char_field(lines89, "lobular inflammation \\(severity\\)")
if (!is.null(lob89))
  pheno$lobular_inflammation[idx89] <- as.integer(str_extract(lob89, "[0-9]+$"))

bal89 <- extract_char_field(lines89, "ballooning \\(intensity\\)")
if (!is.null(bal89))
  pheno$ballooning[idx89] <- as.integer(str_extract(bal89, "[0-9]+$"))

nas89 <- extract_char_field(lines89, "nafld activity score")
if (!is.null(nas89))
  pheno$nas_score[idx89] <- as.integer(str_extract(nas89, "[0-9]+$"))

ste89 <- extract_char_field(lines89, "steatosis \\(%\\)")
if (!is.null(ste89)) {
  ste_pct <- as.numeric(str_extract(ste89, "[0-9.]+$"))
  pheno$steatosis[idx89] <- case_when(
    ste_pct <  5  ~ 0L,
    ste_pct <= 33 ~ 1L,
    ste_pct <= 66 ~ 2L,
    ste_pct >  66 ~ 3L,
    TRUE          ~ NA_integer_
  )
}

gen89 <- extract_char_field(lines89, "gender:")
if (!is.null(gen89))
  pheno$sex[idx89] <- toupper(substr(str_remove(gen89, "^gender:\\s*"), 1, 1))

age89 <- extract_char_field(lines89, "age \\(y\\)")
if (!is.null(age89))
  pheno$age[idx89] <- as.integer(str_extract(age89, "[0-9]+$"))

bmi89 <- extract_char_field(lines89, "body mass index")
if (!is.null(bmi89))
  pheno$bmi[idx89] <- as.numeric(str_extract(bmi89, "[0-9.]+$"))

# ─────────────────────────────────────────────────────────────────────────────
# 8. GSE126848  –  disease from sample_title; gender from series matrix
# ─────────────────────────────────────────────────────────────────────────────
cat("\n── GSE126848 ──\n")
idx126 <- which(pheno$GSE == "GSE126848")

pheno$disease_label[idx126] <- case_when(
  str_detect(pheno$sample_title[idx126], regex("^NAFL_",    ignore_case = TRUE)) ~ "NAFL",
  str_detect(pheno$sample_title[idx126], regex("^NASH_",    ignore_case = TRUE)) ~ "NASH",
  str_detect(pheno$sample_title[idx126], regex("^Healthy_", ignore_case = TRUE)) ~ "Control",
  TRUE ~ "NAFLD"
)
cat("  Disease:\n"); print(table(pheno$disease_label[idx126]))

gen126 <- extract_char_field(read_sm_lines("GSE126848"), "gender:")
if (!is.null(gen126))
  pheno$sex[idx126] <- toupper(substr(str_remove(gen126, "^gender:\\s*"), 1, 1))

# ─────────────────────────────────────────────────────────────────────────────
# 9. GSE167523  –  disease subtype / gender / age
# ─────────────────────────────────────────────────────────────────────────────
cat("\n── GSE167523 ──\n")
lines167 <- read_sm_lines("GSE167523")
idx167   <- which(pheno$GSE == "GSE167523")

sub167 <- extract_char_field(lines167, "disease subtype:")
if (!is.null(sub167)) {
  pheno$disease_label[idx167] <- str_remove(sub167, "^disease subtype:\\s*")
  cat("  Disease:\n"); print(table(pheno$disease_label[idx167]))
}

gen167 <- extract_char_field(lines167, "gender:")
if (!is.null(gen167))
  pheno$sex[idx167] <- toupper(substr(str_remove(gen167, "^gender:\\s*"), 1, 1))

age167 <- extract_char_field(lines167, "age:")
if (!is.null(age167))
  pheno$age[idx167] <- as.integer(str_extract(age167, "[0-9]+$"))

# ─────────────────────────────────────────────────────────────────────────────
# 10. GSE162694  –  fibrosis stage text; disease from normal vs NAFLD
# ─────────────────────────────────────────────────────────────────────────────
cat("\n── GSE162694 ──\n")
lines162 <- read_sm_lines("GSE162694")
idx162   <- which(pheno$GSE == "GSE162694")

fib162_raw <- extract_char_field(lines162, "fibrosis stage:")
if (!is.null(fib162_raw)) {
  fib162_text <- str_remove(fib162_raw, "^fibrosis stage:\\s*")
  
  # fibrosis_stage: "normal liver histology" → F0; "stage 1" → F1, etc.
  pheno$fibrosis_stage[idx162] <- case_when(
    str_detect(fib162_text, regex("normal", ignore_case = TRUE)) ~ "F0",
    str_detect(fib162_text, "1$")                               ~ "F1",
    str_detect(fib162_text, "2$")                               ~ "F2",
    str_detect(fib162_text, "3$")                               ~ "F3",
    str_detect(fib162_text, "4$")                               ~ "F4",
    TRUE                                                        ~ NA_character_
  )
  
  pheno$disease_label[idx162] <- case_when(
    str_detect(fib162_text, regex("normal", ignore_case = TRUE)) ~ "Control",
    TRUE                                                         ~ "NAFLD"
  )
  cat("  Disease:\n"); print(table(pheno$disease_label[idx162]))
  cat("  Fibrosis:\n"); print(table(pheno$fibrosis_stage[idx162], useNA = "ifany"))
}

nas162 <- extract_char_field(lines162, "nafld activity score|nas score")
if (!is.null(nas162))
  pheno$nas_score[idx162] <- as.integer(str_extract(nas162, "[0-9]+$"))

# NAS from steatosis+inflammation+ballooning if not found directly
bal162 <- extract_char_field(lines162, "ballooning")
if (!is.null(bal162))
  pheno$ballooning[idx162] <- as.integer(str_extract(bal162, "[0-9]+$"))

lob162 <- extract_char_field(lines162, "inflammation")
if (!is.null(lob162))
  pheno$lobular_inflammation[idx162] <- as.integer(str_extract(lob162, "[0-9]+$"))

age162 <- extract_char_field(lines162, "age:")
if (!is.null(age162))
  pheno$age[idx162] <- as.integer(str_extract(age162, "[0-9]+$"))

sex162 <- extract_char_field(lines162, "gender:|sex:")
if (!is.null(sex162))
  pheno$sex[idx162] <- toupper(substr(str_remove(sex162, "^[^:]+:\\s*"), 1, 1))

# ─────────────────────────────────────────────────────────────────────────────
# 11. GSE130970  –  NAS / fibrosis from series matrix; disease from scores
# ─────────────────────────────────────────────────────────────────────────────
cat("\n── GSE130970 ──\n")
lines130 <- read_sm_lines("GSE130970")
idx130   <- which(pheno$GSE == "GSE130970")

fib130 <- extract_char_field(lines130, "fibrosis")
if (!is.null(fib130)) {
  fib130_text <- str_remove(fib130, "^[^:]+:\\s*")
  pheno$fibrosis_stage[idx130] <- paste0("F", str_extract(fib130_text, "[0-9.]+$"))
}

nas130 <- extract_char_field(lines130, "nas|nafld activity")
if (!is.null(nas130))
  pheno$nas_score[idx130] <- as.integer(str_extract(nas130, "[0-9]+$"))

bal130 <- extract_char_field(lines130, "ballooning")
if (!is.null(bal130))
  pheno$ballooning[idx130] <- as.integer(str_extract(bal130, "[0-9]+$"))

lob130 <- extract_char_field(lines130, "inflammation")
if (!is.null(lob130))
  pheno$lobular_inflammation[idx130] <- as.integer(str_extract(lob130, "[0-9]+$"))

age130 <- extract_char_field(lines130, "age:")
if (!is.null(age130))
  pheno$age[idx130] <- as.integer(str_extract(age130, "[0-9]+$"))

sex130 <- extract_char_field(lines130, "gender:|sex:")
if (!is.null(sex130))
  pheno$sex[idx130] <- toupper(substr(str_remove(sex130, "^[^:]+:\\s*"), 1, 1))

# Derive disease_label from NAS + fibrosis
pheno$disease_label[idx130] <- case_when(
  pheno$nas_score[idx130] == 0 & pheno$fibrosis_stage[idx130] == "F0" ~ "Control",
  pheno$nas_score[idx130] >= 5                                         ~ "NASH",
  pheno$nas_score[idx130] >= 1                                         ~ "NAFL",
  TRUE                                                                 ~ NA_character_
)
cat("  Disease:\n"); print(table(pheno$disease_label[idx130], useNA = "ifany"))

# ─────────────────────────────────────────────────────────────────────────────
# 12. GSE135251  –  NAS / fibrosis from series matrix; disease from labels
# ─────────────────────────────────────────────────────────────────────────────
cat("\n── GSE135251 ──\n")
lines135 <- read_sm_lines("GSE135251")
idx135   <- which(pheno$GSE == "GSE135251")

fib135 <- extract_char_field(lines135, "fibrosis")
if (!is.null(fib135)) {
  fib135_text <- str_remove(fib135, "^[^:]+:\\s*")
  pheno$fibrosis_stage[idx135] <- paste0("F", str_extract(fib135_text, "[0-9.]+$"))
}

nas135 <- extract_char_field(lines135, "nas|nafld activity")
if (!is.null(nas135))
  pheno$nas_score[idx135] <- as.integer(str_extract(nas135, "[0-9]+$"))

bal135 <- extract_char_field(lines135, "ballooning")
if (!is.null(bal135))
  pheno$ballooning[idx135] <- as.integer(str_extract(bal135, "[0-9]+$"))

lob135 <- extract_char_field(lines135, "inflammation")
if (!is.null(lob135))
  pheno$lobular_inflammation[idx135] <- as.integer(str_extract(lob135, "[0-9]+$"))

age135 <- extract_char_field(lines135, "age:")
if (!is.null(age135))
  pheno$age[idx135] <- as.integer(str_extract(age135, "[0-9]+$"))

sex135 <- extract_char_field(lines135, "gender:|sex:")
if (!is.null(sex135))
  pheno$sex[idx135] <- toupper(substr(str_remove(sex135, "^[^:]+:\\s*"), 1, 1))

# Derive disease_label
diag135 <- extract_char_field(lines135, "diagnosis|disease|group")
if (!is.null(diag135)) {
  pheno$disease_label[idx135] <- str_remove(diag135, "^[^:]+:\\s*")
} else {
  pheno$disease_label[idx135] <- case_when(
    pheno$nas_score[idx135] == 0 & pheno$fibrosis_stage[idx135] == "F0" ~ "Control",
    pheno$nas_score[idx135] >= 5                                         ~ "NASH",
    pheno$nas_score[idx135] >= 1                                         ~ "NAFL",
    TRUE                                                                 ~ NA_character_
  )
}
cat("  Disease:\n"); print(table(pheno$disease_label[idx135], useNA = "ifany"))

# ─────────────────────────────────────────────────────────────────────────────
# 13. REBUILD ALL HARMONIZED GROUP COLUMNS
# ─────────────────────────────────────────────────────────────────────────────
pheno <- pheno %>%
  mutate(
    # disease_harmonized:
    #   – HC / control / healthy / normal → Control
    #   – NAFL, SS, Steatosis → NAFL
    #   – NASH* (any suffix) → NASH
    #   – NAFLD → NAFLD (umbrella, no histology)
    disease_harmonized = case_when(
      str_detect(disease_label,
                 regex("^control$|^HC$|healthy|normal", ignore_case = TRUE)) ~ "Control",
      str_detect(disease_label,
                 regex("^NAFL$|^SS$|^Steatosis$|simple steatosis", ignore_case = TRUE)) ~ "NAFL",
      str_detect(disease_label,
                 regex("^NASH", ignore_case = TRUE)) ~ "NASH",
      str_detect(disease_label,
                 regex("NAFLD", ignore_case = TRUE)) ~ "NAFLD",
      TRUE ~ NA_character_
    ),
    
    # fibrosis_group:  F0–F2 + half-stages → Mild;  F3–F4 + half → Advanced
    fibrosis_group = case_when(
      fibrosis_stage %in% c("F0", "F1", "F2") ~ "Mild",
      fibrosis_stage == "F0.5"                 ~ "Mild",
      fibrosis_stage %in% c("F3", "F4")        ~ "Advanced",
      fibrosis_stage == "F3.5"                 ~ "Advanced",
      TRUE                                     ~ NA_character_    # FNA → NA
    ),
    
    nas_group = case_when(
      nas_score <= 4 ~ "Low_NAS",
      nas_score >= 5 ~ "High_NAS",
      TRUE           ~ NA_character_
    ),
    
    ballooning_group = case_when(
      ballooning == 0          ~ "No_Ballooning",
      ballooning %in% c(1, 2) ~ "Ballooning",
      TRUE                    ~ NA_character_
    ),
    
    inflammation_group = case_when(
      lobular_inflammation == 0       ~ "No_Inflammation",
      lobular_inflammation %in% 1:3  ~ "Inflammation",
      TRUE                           ~ NA_character_
    )
  )

# ─────────────────────────────────────────────────────────────────────────────
# 14. SAVE
# ─────────────────────────────────────────────────────────────────────────────
saveRDS(pheno, file.path(new_rds, "pheno_harmonized.rds"))
write_csv(pheno, file.path(new_idmaps, "pheno_harmonized.csv"))
cat("\n✅ pheno_harmonized.rds saved:", nrow(pheno), "rows\n")

# ─────────────────────────────────────────────────────────────────────────────
# 15. FINAL SUMMARY TABLES
# ─────────────────────────────────────────────────────────────────────────────
cat("\n========== HARMONIZED GROUP COUNTS ==========\n")
pheno %>%
  group_by(GSE) %>%
  summarise(
    n        = n(),
    Mild     = sum(fibrosis_group == "Mild",     na.rm = TRUE),
    Advanced = sum(fibrosis_group == "Advanced", na.rm = TRUE),
    Low_NAS  = sum(nas_group == "Low_NAS",       na.rm = TRUE),
    High_NAS = sum(nas_group == "High_NAS",      na.rm = TRUE),
    Control  = sum(disease_harmonized == "Control", na.rm = TRUE),
    NAFL     = sum(disease_harmonized == "NAFL",    na.rm = TRUE),
    NASH     = sum(disease_harmonized == "NASH",    na.rm = TRUE),
    NAFLD    = sum(disease_harmonized == "NAFLD",   na.rm = TRUE),
    NA_dis   = sum(is.na(disease_harmonized)),
    .groups  = "drop"
  ) %>%
  print(n = Inf)

cat("\n========== CONTRASTS AVAILABLE ==========\n")
pheno %>%
  group_by(GSE) %>%
  summarise(
    n               = n(),
    can_fibrosis    = sum(!is.na(fibrosis_group)) >= 6,
    can_NAS         = sum(!is.na(nas_group))      >= 6,
    can_NASH_ctrl   = sum(disease_harmonized == "NASH",    na.rm = TRUE) >= 3 &
      sum(disease_harmonized == "Control", na.rm = TRUE) >= 3,
    can_NASH_NAFL   = sum(disease_harmonized == "NASH",    na.rm = TRUE) >= 3 &
      sum(disease_harmonized == "NAFL",    na.rm = TRUE) >= 3,
    can_ballooning  = sum(!is.na(ballooning_group))    >= 6,
    can_inflamm     = sum(!is.na(inflammation_group))  >= 6,
    .groups         = "drop"
  ) %>%
  print(n = Inf)

cat("\n========== EXPRESSION STORE CHECK ==========\n")
expr_candidates <- c(
  file.path(new_rds, "expr_store_SYMBOL_all8.rds"),
  "D:/NASH_metaanalysis/old_analysis/NASH_Outputs/RDS_Objects/expr_store_all8.rds"
)
for (f in expr_candidates) {
  cat(sprintf("  %s  →  %s\n", f, if (file.exists(f)) "FOUND ✅" else "missing ❌"))
}

eset_files <- list.files(old_cache,
                         pattern = "_eset\\.rds$",
                         full.names = TRUE)
cat("\nIndividual eset files (for rebuild if needed):\n")
cat(basename(eset_files), sep = "\n")
cat("\n✅ Script 00 complete. Proceed to Script 01.\n")
