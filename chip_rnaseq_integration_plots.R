#!/usr/bin/env Rscript
############################################################
# chip_rnaseq_integration_plots.R
#
# Purpose:
#   Reproducible plotting script for promoter H3Ac (ChIP-seq) and RNA-seq integration:
#     1) Overall promoter H3Ac density (Col-0 vs hda6)
#     2) Heat GO density (GO:0009408)
#     3) ROS/Oxidative GO density (GO:0000302 / GO:0006979)
#     4) Scatter + Quadrant plots: promoter ΔH3Ac vs RNA log2FC (hda6NS vs Col0NS)
#     5) GO-set ΔH3Ac density (Heat vs ROS/OX)
#     6) HSF/HSP-focused panels (slope + dot + key-gene bars)
#
# Notes:
#   - No files are saved automatically (GitHub-friendly). Add ggsave() manually if needed.
#   - Requires annotation package org.At.tair.db (TAIR10).
#
# Inputs (edit if needed):
#   CHiP promoter quant TSV:
#     F:/Chip seq/05_quant/promoter_H3Ac_WT_vs_hda6*.tsv
#   Promoter BED (used if quant TSV lacks gene_id):
#     F:/Chip seq/ref/TAIR10_promoters_TSSpm1kb.clipped.bed
#   RNA-seq DESeq2 results:
#     F:/RNA seq/counts/DEG_hda6NS_vs_Col0NS_DESeq2_full.csv
#     Columns: GeneID, log2FoldChange, padj
############################################################

rm(list=ls())
options(stringsAsFactors=FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(AnnotationDbi)
  library(org.At.tair.db)
  library(ggrepel)
})

# ---------------------------
# PATHS (EDIT IF NEEDED)
# ---------------------------
CHIP_DIR  <- "F:/Chip seq"
QUANT_DIR <- file.path(CHIP_DIR, "05_quant")
PROM_BED  <- file.path(CHIP_DIR, "ref", "TAIR10_promoters_TSSpm1kb.clipped.bed")
RNA_FILE  <- "F:/RNA seq/counts/DEG_hda6NS_vs_Col0NS_DESeq2_full.csv"

stopifnot(dir.exists(QUANT_DIR))
stopifnot(file.exists(RNA_FILE))

# ---------------------------
# GLOBAL SETTINGS
# ---------------------------
PSEUDO <- 1e-6

COL_WT   <- "#2E8B57"  # Col-0 (WT)
COL_HDA6 <- "#6A5ACD"  # hda6
COL_DIR  <- c(
  "UP (padj<0.05)"   = "#6A5ACD",
  "DOWN (padj<0.05)" = "#2E8B57",
  "NS"               = "grey60"
)

BASE_SIZE <- 16
bold_theme <- theme_classic(base_size = BASE_SIZE) +
  theme(
    plot.title    = element_text(face="bold"),
    plot.subtitle = element_text(face="bold"),
    axis.title    = element_text(face="bold"),
    axis.text     = element_text(face="bold"),
    legend.title  = element_text(face="bold"),
    legend.text   = element_text(face="bold")
  )

# ---------------------------
# HELPERS
# ---------------------------
clean_gene_id <- function(x) {
  x <- as.character(x)
  x <- sub("\\..*$", "", x) # strip transcript suffix if present
  toupper(x)
}

pick_chip_quant_tsv <- function(quant_dir) {
  cand <- c(
    list.files(quant_dir, pattern="promoter_H3Ac_WT_vs_hda6.*Gene.*\\.tsv$", full.names=TRUE),
    list.files(quant_dir, pattern="promoter_H3Ac_WT_vs_hda6.*with.*\\.tsv$", full.names=TRUE),
    list.files(quant_dir, pattern="^promoter_H3Ac_WT_vs_hda6\\.tsv$", full.names=TRUE)
  ) %>% unique()
  if (length(cand) == 0) stop("No promoter_H3Ac_WT_vs_hda6*.tsv found in: ", quant_dir)
  cand[1]
}

detect_sig_cols <- function(df) {
  nm <- names(df)

  wt_col <- nm[
    str_detect(nm, regex("^WT(_|\\b)", ignore_case = FALSE)) &
      str_detect(nm, regex("RPGC|sig|signal", ignore_case = TRUE))
  ]
  hda6_col <- nm[
    str_detect(nm, regex("^hda6(_|\\b)", ignore_case = TRUE)) &
      str_detect(nm, regex("RPGC|sig|signal", ignore_case = TRUE))
  ]

  # fallbacks
  if (length(wt_col) == 0 && "WT_RPGC" %in% nm) wt_col <- "WT_RPGC"
  if (length(hda6_col) == 0 && "hda6_RPGC" %in% nm) hda6_col <- "hda6_RPGC"
  if (length(wt_col) == 0 && "WT_sig" %in% nm) wt_col <- "WT_sig"
  if (length(hda6_col) == 0 && "hda6_sig" %in% nm) hda6_col <- "hda6_sig"

  if (length(wt_col) == 0 || length(hda6_col) == 0) {
    stop("Cannot identify WT/hda6 signal columns.\nColumns:\n", paste(nm, collapse=", "))
  }
  list(wt = wt_col[1], hda6 = hda6_col[1])
}

ensure_gene_id <- function(chip_df, prom_bed) {
  if ("gene_id" %in% names(chip_df)) return(chip_df)

  if (!file.exists(prom_bed)) stop("Missing promoter bed: ", prom_bed)

  bed_raw <- read_tsv(prom_bed, col_names = FALSE, show_col_types = FALSE)
  bed <- bed_raw %>%
    transmute(chr = X1, start = as.integer(X2), end = as.integer(X3), gene_id = as.character(X4)) %>%
    distinct()

  chip_df %>%
    rename(
      chr   = any_of(c("chr","chrom","seqnames")),
      start = any_of(c("start","chromStart")),
      end   = any_of(c("end","chromEnd"))
    ) %>%
    mutate(start = as.integer(start), end = as.integer(end)) %>%
    left_join(bed, by = c("chr","start","end"))
}

map_tair_to_symbol <- function(tair_ids) {
  out <- AnnotationDbi::select(
    org.At.tair.db,
    keys = unique(tair_ids),
    keytype = "TAIR",
    columns = c("SYMBOL")
  ) %>%
    filter(!is.na(TAIR)) %>%
    distinct(TAIR, .keep_all = TRUE) %>%
    transmute(gene_id = clean_gene_id(TAIR), symbol = SYMBOL)
  out
}

get_GO_members <- function(go_ids, genes_tair) {
  out <- AnnotationDbi::select(
    org.At.tair.db,
    keys = unique(genes_tair),
    keytype = "TAIR",
    columns = c("GO")
  )
  out %>%
    filter(!is.na(GO), GO %in% go_ids) %>%
    pull(TAIR) %>%
    unique() %>%
    clean_gene_id()
}

symbol_to_tair <- function(symbols) {
  symbols <- unique(symbols[!is.na(symbols) & symbols != ""])
  if (length(symbols) == 0) return(tibble(symbol=character(0), gene_id=character(0)))
  out <- tryCatch(
    AnnotationDbi::select(
      org.At.tair.db,
      keys = symbols,
      keytype = "SYMBOL",
      columns = c("TAIR","SYMBOL")
    ),
    error = function(e) NULL
  )
  if (is.null(out)) return(tibble(symbol=character(0), gene_id=character(0)))
  out %>%
    filter(!is.na(TAIR), !is.na(SYMBOL)) %>%
    distinct(SYMBOL, TAIR) %>%
    transmute(symbol = SYMBOL, gene_id = clean_gene_id(TAIR))
}

# ---------------------------
# LOAD ChIP promoter quant (WT vs hda6)
# ---------------------------
chip_tsv <- pick_chip_quant_tsv(QUANT_DIR)
message("Using ChIP promoter quant TSV: ", chip_tsv)

chip_raw <- read_tsv(chip_tsv, show_col_types = FALSE)
chip_raw <- ensure_gene_id(chip_raw, PROM_BED)

sigcols <- detect_sig_cols(chip_raw)
message("WT signal column: ", sigcols$wt)
message("hda6 signal column: ", sigcols$hda6)

chip_prom <- chip_raw %>%
  transmute(
    gene_id  = clean_gene_id(gene_id),
    WT_sig   = suppressWarnings(as.numeric(.data[[sigcols$wt]])),
    hda6_sig = suppressWarnings(as.numeric(.data[[sigcols$hda6]]))
  ) %>%
  filter(!is.na(gene_id), !is.na(WT_sig), !is.na(hda6_sig)) %>%
  group_by(gene_id) %>%
  summarise(
    WT_sig   = mean(WT_sig,   na.rm=TRUE),
    hda6_sig = mean(hda6_sig, na.rm=TRUE),
    .groups="drop"
  ) %>%
  mutate(
    chip_ratio   = (hda6_sig + PSEUDO) / (WT_sig + PSEUDO),
    chip_log2rat = log2(chip_ratio)
  )

message("ChIP promoters (unique genes): ", nrow(chip_prom))

sym_map <- map_tair_to_symbol(chip_prom$gene_id)
chip_prom <- chip_prom %>%
  left_join(sym_map, by="gene_id") %>%
  mutate(gene_name = ifelse(!is.na(symbol) & symbol != "", symbol, gene_id))

# ---------------------------
# PLOT 1: Overall promoter H3Ac density (WT vs hda6)
# ---------------------------
df_all <- chip_prom %>%
  select(gene_id, WT_sig, hda6_sig) %>%
  pivot_longer(cols = c(WT_sig, hda6_sig), names_to = "Genotype", values_to = "Signal") %>%
  mutate(
    Genotype = recode(Genotype, WT_sig="Col-0", hda6_sig="hda6"),
    Signal = Signal + PSEUDO,
    log10Signal = log10(Signal)
  )

med_all <- df_all %>%
  group_by(Genotype) %>%
  summarise(med = median(log10Signal, na.rm=TRUE), .groups="drop")

p_density_all <- ggplot(df_all, aes(x = log10Signal, color = Genotype)) +
  geom_density(linewidth = 1.3) +
  scale_color_manual(values = c("Col-0"=COL_WT, "hda6"=COL_HDA6)) +
  geom_vline(data = med_all, aes(xintercept = med, color = Genotype),
             linetype="dashed", linewidth=0.9, show.legend = FALSE) +
  labs(
    title = "Overall promoter H3Ac levels (Col-0 vs hda6)",
    subtitle = "Density of log10(promoter H3Ac RPGC). Dashed lines indicate medians.",
    x = "log10 promoter H3Ac (RPGC)",
    y = "Density",
    color = NULL
  ) +
  bold_theme

print(p_density_all)

# ---------------------------
# PLOT 2: Heat GO density (GO:0009408)
# ---------------------------
GO_HEAT <- "GO:0009408"
heat_members <- get_GO_members(c(GO_HEAT), chip_prom$gene_id)

df_heat <- chip_prom %>%
  filter(gene_id %in% heat_members) %>%
  select(gene_id, WT_sig, hda6_sig) %>%
  pivot_longer(cols = c(WT_sig, hda6_sig), names_to = "Genotype", values_to = "Signal") %>%
  mutate(
    Genotype = recode(Genotype, WT_sig="Col-0", hda6_sig="hda6"),
    Signal = Signal + PSEUDO,
    log10Signal = log10(Signal)
  )

med_heat <- df_heat %>%
  group_by(Genotype) %>%
  summarise(med = median(log10Signal, na.rm=TRUE), .groups="drop")

p_density_heat <- ggplot(df_heat, aes(x = log10Signal, color = Genotype)) +
  geom_density(linewidth = 1.3) +
  scale_color_manual(values = c("Col-0"=COL_WT, "hda6"=COL_HDA6)) +
  geom_vline(data = med_heat, aes(xintercept = med, color = Genotype),
             linetype="dashed", linewidth=0.9, show.legend = FALSE) +
  labs(
    title = "Promoter H3Ac levels in heat-response genes (GO:0009408)",
    subtitle = "Density of log10(promoter H3Ac RPGC). Dashed lines indicate medians.",
    x = "log10 promoter H3Ac (RPGC)",
    y = "Density",
    color = NULL
  ) +
  bold_theme

print(p_density_heat)

# ---------------------------
# PLOT 3: ROS/OX GO density (GO:0000302 / GO:0006979)
# ---------------------------
GO_ROS <- "GO:0000302"
GO_OX  <- "GO:0006979"
ros_members <- get_GO_members(c(GO_ROS, GO_OX), chip_prom$gene_id)

df_ros <- chip_prom %>%
  filter(gene_id %in% ros_members) %>%
  select(gene_id, WT_sig, hda6_sig) %>%
  pivot_longer(cols = c(WT_sig, hda6_sig), names_to = "Genotype", values_to = "Signal") %>%
  mutate(
    Genotype = recode(Genotype, WT_sig="Col-0", hda6_sig="hda6"),
    Signal = Signal + PSEUDO,
    log10Signal = log10(Signal)
  )

med_ros <- df_ros %>%
  group_by(Genotype) %>%
  summarise(med = median(log10Signal, na.rm=TRUE), .groups="drop")

p_density_ros <- ggplot(df_ros, aes(x = log10Signal, color = Genotype)) +
  geom_density(linewidth = 1.3) +
  scale_color_manual(values = c("Col-0"=COL_WT, "hda6"=COL_HDA6)) +
  geom_vline(data = med_ros, aes(xintercept = med, color = Genotype),
             linetype="dashed", linewidth=0.9, show.legend = FALSE) +
  labs(
    title = "Promoter H3Ac levels in ROS/oxidative-stress genes (GO:0000302/GO:0006979)",
    subtitle = "Density of log10(promoter H3Ac RPGC). Dashed lines indicate medians.",
    x = "log10 promoter H3Ac (RPGC)",
    y = "Density",
    color = NULL
  ) +
  bold_theme

print(p_density_ros)

# ---------------------------
# LOAD RNA-seq DESeq2 results (hda6NS vs Col0NS)
# ---------------------------
rna_raw <- read_csv(RNA_FILE, show_col_types = FALSE)

req_cols <- c("GeneID","log2FoldChange","padj")
if (!all(req_cols %in% names(rna_raw))) {
  stop("RNA file missing required columns. Found: ", paste(names(rna_raw), collapse=", "))
}

rna <- rna_raw %>%
  transmute(
    gene_id    = clean_gene_id(.data[["GeneID"]]),
    rna_log2fc = suppressWarnings(as.numeric(.data[["log2FoldChange"]])),
    rna_padj   = suppressWarnings(as.numeric(.data[["padj"]]))
  ) %>%
  filter(!is.na(gene_id), !is.na(rna_log2fc)) %>%
  mutate(
    is_DEG = !is.na(rna_padj) & rna_padj < 0.05,
    dir = case_when(
      is_DEG & rna_log2fc > 0 ~ "UP (padj<0.05)",
      is_DEG & rna_log2fc < 0 ~ "DOWN (padj<0.05)",
      TRUE ~ "NS"
    )
  )

message("RNA genes: ", nrow(rna))

# ---------------------------
# MERGE ChIP ∩ RNA
# ---------------------------
m <- chip_prom %>% inner_join(rna, by="gene_id")
message("Overlap genes (ChIP ∩ RNA): ", nrow(m))

pear  <- cor.test(m$chip_log2rat, m$rna_log2fc, method="pearson")
spear <- cor.test(m$chip_log2rat, m$rna_log2fc, method="spearman")

cat("\n==============================\n")
cat("ΔH3Ac (promoter log2 hda6/Col-0) vs RNA log2FC (hda6NS vs Col0NS)\n")
cat("==============================\n")
cat("Pearson r = ", signif(pear$estimate, 4), " | p = ", signif(pear$p.value, 4), "\n", sep="")
cat("Spearman ρ = ", signif(spear$estimate, 4), " | p = ", signif(spear$p.value, 4), "\n", sep="")
cat("==============================\n\n")

lab <- m %>%
  mutate(score = abs(rna_log2fc) + abs(chip_log2rat)) %>%
  arrange(desc(score)) %>%
  slice_head(n = min(15, nrow(.)))

# ---------------------------
# PLOT 4: Scatter ΔH3Ac vs RNA log2FC
# ---------------------------
p_scatter <- ggplot(m, aes(x = chip_log2rat, y = rna_log2fc, color = dir)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.7) +
  geom_point(alpha = 0.75, size = 1.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  scale_color_manual(values = COL_DIR) +
  ggrepel::geom_text_repel(
    data = lab,
    aes(label = gene_name),
    size = 4,
    fontface = "bold",
    max.overlaps = 200,
    box.padding = 0.3,
    min.segment.length = 0
  ) +
  labs(
    title = "Promoter ΔH3Ac vs RNA expression change (hda6NS vs Col-0NS)",
    subtitle = paste0(
      "x = log2(hda6/Col-0) promoter H3Ac; y = RNA log2FC | Pearson r=", signif(pear$estimate,3),
      " (p=", signif(pear$p.value,3), "), Spearman ρ=", signif(spear$estimate,3),
      " (p=", signif(spear$p.value,3), ")"
    ),
    x = "Promoter ΔH3Ac (log2 hda6 / Col-0)",
    y = "RNA log2FC (hda6NS vs Col-0NS)",
    color = NULL
  ) +
  bold_theme

print(p_scatter)

# ---------------------------
# PLOT 5: Quadrant view (ΔH3Ac vs Δexpression; DEG emphasis)
# ---------------------------
p_quad <- ggplot(m, aes(x = chip_log2rat, y = rna_log2fc)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.7) +
  geom_point(aes(color = is_DEG), alpha = 0.75, size = 1.8) +
  scale_color_manual(values = c("TRUE"=COL_HDA6, "FALSE"="grey65")) +
  ggrepel::geom_text_repel(
    data = lab,
    aes(label = gene_name),
    size = 4,
    fontface = "bold",
    max.overlaps = 200,
    box.padding = 0.3,
    min.segment.length = 0
  ) +
  labs(
    title = "Quadrant view: promoter ΔH3Ac vs RNA log2FC",
    subtitle = "Top-right: increased H3Ac + increased expression in hda6; bottom-left: decreased H3Ac + decreased expression",
    x = "Promoter ΔH3Ac (log2 hda6 / Col-0)",
    y = "RNA log2FC (hda6NS vs Col-0NS)",
    color = "DEG (padj<0.05)"
  ) +
  bold_theme

print(p_quad)

# ---------------------------
# PLOT 6: ΔH3Ac density for stress-related GO sets (Heat vs ROS/OX)
# ---------------------------
heat_members_m <- get_GO_members(c(GO_HEAT), m$gene_id)
ros_members_m  <- get_GO_members(c(GO_ROS, GO_OX), m$gene_id)

m_go <- m %>%
  mutate(
    set = case_when(
      gene_id %in% heat_members_m ~ "HEAT (GO:0009408)",
      gene_id %in% ros_members_m  ~ "ROS/OX (GO:0000302/GO:0006979)",
      TRUE                        ~ "BACKGROUND"
    )
  )

w_heat <- wilcox.test(chip_log2rat ~ (set == "HEAT (GO:0009408)"), data = m_go)
w_ros  <- wilcox.test(chip_log2rat ~ (set == "ROS/OX (GO:0000302/GO:0006979)"), data = m_go)

cat("Wilcoxon ΔH3Ac shift vs background:\n")
cat(" - HEAT vs background: p=", signif(w_heat$p.value,4), "\n", sep="")
cat(" - ROS/OX vs background: p=", signif(w_ros$p.value,4), "\n\n", sep="")

p_den_go <- ggplot(m_go %>% filter(set != "BACKGROUND"), aes(x = chip_log2rat, color = set)) +
  geom_density(linewidth = 1.2) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.7) +
  labs(
    title = "ΔH3Ac distribution in stress-related gene sets",
    subtitle = "x = log2(hda6/Col-0) promoter H3Ac; dashed line at 0",
    x = "Promoter ΔH3Ac (log2 hda6 / Col-0)",
    y = "Density",
    color = NULL
  ) +
  bold_theme

print(p_den_go)

# ---------------------------
# HSF/HSP-FOCUSED PANELS
# ---------------------------
force_symbols <- c(
  "HSFA2","HSFA3","HSFA1A","HSFA1B","HSFA1D","HSFA1E",
  "HSFA4A","HSFA4C","HSFA5","HSFA6A","HSFA6B",
  "HSFA7A","HSFA7B","HSFA8",
  "HSFB1","HSFB2A","HSFB2B","HSFB3","HSFB4",
  "HSP101","HSP70","HSP90","HSP90.1","HSP90.2","HSP90.3","HSP90.4","HSP90.5",
  "HSP17.6","HSP17.6A","HSP17.6B","HSP17.6C","HSP17.4","HSP18.2","HSP23.5"
)

rx_family <- "^(HSF|HSP|HSC)"
force_map <- symbol_to_tair(force_symbols)

focus <- m %>%
  mutate(
    is_family = !is.na(gene_name) & str_detect(gene_name, rx_family),
    is_forced_symbol = !is.na(gene_name) & gene_name %in% force_symbols,
    is_forced_id = gene_id %in% force_map$gene_id
  ) %>%
  filter(is_family | is_forced_symbol | is_forced_id) %>%
  mutate(mean_sig = (WT_sig + hda6_sig)/2) %>%
  arrange(desc(mean_sig)) %>%
  group_by(gene_name) %>%
  slice_head(n=1) %>%
  ungroup()

message("HSF/HSP focus genes found in merged data: ", nrow(focus))

present_force <- sort(unique(as.character(focus$gene_name[focus$gene_name %in% force_symbols])))
missing_force <- setdiff(force_symbols, present_force)
cat("\n---- Forced symbols present (first 50) ----\n")
print(head(present_force, 50))
cat("\n---- Forced symbols MISSING (first 80) ----\n")
print(head(missing_force, 80))

if (nrow(focus) > 0) {

  focus <- focus %>%
    arrange(chip_log2rat) %>%
    mutate(gene_name = factor(gene_name, levels = unique(gene_name)))

  # Plot 7: Slope plot (promoter signal WT -> hda6)
  slope_long <- focus %>%
    select(gene_name, WT_sig, hda6_sig, chip_log2rat) %>%
    pivot_longer(cols = c(WT_sig, hda6_sig), names_to = "Genotype", values_to = "Signal") %>%
    mutate(
      Genotype = recode(Genotype, WT_sig="Col-0", hda6_sig="hda6"),
      Genotype = factor(Genotype, levels = c("Col-0","hda6")),
      Signal = Signal + PSEUDO
    )

  labA <- focus %>%
    arrange(desc(abs(chip_log2rat))) %>%
    slice_head(n = min(12, nrow(focus))) %>%
    pull(gene_name) %>% as.character()

  p_slope <- ggplot(slope_long, aes(x = Genotype, y = Signal, group = gene_name)) +
    geom_line(alpha = 0.55, linewidth = 0.7, color = "grey30") +
    geom_point(aes(color = Genotype), size = 2.4) +
    scale_color_manual(values = c("Col-0"=COL_WT, "hda6"=COL_HDA6)) +
    scale_y_log10(labels = label_number()) +
    ggrepel::geom_text_repel(
      data = slope_long %>% filter(Genotype=="hda6", gene_name %in% labA),
      aes(label = gene_name),
      size = 4,
      fontface = "bold",
      min.segment.length = 0,
      box.padding = 0.3,
      max.overlaps = 200
    ) +
    labs(
      title = "HSF/HSP genes: promoter H3Ac (slope plot, Col-0 → hda6)",
      subtitle = "Y-axis log10 (RPGC). Upward lines indicate higher promoter H3Ac in hda6.",
      x = NULL,
      y = "Promoter H3Ac (RPGC, log10)",
      color = NULL
    ) +
    bold_theme

  print(p_slope)

  # Plot 8: Two-panel dot plot (ΔH3Ac and RNA log2FC)
  dot_long <- focus %>%
    transmute(
      gene_name,
      `Promoter ΔH3Ac (log2 hda6/Col-0)` = chip_log2rat,
      `RNA log2FC (hda6NS vs Col-0NS)`  = rna_log2fc,
      is_DEG
    ) %>%
    pivot_longer(
      cols = c(`Promoter ΔH3Ac (log2 hda6/Col-0)`, `RNA log2FC (hda6NS vs Col-0NS)`),
      names_to = "Metric",
      values_to = "Value"
    ) %>%
    mutate(
      Metric = factor(Metric, levels = c("Promoter ΔH3Ac (log2 hda6/Col-0)", "RNA log2FC (hda6NS vs Col-0NS)"))
    )

  p_dot2 <- ggplot(dot_long, aes(x = Value, y = gene_name)) +
    geom_vline(
      data = data.frame(Metric=levels(dot_long$Metric)),
      aes(xintercept = 0),
      linetype="dashed",
      linewidth=0.8,
      inherit.aes = FALSE
    ) +
    geom_point(aes(shape = is_DEG), size = 2.6, color = COL_HDA6) +
    scale_shape_manual(values = c("TRUE"=16, "FALSE"=1)) +
    facet_wrap(~Metric, scales = "free_x", nrow = 1) +
    labs(
      title = "HSF/HSP genes: ΔH3Ac and Δexpression (same gene order)",
      subtitle = "Filled points = DEG (padj<0.05).",
      x = NULL,
      y = NULL,
      shape = "DEG"
    ) +
    bold_theme

  print(p_dot2)

  # Plot 9: Absolute promoter signal bars for key genes
  key_symbols <- c("HSFA2","HSFA3","HSFA7A","HSFA7B","HSP101")
  key_df <- focus %>%
    filter(as.character(gene_name) %in% key_symbols) %>%
    select(gene_name, WT_sig, hda6_sig) %>%
    pivot_longer(cols = c(WT_sig, hda6_sig), names_to = "Genotype", values_to = "Signal") %>%
    mutate(
      Genotype = recode(Genotype, WT_sig="Col-0", hda6_sig="hda6"),
      gene_name = factor(gene_name, levels = key_symbols),
      Signal = Signal + PSEUDO
    )

  if (nrow(key_df) > 0) {
    p_key <- ggplot(key_df, aes(x = gene_name, y = Signal, fill = Genotype)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7) +
      scale_fill_manual(values = c("Col-0"=COL_WT, "hda6"=COL_HDA6)) +
      scale_y_log10(labels = label_number()) +
      labs(
        title = "Key heat regulators: promoter H3Ac signal (Col-0 vs hda6)",
        subtitle = "Y-axis log10 (RPGC).",
        x = NULL,
        y = "Promoter H3Ac (RPGC, log10)",
        fill = NULL
      ) +
      bold_theme
    print(p_key)
  } else {
    message("None of the key symbols (HSFA2/HSFA3/HSFA7A/HSFA7B/HSP101) were found in the merged data.")
  }

} else {
  message("No HSF/HSP genes detected in merged data. Check SYMBOL mapping or promoter quant coverage.")
}

############################################################
# END
############################################################
