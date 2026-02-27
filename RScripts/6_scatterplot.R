## =========================
## 6_scatterplots.R
## =========================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
})

## -------------------- Parameters --------------------
analysis_level <- "gene"   # or "isoform" (must match your 4_categories output)
padj_cutoff <- 0.05

## -------------------- Paths --------------------
base_root <- getwd()
de_tbl   <- file.path(base_root, "results/1_DEAnalysis/tables/shrunk")
cat_tbl  <- file.path(base_root, "results/4_categories/tables")
fig_out  <- file.path(base_root, "results/6_scatterplots/figures")
dir.create(fig_out, recursive = TRUE, showWarnings = FALSE)

## -------------------- Drug list --------------------
drug_list <- c("11j", "KVS", "11j_PlaB", "KVS_PlaB", "PlaB")

## -------------------- Colors --------------------
col_gray <- "gray88"

cols_main <- c(
  enhanced_up       = "#D62728",
  enhanced_down     = "#1F77B4",
  suppressed_up     = "#FF7F0E",
  suppressed_down   = "#9467BD",
  switched_positive = "#2CA02C",
  switched_negative = "#17BECF"
)

cols_ind <- c(
  independent_up   = "#E31A1C",
  independent_down = "#1F78B4"
)

cols_shift <- c(
  shifted_baseline_enhanced_up     = "#D62728",
  shifted_baseline_enhanced_down   = "#1F77B4",
  shifted_baseline_suppressed_up   = "#FF7F0E",
  shifted_baseline_suppressed_down = "#9467BD"
)

## -------------------- Sizes (publication) --------------------
pt_other <- 1.2
pt_cat   <- 1.6

safe_num <- function(x) x[is.finite(x) & !is.na(x)]

pick_id_col <- function(df) {
  cand <- intersect(c("isoform_id", "gene_id", "feature_id"), names(df))
  if (length(cand) == 0) stop("No ID column found in DE table.")
  cand[1]
}

load_tt <- function(name) {
  fp <- file.path(de_tbl, paste0("tT_", name, ".tsv"))
  if (!file.exists(fp)) stop("Missing DE table: ", fp)
  read.delim(fp, check.names = FALSE)
}

load_ids <- function(tag, drug, analysis_level) {
  fn <- file.path(cat_tbl, paste0(tag, "_", drug, ".tsv"))
  if (!file.exists(fn)) return(character(0))
  df <- read.delim(fn, check.names = FALSE)
  
  if (analysis_level == "gene" && "gene_id" %in% names(df)) {
    unique(na.omit(df$gene_id))
  } else if (analysis_level == "isoform" && "isoform_id" %in% names(df)) {
    unique(na.omit(df$isoform_id))
  } else {
    character(0)
  }
}

make_merged_fc <- function(drug) {
  t_no  <- load_tt(drug)
  t_yes <- load_tt(paste0(drug, "_OHT"))
  
  id_no  <- pick_id_col(t_no)
  id_yes <- pick_id_col(t_yes)
  t_no  <- dplyr::rename(t_no,  feature_id = !!id_no)
  t_yes <- dplyr::rename(t_yes, feature_id = !!id_yes)
  
  full_join(
    t_no  %>% dplyr::select(feature_id, log2FoldChange, padj) %>%
      dplyr::rename(log2FC_noOHT = log2FoldChange, padj_noOHT = padj),
    t_yes %>% dplyr::select(feature_id, log2FoldChange, padj) %>%
      dplyr::rename(log2FC_withOHT = log2FoldChange, padj_withOHT = padj),
    by = "feature_id"
  )
}

axis_x <- function(drug) {
  bquote(bold(atop(.(drug) ~ "/ DMSO", "(" * log[2] * "FC" * ")")))
}
axis_y <- function(drug) {
  bquote(bold(atop(.(drug) ~ "+ OHT / DMSO + OHT", "(" * log[2] * "FC" * ")")))
}

theme_scatter <- function() {
  theme_bw(base_size = 14) +
    theme(
      panel.grid   = element_blank(),
      plot.title   = element_text(face = "bold", hjust = 0.5, size = 16),
      axis.title   = element_text(face = "bold", size = 15),
      axis.text    = element_text(face = "bold", size = 13),
      legend.title = element_text(face = "bold", size = 14),
      legend.text  = element_text(face = "bold", size = 12),
      legend.key.size   = unit(0.6, "cm"),
      legend.spacing.x  = unit(0.4, "cm"),
      legend.position   = "right"
    )
}

save_both <- function(filename_base, plot_obj, width, height, dpi = 600) {
  ggsave(paste0(filename_base, ".pdf"),  plot_obj, width = width, height = height)
  ggsave(paste0(filename_base, ".tiff"), plot_obj,
         width = width, height = height, dpi = dpi,
         device = "tiff", compression = "lzw")
}

get_legend <- function(p) {
  g <- ggplotGrob(p)
  idx <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
  if (length(idx) == 0) return(NULL)
  g$grobs[[idx[1]]]
}

make_vertical_legend <- function(p, title = NULL) {
  p2 <- p +
    guides(color = guide_legend(ncol = 1, byrow = TRUE)) +
    theme(
      legend.position  = "right",
      legend.direction = "vertical",
      legend.box       = "vertical"
    )
  if (!is.null(title)) p2 <- p2 + labs(color = title)
  p2
}
## -------------------- Global limits --------------------
scatter_list <- lapply(drug_list, make_merged_fc)
names(scatter_list) <- drug_list

all_x <- safe_num(unlist(lapply(scatter_list, function(d) d$log2FC_noOHT),  use.names = FALSE))
all_y <- safe_num(unlist(lapply(scatter_list, function(d) d$log2FC_withOHT), use.names = FALSE))

lim_max <- max(abs(c(all_x, all_y)), na.rm = TRUE)
xy_lim <- c(-lim_max, lim_max)

## -------------------- Collect for combined panels --------------------
p_main_list  <- list()
p_ind_list   <- list()
p_shift_list <- list()

for (drug in drug_list) {
  message("Plotting: ", drug)
  
  df <- scatter_list[[drug]]
  df$category_main  <- "Other"
  df$category_ind   <- "Other"
  df$category_shift <- "Other"
  
  enh_up   <- load_ids("enhanced_up",       drug, analysis_level)
  enh_down <- load_ids("enhanced_down",     drug, analysis_level)
  sup_up   <- load_ids("suppressed_up",     drug, analysis_level)
  sup_down <- load_ids("suppressed_down",   drug, analysis_level)
  swi_pos  <- load_ids("switched_positive", drug, analysis_level)
  swi_neg  <- load_ids("switched_negative", drug, analysis_level)
  
  ind_up   <- load_ids("independent_up",    drug, analysis_level)
  ind_down <- load_ids("independent_down",  drug, analysis_level)
  
  sh_enh_up   <- load_ids("shifted_baseline_enhanced_up",     drug, analysis_level)
  sh_enh_down <- load_ids("shifted_baseline_enhanced_down",   drug, analysis_level)
  sh_sup_up   <- load_ids("shifted_baseline_suppressed_up",   drug, analysis_level)
  sh_sup_down <- load_ids("shifted_baseline_suppressed_down", drug, analysis_level)
  
  df$category_main[df$feature_id %in% enh_up]   <- "enhanced_up"
  df$category_main[df$feature_id %in% enh_down] <- "enhanced_down"
  df$category_main[df$feature_id %in% sup_up]   <- "suppressed_up"
  df$category_main[df$feature_id %in% sup_down] <- "suppressed_down"
  df$category_main[df$feature_id %in% swi_pos]  <- "switched_positive"
  df$category_main[df$feature_id %in% swi_neg]  <- "switched_negative"
  df$category_main <- factor(df$category_main, levels = c("Other", names(cols_main)))
  
  df$category_ind[df$feature_id %in% ind_up]   <- "independent_up"
  df$category_ind[df$feature_id %in% ind_down] <- "independent_down"
  df$category_ind <- factor(df$category_ind, levels = c("Other", names(cols_ind)))
  
  df$category_shift[df$feature_id %in% sh_enh_up]   <- "shifted_baseline_enhanced_up"
  df$category_shift[df$feature_id %in% sh_enh_down] <- "shifted_baseline_enhanced_down"
  df$category_shift[df$feature_id %in% sh_sup_up]   <- "shifted_baseline_suppressed_up"
  df$category_shift[df$feature_id %in% sh_sup_down] <- "shifted_baseline_suppressed_down"
  df$category_shift <- factor(df$category_shift, levels = c("Other", names(cols_shift)))
  
  p1 <- ggplot(df, aes(x = log2FC_noOHT, y = log2FC_withOHT)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) +
    geom_point(
      data  = df[df$category_main == "Other" | is.na(df$category_main), , drop = FALSE],
      color = col_gray, alpha = 0.55, size = pt_other
    ) +
    geom_point(
      data = df[df$category_main != "Other" & !is.na(df$category_main), , drop = FALSE],
      aes(color = category_main),
      alpha = 0.9, size = pt_cat
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    scale_color_manual(values = cols_main, drop = FALSE, name = "Category") +
    labs(title = drug, x = axis_x(drug), y = axis_y(drug)) +
    coord_fixed(ratio = 1) +
    coord_cartesian(xlim = xy_lim, ylim = xy_lim) +
    theme_scatter()
  
  p2 <- ggplot(df, aes(x = log2FC_noOHT, y = log2FC_withOHT)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) +
    geom_point(
      data  = df[df$category_ind == "Other" | is.na(df$category_ind), , drop = FALSE],
      color = col_gray, alpha = 0.55, size = pt_other
    ) +
    geom_point(
      data = df[df$category_ind != "Other" & !is.na(df$category_ind), , drop = FALSE],
      aes(color = category_ind),
      alpha = 0.9, size = pt_cat
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    scale_color_manual(values = cols_ind, drop = FALSE, name = "Independent") +
    labs(title = drug, x = axis_x(drug), y = axis_y(drug)) +
    coord_fixed(ratio = 1) +
    coord_cartesian(xlim = xy_lim, ylim = xy_lim) +
    theme_scatter()
  
  p3 <- ggplot(df, aes(x = log2FC_noOHT, y = log2FC_withOHT)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) +
    geom_point(
      data  = df[df$category_shift == "Other" | is.na(df$category_shift), , drop = FALSE],
      color = col_gray, alpha = 0.55, size = pt_other
    ) +
    geom_point(
      data = df[df$category_shift != "Other" & !is.na(df$category_shift), , drop = FALSE],
      aes(color = category_shift),
      alpha = 0.9, size = pt_cat
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    scale_color_manual(values = cols_shift, drop = FALSE, name = "Shifted") +
    labs(title = drug, x = axis_x(drug), y = axis_y(drug)) +
    coord_fixed(ratio = 1) +
    coord_cartesian(xlim = xy_lim, ylim = xy_lim) +
    theme_scatter()
  
  ## per-drug true squares (no legend)
  save_both(file.path(fig_out, drug),                 p1 + theme(legend.position = "none"), 6, 6)
  save_both(file.path(fig_out, paste0(drug, "_independent")), p2 + theme(legend.position = "none"), 6, 6)
  save_both(file.path(fig_out, paste0(drug, "_shifted")),     p3 + theme(legend.position = "none"), 6, 6)
  
  ## keep for combined panels
  p_main_list[[drug]]  <- p1
  p_ind_list[[drug]]   <- p2
  p_shift_list[[drug]] <- p3
}


comb_main <- wrap_plots(p_main_list, ncol = 3) &
  theme(legend.position = "none")

save_both(file.path(fig_out, "ALL_main_categories"), comb_main, 18, 12)

comb_ind <- wrap_plots(p_ind_list, ncol = 3) &
  theme(legend.position = "none")

save_both(file.path(fig_out, "ALL_independent"), comb_ind, 18, 12)

comb_shift <- wrap_plots(p_shift_list, ncol = 3) &
  theme(legend.position = "none")

save_both(file.path(fig_out, "ALL_shifted_baseline"), comb_shift, 18, 12)

## ---- MAIN legend only (vertical) ----
p_rep_main <- make_vertical_legend(p_main_list[[drug_list[1]]], title = "Category")
leg_main <- get_legend(p_rep_main)
if (!is.null(leg_main)) {
  save_both(
    file.path(fig_out, "LEGEND_main_categories"),
    patchwork::wrap_elements(full = leg_main),
    width = 4, height = 8, dpi = 600
  )
}

## ---- INDEPENDENT legend only (vertical) ----
p_rep_ind <- make_vertical_legend(p_ind_list[[drug_list[1]]], title = "Independent")
leg_ind <- get_legend(p_rep_ind)
if (!is.null(leg_ind)) {
  save_both(
    file.path(fig_out, "LEGEND_independent"),
    patchwork::wrap_elements(full = leg_ind),
    width = 4, height = 5, dpi = 600
  )
}

## ---- SHIFTED legend only (vertical) ----
p_rep_shift <- make_vertical_legend(p_shift_list[[drug_list[1]]], title = "Shifted")
leg_shift <- get_legend(p_rep_shift)
if (!is.null(leg_shift)) {
  save_both(
    file.path(fig_out, "LEGEND_shifted_baseline"),
    patchwork::wrap_elements(full = leg_shift),
    width = 4, height = 7, dpi = 600
  )
}
message("Done: ", fig_out)