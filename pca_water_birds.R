# pca_water_birds.R
# ============================================================================
# PCA to identify which water variables drive bird counts
# Runs for:
#   1. Total count (all species combined)
#   2. Individual species
#   3. By region
# ============================================================================

# =============================================================================
# PACKAGES
# =============================================================================

required_packages <- c("dplyr", "tidyr", "ggplot2",
                       "scales", "factoextra",
                       "ggrepel", "patchwork",
                       "tsibble", "readr")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

filter    <- dplyr::filter
mutate    <- dplyr::mutate
summarise <- dplyr::summarise
select    <- dplyr::select

# =============================================================================
# LOAD DATA
# =============================================================================

if (!exists("data_sub")) {
  if (!exists("get_wading_bird_data")) source("data_functions.R")
  
  piper_config <- list(
    spatial = list(
      level              = "subregion",
      fill_missing       = TRUE,
      fill_value         = 0,
      min_years_required = 10,
      exclude_regions    = list(),
      include_regions    = list(),
      exclude_colonies   = list(),
      include_colonies   = list(),
      run_by_region      = FALSE
    )
  )
  
  if (!dir.exists("SiteandMethods")) {
    library(wader)
    download_observations(".")
  }
  data_sub <- get_wading_bird_data(piper_config)
}

# =============================================================================
# WATER VARIABLES
# =============================================================================

water_vars <- c(
  "breed_season_depth",
  "dry_days",
  "recession",
  "init_depth",
  "pre_recession",
  "post_recession",
  "reversals"
)

species_list <- c("gbhe", "greg", "rosp", "sneg", "whib", "wost")

# Species common names for labels
species_names <- c(
  gbhe = "Great Blue Heron",
  greg = "Great Egret",
  rosp = "Roseate Spoonbill",
  sneg = "Snowy Egret",
  whib = "White Ibis",
  wost = "Wood Stork"
)

# =============================================================================
# PREPARE DATA
# =============================================================================

# Base: region x year with water covariates + counts per species
base_data <- data_sub |>
  as_tibble() |>
  dplyr::group_by(year, region, species) |>
  dplyr::summarise(
    count              = sum(count, na.rm = TRUE),
    across(all_of(water_vars), ~mean(., na.rm = TRUE)),
    .groups = "drop"
  )

# Wide format: one row per year x region, species as columns
wide_data <- base_data |>
  dplyr::select(year, region, species, count,
                all_of(water_vars)) |>
  tidyr::pivot_wider(
    names_from  = species,
    values_from = count,
    names_prefix = "count_"
  ) |>
  dplyr::mutate(
    total_count = rowSums(
      dplyr::select(dplyr::cur_data(),
                    starts_with("count_")),
      na.rm = TRUE
    )
  )

cat("Data prepared:", nrow(wide_data), "rows\n")
cat("Regions:", paste(unique(wide_data$region), collapse = ", "), "\n")
cat("Years:", min(wide_data$year), "-", max(wide_data$year), "\n\n")

# =============================================================================
# PCA FUNCTION
# =============================================================================

run_pca <- function(df,
                    response_var,
                    label        = NULL,
                    scale_vars   = TRUE) {
  
  # Remove rows with NA in water vars or response
  df_clean <- df |>
    dplyr::select(all_of(c(water_vars, response_var,
                           "year", "region"))) |>
    tidyr::drop_na()
  
  if (nrow(df_clean) < 5) {
    cat("  Not enough data for:", label %||% response_var, "\n")
    return(NULL)
  }
  
  # Scale water variables
  X <- df_clean |>
    dplyr::select(all_of(water_vars)) |>
    as.matrix()
  
  y <- df_clean[[response_var]]
  
  # Run PCA
  pca_result <- prcomp(X, scale. = scale_vars, center = TRUE)
  
  # Variance explained
  var_exp <- summary(pca_result)$importance
  pc1_var <- round(var_exp[2, 1] * 100, 1)
  pc2_var <- round(var_exp[2, 2] * 100, 1)
  pc3_var <- round(var_exp[2, 3] * 100, 1)
  
  cat(sprintf("  %s: PC1=%.1f%% PC2=%.1f%% PC3=%.1f%%\n",
              label %||% response_var,
              pc1_var, pc2_var, pc3_var))
  
  # Loadings
  loadings <- as.data.frame(pca_result$rotation) |>
    tibble::rownames_to_column("variable") |>
    dplyr::mutate(
      variable = dplyr::recode(variable,
                               breed_season_depth = "Breed Depth",
                               dry_days           = "Dry Days",
                               recession          = "Recession",
                               init_depth         = "Init Depth",
                               pre_recession      = "Pre-Recession",
                               post_recession     = "Post-Recession",
                               reversals          = "Reversals"
      )
    )
  
  # Scores
  scores <- as.data.frame(pca_result$x) |>
    dplyr::mutate(
      response = y,
      year     = df_clean$year,
      region   = df_clean$region
    )
  
  # Correlation of PCs with response
  pc_cors <- sapply(1:min(5, ncol(scores) - 3), function(i) {
    cor(scores[[paste0("PC", i)]], scores$response,
        use = "complete.obs")
  })
  
  list(
    pca      = pca_result,
    loadings = loadings,
    scores   = scores,
    var_exp  = var_exp,
    pc_cors  = pc_cors,
    label    = label %||% response_var,
    pc1_var  = pc1_var,
    pc2_var  = pc2_var,
    n        = nrow(df_clean)
  )
}

# =============================================================================
# PLOT FUNCTIONS
# =============================================================================

plot_biplot <- function(pca_result, title = NULL) {
  
  scores   <- pca_result$scores
  loadings <- pca_result$loadings
  pc1_var  <- pca_result$pc1_var
  pc2_var  <- pca_result$pc2_var
  
  # Scale loadings for overlay
  scale_factor <- max(abs(scores[, c("PC1", "PC2")])) /
    max(abs(loadings[, c("PC1", "PC2")])) * 0.7
  
  loading_arrows <- loadings |>
    dplyr::mutate(
      x    = 0,
      y    = 0,
      xend = PC1 * scale_factor,
      yend = PC2 * scale_factor
    )
  
  p <- ggplot2::ggplot() +
    
    # Points colored by response
    ggplot2::geom_point(
      data = scores,
      ggplot2::aes(x = PC1, y = PC2, colour = response),
      size = 2.5, alpha = 0.7
    ) +
    ggplot2::scale_colour_gradient2(
      low      = "#3498db",
      mid      = "#f7dc6f",
      high     = "#e74c3c",
      midpoint = median(scores$response, na.rm = TRUE),
      name     = "Count",
      labels   = scales::comma
    ) +
    
    # Loading arrows
    ggplot2::geom_segment(
      data = loading_arrows,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      arrow     = ggplot2::arrow(length = ggplot2::unit(0.25, "cm"),
                                 type   = "closed"),
      colour    = "grey20",
      linewidth = 0.8
    ) +
    
    # Loading labels
    ggrepel::geom_text_repel(
      data = loading_arrows,
      ggplot2::aes(x = xend * 1.1, y = yend * 1.1,
                   label = variable),
      size        = 3.5,
      colour      = "grey20",
      fontface    = "bold",
      box.padding = 0.3,
      max.overlaps = 20
    ) +
    
    # Zero lines
    ggplot2::geom_hline(yintercept = 0,
                        linetype = "dashed", colour = "grey70") +
    ggplot2::geom_vline(xintercept = 0,
                        linetype = "dashed", colour = "grey70") +
    
    ggplot2::labs(
      title = title %||% pca_result$label,
      x     = sprintf("PC1 (%.1f%% variance)", pc1_var),
      y     = sprintf("PC2 (%.1f%% variance)", pc2_var)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(face = "bold", size = 11),
      legend.position = "right"
    )
  
  return(p)
}

plot_loadings <- function(pca_result, n_pcs = 3, title = NULL) {
  
  loadings_long <- pca_result$loadings |>
    tidyr::pivot_longer(
      cols      = starts_with("PC"),
      names_to  = "PC",
      values_to = "loading"
    ) |>
    dplyr::filter(PC %in% paste0("PC", 1:n_pcs)) |>
    dplyr::mutate(
      PC = factor(PC, levels = paste0("PC", 1:n_pcs)),
      direction = ifelse(loading > 0, "Positive", "Negative")
    )
  
  # Add variance explained to PC labels
  var_exp <- pca_result$var_exp
  pc_labs <- sapply(1:n_pcs, function(i) {
    sprintf("PC%d\n(%.1f%%)", i, var_exp[2, i] * 100)
  })
  names(pc_labs) <- paste0("PC", 1:n_pcs)
  
  p <- ggplot2::ggplot(
    loadings_long,
    ggplot2::aes(x    = reorder(variable, abs(loading)),
                 y    = loading,
                 fill = direction)
  ) +
    ggplot2::geom_col(alpha = 0.85) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40") +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~PC,
                        ncol    = n_pcs,
                        labeller = ggplot2::as_labeller(pc_labs)) +
    ggplot2::scale_fill_manual(
      values = c("Positive" = "#e74c3c",
                 "Negative" = "#3498db"),
      name   = "Direction"
    ) +
    ggplot2::labs(
      title = title %||% paste(pca_result$label, "- Loadings"),
      x     = "Water Variable",
      y     = "Loading"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(face = "bold", size = 11),
      strip.text      = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )
  
  return(p)
}

plot_pc_correlations <- function(pca_list, title = NULL) {
  
  # Combine PC-response correlations across all analyses
  cor_df <- do.call(rbind, lapply(names(pca_list), function(nm) {
    res <- pca_list[[nm]]
    if (is.null(res)) return(NULL)
    data.frame(
      analysis = res$label,
      PC       = paste0("PC", seq_along(res$pc_cors)),
      cor      = res$pc_cors,
      stringsAsFactors = FALSE
    )
  }))
  
  cor_df <- cor_df |>
    dplyr::filter(!is.na(cor)) |>
    dplyr::mutate(
      abs_cor   = abs(cor),
      direction = ifelse(cor > 0, "Positive", "Negative")
    )
  
  p <- ggplot2::ggplot(
    cor_df,
    ggplot2::aes(x    = PC,
                 y    = reorder(analysis, abs_cor),
                 fill = cor)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
    ggplot2::geom_text(
      ggplot2::aes(label = round(cor, 2)),
      size = 3, colour = "black"
    ) +
    ggplot2::scale_fill_gradient2(
      low      = "#3498db",
      mid      = "white",
      high     = "#e74c3c",
      midpoint = 0,
      limits   = c(-1, 1),
      name     = "Correlation\nwith Count"
    ) +
    ggplot2::labs(
      title = title %||% "PC Correlation with Bird Counts",
      x     = "Principal Component",
      y     = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(face = "bold", size = 12),
      axis.text.y     = ggplot2::element_text(size = 9),
      legend.position = "right"
    )
  
  return(p)
}

plot_variance_explained <- function(pca_list, title = NULL) {
  
  var_df <- do.call(rbind, lapply(names(pca_list), function(nm) {
    res <- pca_list[[nm]]
    if (is.null(res)) return(NULL)
    ve  <- res$var_exp[2, 1:min(5, ncol(res$var_exp))] * 100
    data.frame(
      analysis    = res$label,
      PC          = paste0("PC", seq_along(ve)),
      var_exp_pct = as.numeric(ve),
      stringsAsFactors = FALSE
    )
  }))
  
  p <- ggplot2::ggplot(
    var_df,
    ggplot2::aes(x    = PC,
                 y    = var_exp_pct,
                 fill = analysis,
                 group = analysis)
  ) +
    ggplot2::geom_line(
      ggplot2::aes(colour = analysis),
      linewidth = 1
    ) +
    ggplot2::geom_point(
      ggplot2::aes(colour = analysis),
      size = 3
    ) +
    ggplot2::scale_colour_brewer(palette = "Set2", name = NULL) +
    ggplot2::labs(
      title = title %||% "Variance Explained by PC",
      x     = "Principal Component",
      y     = "% Variance Explained"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(face = "bold", size = 12),
      legend.position = "bottom",
      legend.text     = ggplot2::element_text(size = 8)
    )
  
  return(p)
}

# =============================================================================
# RUN PCA
# =============================================================================

out_dir <- "results/pca_plots"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cat("=== Running PCA analyses ===\n\n")

# ==========================================================================
# 1. TOTAL COUNT
# ==========================================================================

cat("1. Total count...\n")

pca_total <- run_pca(
  df           = wide_data,
  response_var = "total_count",
  label        = "Total (All Species)"
)

if (!is.null(pca_total)) {
  p_bi <- plot_biplot(pca_total,
                      title = "PCA Biplot - Total Bird Count")
  p_lo <- plot_loadings(pca_total,
                        title = "Water Variable Loadings - Total Count")
  
  ggplot2::ggsave(file.path(out_dir, "pca_total_biplot.png"),
                  p_bi, width = 10, height = 7, dpi = 300, bg = "white")
  ggplot2::ggsave(file.path(out_dir, "pca_total_loadings.png"),
                  p_lo, width = 12, height = 5, dpi = 300, bg = "white")
  cat("  Saved\n\n")
}

# ==========================================================================
# 2. INDIVIDUAL SPECIES
# ==========================================================================

cat("2. Individual species...\n")

pca_species <- list()

for (sp in species_list) {
  resp <- paste0("count_", sp)
  
  if (!resp %in% names(wide_data)) next
  
  pca_species[[sp]] <- run_pca(
    df           = wide_data,
    response_var = resp,
    label        = species_names[sp]
  )
}

# Biplot grid - all species
cat("  Saving species biplot grid...\n")
biplot_list <- lapply(species_list, function(sp) {
  res <- pca_species[[sp]]
  if (is.null(res)) return(NULL)
  plot_biplot(res, title = species_names[sp]) +
    ggplot2::theme(legend.position = "none",
                   plot.title = ggplot2::element_text(size = 9))
})
biplot_list <- Filter(Negate(is.null), biplot_list)

if (length(biplot_list) > 0) {
  grid_bi <- gridExtra::grid.arrange(
    grobs = biplot_list,
    ncol  = 3,
    top   = grid::textGrob(
      "PCA Biplots by Species",
      gp = grid::gpar(fontsize = 14, fontface = "bold")
    )
  )
  ggplot2::ggsave(
    file.path(out_dir, "pca_species_biplots.png"),
    grid_bi, width = 18, height = 12, dpi = 250, bg = "white"
  )
}

# Loading grid - all species
cat("  Saving species loadings grid...\n")
loading_list <- lapply(species_list, function(sp) {
  res <- pca_species[[sp]]
  if (is.null(res)) return(NULL)
  plot_loadings(res, n_pcs = 2, title = species_names[sp]) +
    ggplot2::theme(legend.position = "none",
                   plot.title = ggplot2::element_text(size = 9))
})
loading_list <- Filter(Negate(is.null), loading_list)

if (length(loading_list) > 0) {
  grid_lo <- gridExtra::grid.arrange(
    grobs = loading_list,
    ncol  = 3,
    top   = grid::textGrob(
      "Water Variable Loadings by Species",
      gp = grid::gpar(fontsize = 14, fontface = "bold")
    )
  )
  ggplot2::ggsave(
    file.path(out_dir, "pca_species_loadings.png"),
    grid_lo, width = 18, height = 14, dpi = 250, bg = "white"
  )
}
cat("  Saved\n\n")

# ==========================================================================
# 3. BY REGION
# ==========================================================================

cat("3. By region...\n")

regions      <- sort(unique(wide_data$region))
pca_regions  <- list()

for (reg in regions) {
  df_reg <- wide_data |>
    dplyr::filter(region == reg)
  
  pca_regions[[reg]] <- run_pca(
    df           = df_reg,
    response_var = "total_count",
    label        = reg
  )
}

# Biplot grid - all regions
cat("  Saving region biplot grid...\n")
region_biplot_list <- lapply(regions, function(reg) {
  res <- pca_regions[[reg]]
  if (is.null(res)) return(NULL)
  plot_biplot(res, title = reg) +
    ggplot2::theme(legend.position = "none",
                   plot.title = ggplot2::element_text(size = 9))
})
region_biplot_list <- Filter(Negate(is.null), region_biplot_list)

if (length(region_biplot_list) > 0) {
  grid_reg <- gridExtra::grid.arrange(
    grobs = region_biplot_list,
    ncol  = 4,
    top   = grid::textGrob(
      "PCA Biplots by Region (Total Count)",
      gp = grid::gpar(fontsize = 14, fontface = "bold")
    )
  )
  ggplot2::ggsave(
    file.path(out_dir, "pca_region_biplots.png"),
    grid_reg, width = 20, height = 12, dpi = 250, bg = "white"
  )
}

# Loading grid - all regions
cat("  Saving region loadings grid...\n")
region_loading_list <- lapply(regions, function(reg) {
  res <- pca_regions[[reg]]
  if (is.null(res)) return(NULL)
  plot_loadings(res, n_pcs = 2, title = reg) +
    ggplot2::theme(legend.position = "none",
                   plot.title = ggplot2::element_text(size = 9))
})
region_loading_list <- Filter(Negate(is.null), region_loading_list)

if (length(region_loading_list) > 0) {
  grid_reg_lo <- gridExtra::grid.arrange(
    grobs = region_loading_list,
    ncol  = 4,
    top   = grid::textGrob(
      "Water Variable Loadings by Region",
      gp = grid::gpar(fontsize = 14, fontface = "bold")
    )
  )
  ggplot2::ggsave(
    file.path(out_dir, "pca_region_loadings.png"),
    grid_reg_lo, width = 20, height = 12, dpi = 250, bg = "white"
  )
}
cat("  Saved\n\n")

# ==========================================================================
# 4. SUMMARY PLOTS
# ==========================================================================

cat("4. Summary plots...\n")

# Combine all PCAs for comparison
all_pcas <- c(
  list(total = pca_total),
  pca_species,
  pca_regions
)
all_pcas <- Filter(Negate(is.null), all_pcas)

# PC correlation heatmap
p_cor <- plot_pc_correlations(
  all_pcas,
  title = "PC Correlation with Bird Counts (All Analyses)"
)
ggplot2::ggsave(
  file.path(out_dir, "pca_correlation_heatmap.png"),
  p_cor, width = 10, height = 12, dpi = 300, bg = "white"
)

# Variance explained comparison
p_var <- plot_variance_explained(
  c(list(total = pca_total), pca_species),
  title = "Variance Explained: Total vs Species"
)
ggplot2::ggsave(
  file.path(out_dir, "pca_variance_explained.png"),
  p_var, width = 10, height = 6, dpi = 300, bg = "white"
)

cat("  Saved\n\n")

# ==========================================================================
# 5. SUMMARY TABLE
# ==========================================================================

cat("=== Summary: Most Important Water Variables ===\n\n")

# For each analysis, which variable has the highest PC1 loading?
summary_df <- do.call(rbind, lapply(names(all_pcas), function(nm) {
  res <- all_pcas[[nm]]
  if (is.null(res)) return(NULL)
  
  # Top variables by absolute PC1 loading
  top_vars <- res$loadings |>
    dplyr::arrange(desc(abs(PC1))) |>
    dplyr::slice_head(n = 3)
  
  data.frame(
    analysis    = res$label,
    top_var_1   = top_vars$variable[1],
    loading_1   = round(top_vars$PC1[1], 3),
    top_var_2   = top_vars$variable[2],
    loading_2   = round(top_vars$PC1[2], 3),
    top_var_3   = top_vars$variable[3],
    loading_3   = round(top_vars$PC1[3], 3),
    pc1_var_pct = res$pc1_var,
    pc1_cor     = round(res$pc_cors[1], 3),
    n_obs       = res$n,
    stringsAsFactors = FALSE
  )
}))

cat("Top water variables driving PC1 by analysis:\n\n")
print(summary_df, row.names = FALSE)

# Most consistently important variables across all analyses
cat("\n=== Most Consistently Important Variables ===\n")
var_importance <- summary_df |>
  dplyr::select(top_var_1, top_var_2, top_var_3) |>
  unlist() |>
  table() |>
  sort(decreasing = TRUE)

print(var_importance)

cat("\n=== Complete ===\n")
cat("Plots saved to:", out_dir, "\n\n")
print(list.files(out_dir, pattern = "\\.png$"))






# variable selector -------------------------------------------------------



# pca_variable_selection.R
# ============================================================================
# Extract and visualise the most important water variables
# from PCA results across total count, species, and regions
# ============================================================================

# =============================================================================
# PACKAGES
# =============================================================================

required_packages <- c("dplyr", "tidyr", "ggplot2",
                       "scales", "ggrepel", "forcats",
                       "tsibble", "readr")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

filter    <- dplyr::filter
mutate    <- dplyr::mutate
summarise <- dplyr::summarise
select    <- dplyr::select

# =============================================================================
# LOAD DATA AND RUN PCA IF NOT ALREADY DONE
# =============================================================================

if (!exists("pca_total") || !exists("pca_species") ||
    !exists("pca_regions")) {
  cat("Running PCA first...\n")
  source("pca_water_birds.R")
}

# =============================================================================
# EXTRACT VARIABLE IMPORTANCE
# Uses two criteria:
#   1. Loading magnitude on PC1 and PC2
#   2. Correlation of that PC with bird count
#   Combined = weighted importance score
# =============================================================================

extract_importance <- function(pca_result) {
  
  if (is.null(pca_result)) return(NULL)
  
  # Weight loadings by PC-count correlation
  # If PC1 correlates strongly with count, PC1 loadings matter more
  pc1_cor <- abs(pca_result$pc_cors[1])
  pc2_cor <- abs(pca_result$pc_cors[2])
  pc3_cor <- if (length(pca_result$pc_cors) >= 3) {
    abs(pca_result$pc_cors[3])
  } else 0
  
  pca_result$loadings |>
    dplyr::select(variable, PC1, PC2, PC3) |>
    dplyr::mutate(
      # Weighted importance score
      importance = abs(PC1) * pc1_cor +
        abs(PC2) * pc2_cor +
        abs(PC3) * pc3_cor,
      # Raw loadings weighted by PC variance explained
      pc1_wt = abs(PC1) * pca_result$var_exp[2, 1],
      pc2_wt = abs(PC2) * pca_result$var_exp[2, 2],
      combined_loading = pc1_wt + pc2_wt,
      # Sign from PC1
      direction  = ifelse(PC1 > 0, "Positive", "Negative"),
      analysis   = pca_result$label,
      pc1_cor    = pca_result$pc_cors[1],
      pc2_cor    = pca_result$pc_cors[2]
    ) |>
    dplyr::arrange(desc(importance))
}

# =============================================================================
# COLLECT IMPORTANCE FROM ALL ANALYSES
# =============================================================================

cat("Extracting variable importance...\n\n")

# Total count
imp_total <- extract_importance(pca_total) |>
  dplyr::mutate(group = "Total", subgroup = "All Species")

# Species
imp_species <- do.call(rbind, lapply(names(pca_species), function(sp) {
  res <- extract_importance(pca_species[[sp]])
  if (is.null(res)) return(NULL)
  res |> dplyr::mutate(group = "Species", subgroup = sp)
}))

# Regions
imp_regions <- do.call(rbind, lapply(names(pca_regions), function(reg) {
  res <- extract_importance(pca_regions[[reg]])
  if (is.null(res)) return(NULL)
  res |> dplyr::mutate(group = "Region", subgroup = reg)
}))

all_importance <- rbind(imp_total, imp_species, imp_regions)

# =============================================================================
# RANK VARIABLES OVERALL
# =============================================================================

overall_rank <- all_importance |>
  dplyr::group_by(variable) |>
  dplyr::summarise(
    mean_importance     = mean(importance,       na.rm = TRUE),
    median_importance   = median(importance,     na.rm = TRUE),
    mean_loading        = mean(combined_loading, na.rm = TRUE),
    n_top3              = sum(rank(desc(importance)) <= 3,
                              na.rm = TRUE),
    pct_positive        = mean(direction == "Positive",
                               na.rm = TRUE) * 100,
    .groups = "drop"
  ) |>
  dplyr::arrange(desc(mean_importance)) |>
  dplyr::mutate(
    rank        = dplyr::row_number(),
    importance_class = dplyr::case_when(
      rank == 1             ~ "Most Important",
      rank <= 3             ~ "High",
      rank <= 5             ~ "Medium",
      TRUE                  ~ "Low"
    ),
    importance_class = factor(
      importance_class,
      levels = c("Most Important", "High", "Medium", "Low")
    )
  )

cat("=== Overall Variable Ranking ===\n")
print(overall_rank |>
        dplyr::select(rank, variable, mean_importance,
                      n_top3, pct_positive, importance_class),
      row.names = FALSE)

# Top variables
top_vars <- overall_rank |>
  dplyr::filter(rank <= 3) |>
  dplyr::pull(variable)

cat("\nTop 3 most important water variables:\n")
for (i in seq_along(top_vars)) {
  cat(sprintf("  %d. %s\n", i, top_vars[i]))
}

# =============================================================================
# PLOTS
# =============================================================================

out_dir <- "results/pca_plots"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -----------------------------------------------------------------------
# PLOT 1: Overall ranking bar chart
# -----------------------------------------------------------------------

p1 <- ggplot2::ggplot(
  overall_rank,
  ggplot2::aes(
    x    = forcats::fct_reorder(variable, mean_importance),
    y    = mean_importance,
    fill = importance_class
  )
) +
  ggplot2::geom_col(alpha = 0.85) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = mean_importance - median_importance * 0.1,
      ymax = mean_importance + median_importance * 0.1
    ),
    width = 0.3, colour = "grey40"
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = sprintf("#%d", rank),
      y     = mean_importance + 0.005
    ),
    hjust    = 0,
    size     = 3.5,
    fontface = "bold"
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "Most Important" = "#e74c3c",
      "High"           = "#e67e22",
      "Medium"         = "#f7dc6f",
      "Low"            = "#aed6f1"
    ),
    name = "Importance"
  ) +
  ggplot2::coord_flip() +
  ggplot2::labs(
    title    = "Water Variable Importance for Bird Counts",
    subtitle = "Averaged across total count, all species, and all regions",
    x        = NULL,
    y        = "Mean Weighted Importance Score"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    plot.title    = ggplot2::element_text(face = "bold", size = 13),
    plot.subtitle = ggplot2::element_text(size = 10, colour = "grey50"),
    legend.position = "right"
  )

ggplot2::ggsave(
  file.path(out_dir, "variable_importance_overall.png"),
  p1, width = 10, height = 6, dpi = 300, bg = "white"
)
cat("\nSaved: variable_importance_overall.png\n")

# -----------------------------------------------------------------------
# PLOT 2: Heatmap - importance by species
# -----------------------------------------------------------------------

imp_species_wide <- imp_species |>
  dplyr::select(variable, subgroup, importance) |>
  dplyr::mutate(
    subgroup = dplyr::recode(subgroup, !!!species_names)
  )

p2 <- ggplot2::ggplot(
  imp_species_wide,
  ggplot2::aes(
    x    = forcats::fct_reorder(variable, importance,
                                .fun = mean, .desc = TRUE),
    y    = subgroup,
    fill = importance
  )
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
  ggplot2::geom_text(
    ggplot2::aes(label = round(importance, 2)),
    size = 3, colour = "grey20"
  ) +
  ggplot2::scale_fill_gradient(
    low  = "#f7f9ff",
    high = "#e74c3c",
    name = "Importance"
  ) +
  ggplot2::labs(
    title    = "Water Variable Importance by Species",
    subtitle = "Darker = more important for that species",
    x        = "Water Variable",
    y        = NULL
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    plot.title      = ggplot2::element_text(face = "bold", size = 13),
    plot.subtitle   = ggplot2::element_text(size = 10, colour = "grey50"),
    axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

ggplot2::ggsave(
  file.path(out_dir, "variable_importance_by_species.png"),
  p2, width = 12, height = 6, dpi = 300, bg = "white"
)
cat("Saved: variable_importance_by_species.png\n")

# -----------------------------------------------------------------------
# PLOT 3: Heatmap - importance by region
# -----------------------------------------------------------------------

imp_regions_wide <- imp_regions |>
  dplyr::select(variable, subgroup, importance)

p3 <- ggplot2::ggplot(
  imp_regions_wide,
  ggplot2::aes(
    x    = forcats::fct_reorder(variable, importance,
                                .fun = mean, .desc = TRUE),
    y    = subgroup,
    fill = importance
  )
) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
  ggplot2::geom_text(
    ggplot2::aes(label = round(importance, 2)),
    size = 3, colour = "grey20"
  ) +
  ggplot2::scale_fill_gradient(
    low  = "#f7f9ff",
    high = "#2980b9",
    name = "Importance"
  ) +
  ggplot2::labs(
    title    = "Water Variable Importance by Region",
    subtitle = "Darker = more important for that region",
    x        = "Water Variable",
    y        = "Region"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    plot.title      = ggplot2::element_text(face = "bold", size = 13),
    plot.subtitle   = ggplot2::element_text(size = 10, colour = "grey50"),
    axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

ggplot2::ggsave(
  file.path(out_dir, "variable_importance_by_region.png"),
  p3, width = 12, height = 6, dpi = 300, bg = "white"
)
cat("Saved: variable_importance_by_region.png\n")

# -----------------------------------------------------------------------
# PLOT 4: Dot plot - top 3 variables highlighted across all analyses
# -----------------------------------------------------------------------

top3_highlight <- all_importance |>
  dplyr::mutate(
    is_top3    = variable %in% top_vars,
    var_label  = ifelse(is_top3, variable, "Other"),
    group_label = dplyr::case_when(
      group == "Species" ~
        dplyr::recode(subgroup, !!!species_names),
      TRUE ~ subgroup
    )
  )

p4 <- ggplot2::ggplot(
  top3_highlight,
  ggplot2::aes(
    x     = forcats::fct_reorder(group_label,
                                 importance, .fun = max),
    y     = importance,
    colour = variable,
    size  = is_top3,
    alpha = is_top3
  )
) +
  ggplot2::geom_point() +
  ggplot2::scale_size_manual(
    values = c("TRUE" = 4, "FALSE" = 1.5),
    guide  = "none"
  ) +
  ggplot2::scale_alpha_manual(
    values = c("TRUE" = 1.0, "FALSE" = 0.3),
    guide  = "none"
  ) +
  ggplot2::scale_colour_brewer(
    palette = "Set1",
    name    = "Variable"
  ) +
  ggplot2::coord_flip() +
  ggplot2::facet_wrap(~group, scales = "free_y", ncol = 1) +
  ggplot2::labs(
    title    = "Variable Importance Across All Analyses",
    subtitle = paste("Top 3 highlighted:",
                     paste(top_vars, collapse = ", ")),
    x        = NULL,
    y        = "Importance Score"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    plot.title      = ggplot2::element_text(face = "bold", size = 13),
    plot.subtitle   = ggplot2::element_text(size = 10, colour = "grey50"),
    strip.text      = ggplot2::element_text(face = "bold"),
    legend.position = "right"
  )

ggplot2::ggsave(
  file.path(out_dir, "variable_importance_dotplot.png"),
  p4, width = 12, height = 14, dpi = 300, bg = "white"
)
cat("Saved: variable_importance_dotplot.png\n")

# -----------------------------------------------------------------------
# PLOT 5: Direction - does more = more birds or less?
# -----------------------------------------------------------------------

direction_summary <- all_importance |>
  dplyr::filter(variable %in% top_vars) |>
  dplyr::group_by(variable, group) |>
  dplyr::summarise(
    pct_positive = mean(direction == "Positive") * 100,
    mean_pc1_loading = mean(PC1),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    interpretation = dplyr::case_when(
      pct_positive >= 75 ~ "More = More Birds",
      pct_positive <= 25 ~ "More = Fewer Birds",
      TRUE               ~ "Mixed / Nonlinear"
    )
  )

p5 <- ggplot2::ggplot(
  direction_summary,
  ggplot2::aes(
    x    = variable,
    y    = mean_pc1_loading,
    fill = interpretation
  )
) +
  ggplot2::geom_col(alpha = 0.85) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey40") +
  ggplot2::facet_wrap(~group, ncol = 3) +
  ggplot2::scale_fill_manual(
    values = c(
      "More = More Birds"   = "#2ecc71",
      "More = Fewer Birds"  = "#e74c3c",
      "Mixed / Nonlinear"   = "#f7dc6f"
    ),
    name = NULL
  ) +
  ggplot2::labs(
    title    = paste("Direction of Effect - Top Variables:",
                     paste(top_vars, collapse = ", ")),
    subtitle = "Does more of this variable = more or fewer birds?",
    x        = "Variable",
    y        = "Mean PC1 Loading"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    plot.title      = ggplot2::element_text(face = "bold", size = 12),
    plot.subtitle   = ggplot2::element_text(size = 10, colour = "grey50"),
    strip.text      = ggplot2::element_text(face = "bold"),
    legend.position = "bottom",
    axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1)
  )

ggplot2::ggsave(
  file.path(out_dir, "variable_direction_top3.png"),
  p5, width = 12, height = 8, dpi = 300, bg = "white"
)
cat("Saved: variable_direction_top3.png\n")

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("\n")
cat("============================================================\n")
cat("VARIABLE SELECTION SUMMARY\n")
cat("============================================================\n\n")

cat("Most important water variables (ranked):\n\n")
for (i in 1:nrow(overall_rank)) {
  row <- overall_rank[i, ]
  cat(sprintf(
    "%d. %-20s | Importance: %.3f | Top-3 appearances: %d | Direction: %s%%+\n",
    row$rank,
    row$variable,
    row$mean_importance,
    row$n_top3,
    round(row$pct_positive)
  ))
}

cat("\n------------------------------------------------------------\n")
cat("RECOMMENDATION: Use these variables in forecasting models:\n\n")
for (i in 1:min(3, nrow(overall_rank))) {
  row <- overall_rank[i, ]
  dir_text <- ifelse(row$pct_positive >= 75,
                     "positive effect on birds",
                     ifelse(row$pct_positive <= 25,
                            "negative effect on birds",
                            "mixed/nonlinear effect"))
  cat(sprintf("  %d. %s (%s)\n", i, row$variable, dir_text))
}
cat("------------------------------------------------------------\n\n")

cat("Plots saved to:", out_dir, "\n")
print(list.files(out_dir,
                 pattern = "variable_importance|variable_direction"))
