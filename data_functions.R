# piper_ecology.R
# ============================================================================
# Piper-style ternary diagrams for Everglades wading bird ecology
# Can be run independently of main.R
# Each point = one region x one year
# Each species has independent color scale
# ============================================================================

# =============================================================================
# DEPENDENCIES
# =============================================================================

`%||%` <- function(x, y) if (is.null(x)) y else x

required_packages <- c("ggtern", "grid", "gridExtra",
                       "dplyr", "tidyr", "ggplot2",
                       "scales", "tsibble", "readr",
                       "viridis")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# Resolve ggtern/ggplot2 conflicts
aes                      <- ggplot2::aes
theme_bw                 <- ggplot2::theme_bw
theme                    <- ggplot2::theme
labs                     <- ggplot2::labs
ggsave                   <- ggplot2::ggsave
geom_point               <- ggplot2::geom_point
geom_path                <- ggplot2::geom_path
geom_text                <- ggplot2::geom_text
facet_wrap               <- ggplot2::facet_wrap
scale_color_gradient     <- ggplot2::scale_color_gradient
scale_color_gradient2    <- ggplot2::scale_color_gradient2
scale_color_viridis_c    <- ggplot2::scale_color_viridis_c
scale_size_continuous    <- ggplot2::scale_size_continuous
scale_fill_manual        <- ggplot2::scale_fill_manual
element_text             <- ggplot2::element_text
annotate                 <- ggplot2::annotate
unit                     <- ggplot2::unit

# Resolve dplyr conflicts
filter    <- dplyr::filter
select    <- dplyr::select
arrange   <- dplyr::arrange
mutate    <- dplyr::mutate
summarise <- dplyr::summarise
group_by  <- dplyr::group_by
left_join <- dplyr::left_join
slice_min <- dplyr::slice_min
slice_max <- dplyr::slice_max

# =============================================================================
# LOAD PROJECT FUNCTIONS
# =============================================================================

if (!exists("get_wading_bird_data")) {
  source("data_functions.R")
}

# =============================================================================
# CONFIGURATION
# =============================================================================

spatial_level  <- "subregion"
species_to_use <- c("gbhe", "greg", "rosp", "sneg", "whib", "wost")
year_min       <- 1993
year_max       <- 2025

cat("=== Piper Ecology Setup ===\n")
cat("Spatial level:", spatial_level, "\n")
cat("Species:", paste(species_to_use, collapse = ", "), "\n")
cat("Years:", year_min, "-", year_max, "\n\n")

# =============================================================================
# LOAD DATA
# =============================================================================

cat("Loading data...\n")

piper_config <- list(
  spatial = list(
    level              = spatial_level,
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
  cat("Downloading wader observation data...\n")
  library(wader)
  download_observations(".")
}

data_sub <- get_wading_bird_data(piper_config)

cat("Data loaded:", nrow(data_sub), "observations\n")
cat("Regions:", paste(unique(as_tibble(data_sub)$region),
                      collapse = ", "), "\n\n")

# =============================================================================
# PREPARE TERNARY DATA
# Each row = one region x one year (NOT averaged)
# =============================================================================

prepare_ternary_data <- function(data,
                                 species_filter = NULL,
                                 year_min       = NULL,
                                 year_max       = NULL) {
  
  df <- data |> as_tibble()
  
  if (!is.null(species_filter)) {
    df <- df |> dplyr::filter(species %in% species_filter)
  }
  if (!is.null(year_min)) df <- df |> dplyr::filter(year >= year_min)
  if (!is.null(year_max)) df <- df |> dplyr::filter(year <= year_max)
  
  region_year <- df |>
    dplyr::group_by(year, region) |>
    dplyr::summarise(
      total_count        = sum(count,                na.rm = TRUE),
      breed_season_depth = mean(breed_season_depth,  na.rm = TRUE),
      dry_days           = mean(dry_days,            na.rm = TRUE),
      recession          = mean(recession,           na.rm = TRUE),
      init_depth         = mean(init_depth,          na.rm = TRUE),
      reversals          = mean(reversals,           na.rm = TRUE),
      .groups            = "drop"
    )
  
  # ==========================================================
  # TRIANGLE 1: WATER REGIME
  # Dry     = high dry_days
  # Wet     = high breed_season_depth + init_depth
  # Dynamic = high recession + reversals
  # ==========================================================
  
  region_year <- region_year |>
    dplyr::mutate(
      raw_dry     = dry_days,
      raw_wet     = breed_season_depth + init_depth,
      raw_dynamic = abs(recession) * 100 + reversals,
      total_water = pmax(raw_dry + raw_wet + raw_dynamic, 0.001),
      T1_dry      = raw_dry     / total_water,
      T1_wet      = raw_wet     / total_water,
      T1_dynamic  = raw_dynamic / total_water
    )
  
  # ==========================================================
  # TRIANGLE 2: FORAGING OPPORTUNITY
  # Concentration = dry_days x recession (prey trapped)
  # Extent        = breed_depth x init_depth (habitat volume)
  # Access        = reversals + recession (rate of change)
  # ==========================================================
  
  region_year <- region_year |>
    dplyr::mutate(
      recession_pos     = pmax(recession, 0),
      raw_concentration = dry_days * recession_pos * 100,
      raw_extent        = breed_season_depth * init_depth,
      raw_access        = reversals + abs(recession) * 50,
      total_forage      = pmax(
        raw_concentration + raw_extent + raw_access,
        0.001
      ),
      T2_concentration  = raw_concentration / total_forage,
      T2_extent         = raw_extent        / total_forage,
      T2_access         = raw_access        / total_forage
    )
  
  cat("Ternary data prepared:", nrow(region_year), "rows\n")
  cat("  Each row = 1 region x 1 year (no averaging)\n")
  cat("Regions:", paste(sort(unique(region_year$region)),
                        collapse = ", "), "\n")
  cat("Years:", min(region_year$year), "-",
      max(region_year$year), "\n\n")
  
  return(region_year)
}

# =============================================================================
# PLOT 1: SPATIAL - all region x year combinations
# =============================================================================

plot_spatial_ternary <- function(ternary_data,
                                 triangle   = 1,
                                 facet_by   = NULL,
                                 year_range = NULL,
                                 title      = NULL) {
  
  if (!is.null(year_range)) {
    ternary_data <- ternary_data |>
      dplyr::filter(year >= year_range[1],
                    year <= year_range[2])
  }
  
  midpoint <- median(ternary_data$total_count, na.rm = TRUE)
  
  if (triangle == 1) {
    p <- ggtern(
      data = ternary_data,
      aes(x = T1_dry, y = T1_wet, z = T1_dynamic,
          color = total_count, size = total_count)
    ) +
      Tlab("Dynamic\n(Recession+Reversals)") +
      Llab("Dry\n(Dry Days)") +
      Rlab("Wet\n(Water Depth)") +
      labs(
        title = title %||%
          "Water Regime: Each point = 1 region x 1 year",
        color = "Bird Count", size = "Bird Count"
      )
  } else {
    p <- ggtern(
      data = ternary_data,
      aes(x = T2_concentration, y = T2_extent, z = T2_access,
          color = total_count, size = total_count)
    ) +
      Tlab("Access\n(Recession+Reversals)") +
      Llab("Concentration\n(Dry x Drawdown)") +
      Rlab("Extent\n(Depth x Volume)") +
      labs(
        title = title %||%
          "Foraging Opportunity: Each point = 1 region x 1 year",
        color = "Bird Count", size = "Bird Count"
      )
  }
  
  p <- p +
    geom_point(alpha = 0.7) +
    scale_color_gradient2(
      low = "#3498db", mid = "#f39c12", high = "#e74c3c",
      midpoint = midpoint, labels = scales::comma
    ) +
    scale_size_continuous(range = c(1, 8),
                          labels = scales::comma) +
    theme_bw() +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold", size = 11))
  
  if (!is.null(facet_by)) {
    p <- p + facet_wrap(as.formula(paste("~", facet_by)),
                        ncol = 4)
  }
  
  return(p)
}

# =============================================================================
# PLOT 2: BY SPECIES - independent color scale per species
# =============================================================================

plot_species_ternary <- function(ternary_data,
                                 data,
                                 triangle = 1,
                                 title    = NULL) {
  
  species_counts <- data |>
    as_tibble() |>
    dplyr::group_by(year, region, species) |>
    dplyr::summarise(count = sum(count, na.rm = TRUE),
                     .groups = "drop")
  
  plot_data <- ternary_data |>
    dplyr::left_join(species_counts, by = c("year", "region"))
  
  species_list <- sort(unique(plot_data$species))
  plot_list    <- list()
  
  for (sp in species_list) {
    
    sp_data   <- plot_data |> dplyr::filter(species == sp)
    sp_median <- median(sp_data$count, na.rm = TRUE)
    sp_max    <- max(sp_data$count,    na.rm = TRUE)
    
    if (triangle == 1) {
      p <- ggtern(
        data = sp_data,
        aes(x = T1_dry, y = T1_wet, z = T1_dynamic,
            color = count, size = count)
      ) +
        Tlab("Dynamic") +
        Llab("Dry") +
        Rlab("Wet")
    } else {
      p <- ggtern(
        data = sp_data,
        aes(x = T2_concentration, y = T2_extent, z = T2_access,
            color = count, size = count)
      ) +
        Tlab("Access") +
        Llab("Concentration") +
        Rlab("Extent")
    }
    
    p <- p +
      geom_point(alpha = 0.7) +
      scale_color_gradient2(
        low      = "#3498db",
        mid      = "#f39c12",
        high     = "#e74c3c",
        midpoint = sp_median,
        limits   = c(0, sp_max),
        labels   = scales::comma,
        name     = "Count"
      ) +
      scale_size_continuous(
        range  = c(1, 6),
        limits = c(0, sp_max),
        labels = scales::comma,
        name   = "Count"
      ) +
      labs(title = toupper(sp)) +
      theme_bw() +
      theme(
        legend.position  = "bottom",
        legend.key.width = ggplot2::unit(1.5, "cm"),
        plot.title       = element_text(face  = "bold",
                                        size  = 10,
                                        hjust = 0.5),
        legend.text      = element_text(size = 7),
        legend.title     = element_text(size = 8)
      )
    
    plot_list[[sp]] <- p
  }
  
  main_title <- grid::textGrob(
    label = if (!is.null(title)) title else {
      if (triangle == 1) {
        "Water Regime by Species\n(independent color scales per species)"
      } else {
        "Foraging Opportunity by Species\n(independent color scales per species)"
      }
    },
    gp = grid::gpar(fontsize = 13, fontface = "bold")
  )
  
  gridExtra::grid.arrange(
    grobs = plot_list,
    ncol  = 3,
    top   = main_title
  )
}

# =============================================================================
# PLOT 3: INDIVIDUAL SPECIES - one at a time
# =============================================================================

plot_one_species_ternary <- function(ternary_data,
                                     data,
                                     species_name,
                                     triangle = 1,
                                     title    = NULL) {
  
  sp_counts <- data |>
    as_tibble() |>
    dplyr::filter(species == species_name) |>
    dplyr::group_by(year, region) |>
    dplyr::summarise(count = sum(count, na.rm = TRUE),
                     .groups = "drop")
  
  plot_data <- ternary_data |>
    dplyr::left_join(sp_counts, by = c("year", "region"))
  
  sp_median <- median(plot_data$count, na.rm = TRUE)
  sp_max    <- max(plot_data$count,    na.rm = TRUE)
  
  cat(sprintf("Species: %s | Max: %.0f | Median: %.0f\n",
              species_name, sp_max, sp_median))
  
  if (triangle == 1) {
    p <- ggtern(
      data = plot_data,
      aes(x = T1_dry, y = T1_wet, z = T1_dynamic,
          color = count, size = count)
    ) +
      Tlab("Dynamic\n(Recession+Reversals)") +
      Llab("Dry\n(Dry Days)") +
      Rlab("Wet\n(Water Depth)") +
      labs(
        title = title %||%
          sprintf("Water Regime - %s\n(each point = 1 region x 1 year)",
                  toupper(species_name)),
        color = "Count", size = "Count"
      )
  } else {
    p <- ggtern(
      data = plot_data,
      aes(x = T2_concentration, y = T2_extent, z = T2_access,
          color = count, size = count)
    ) +
      Tlab("Access\n(Recession+Reversals)") +
      Llab("Concentration\n(Dry x Drawdown)") +
      Rlab("Extent\n(Depth x Volume)") +
      labs(
        title = title %||%
          sprintf("Foraging - %s\n(each point = 1 region x 1 year)",
                  toupper(species_name)),
        color = "Count", size = "Count"
      )
  }
  
  p <- p +
    geom_point(alpha = 0.7) +
    scale_color_gradient2(
      low      = "#3498db",
      mid      = "#f39c12",
      high     = "#e74c3c",
      midpoint = sp_median,
      limits   = c(0, sp_max),
      labels   = scales::comma
    ) +
    scale_size_continuous(
      range  = c(1, 8),
      limits = c(0, sp_max),
      labels = scales::comma
    ) +
    theme_bw() +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold", size = 12))
  
  return(p)
}

# =============================================================================
# PLOT 4: HIGH COUNT YEARS BY SPECIES - independent threshold per species
# =============================================================================

plot_high_count_by_species <- function(ternary_data,
                                       data,
                                       triangle = 1,
                                       top_pct  = 0.25,
                                       title    = NULL) {
  
  species_counts <- data |>
    as_tibble() |>
    dplyr::group_by(year, region, species) |>
    dplyr::summarise(count = sum(count, na.rm = TRUE),
                     .groups = "drop")
  
  plot_data <- ternary_data |>
    dplyr::left_join(species_counts, by = c("year", "region"))
  
  # Independent threshold per species
  plot_data <- plot_data |>
    dplyr::group_by(species) |>
    dplyr::mutate(
      sp_threshold = quantile(count, 1 - top_pct, na.rm = TRUE),
      count_class  = dplyr::if_else(
        count >= sp_threshold, "High Count", "Normal"
      ),
      count_class  = factor(count_class,
                            levels = c("Normal", "High Count"))
    ) |>
    dplyr::ungroup()
  
  species_list <- sort(unique(plot_data$species))
  plot_list    <- list()
  
  for (sp in species_list) {
    
    sp_data      <- plot_data |> dplyr::filter(species == sp)
    sp_threshold <- round(unique(sp_data$sp_threshold), 0)
    
    if (triangle == 1) {
      p <- ggtern(
        data = sp_data,
        aes(x = T1_dry, y = T1_wet, z = T1_dynamic,
            color = count_class, size = count,
            alpha = count_class)
      ) +
        Tlab("Dynamic") +
        Llab("Dry") +
        Rlab("Wet")
    } else {
      p <- ggtern(
        data = sp_data,
        aes(x = T2_concentration, y = T2_extent, z = T2_access,
            color = count_class, size = count,
            alpha = count_class)
      ) +
        Tlab("Access") +
        Llab("Concentration") +
        Rlab("Extent")
    }
    
    p <- p +
      geom_point() +
      scale_color_manual(
        values = c("Normal"     = "grey70",
                   "High Count" = "#e74c3c"),
        name   = NULL
      ) +
      ggplot2::scale_alpha_manual(
        values = c("Normal"     = 0.35,
                   "High Count" = 0.9),
        guide  = "none"
      ) +
      scale_size_continuous(
        range = c(1, 6),
        guide = "none"
      ) +
      labs(
        title = sprintf("%s\n(top %.0f%% >= %.0f)",
                        toupper(sp),
                        top_pct * 100,
                        sp_threshold)
      ) +
      theme_bw() +
      theme(
        legend.position = "bottom",
        plot.title      = element_text(face  = "bold",
                                       size  = 9,
                                       hjust = 0.5),
        legend.text     = element_text(size = 8)
      )
    
    plot_list[[sp]] <- p
  }
  
  main_title <- grid::textGrob(
    label = if (!is.null(title)) title else {
      sprintf(
        "High Count Years by Species (top %.0f%%, independent thresholds)\n%s",
        top_pct * 100,
        if (triangle == 1) "Water Regime" else "Foraging Opportunity"
      )
    },
    gp = grid::gpar(fontsize = 12, fontface = "bold")
  )
  
  gridExtra::grid.arrange(
    grobs = plot_list,
    ncol  = 3,
    top   = main_title
  )
}

# =============================================================================
# PLOT 5: TEMPORAL TRAJECTORIES
# =============================================================================

plot_temporal_ternary <- function(ternary_data,
                                  triangle = 1,
                                  regions  = NULL,
                                  title    = NULL) {
  
  if (!is.null(regions)) {
    ternary_data <- ternary_data |>
      dplyr::filter(region %in% regions)
  }
  
  ternary_data <- ternary_data |>
    dplyr::arrange(region, year)
  
  midpoint <- median(ternary_data$total_count, na.rm = TRUE)
  
  if (triangle == 1) {
    p <- ggtern(
      data = ternary_data,
      aes(x = T1_dry, y = T1_wet, z = T1_dynamic,
          color = total_count, group = region)
    ) +
      Tlab("Dynamic") +
      Llab("Dry") +
      Rlab("Wet") +
      labs(title = title %||% "Water Regime Trajectories Over Time",
           color = "Bird Count")
  } else {
    p <- ggtern(
      data = ternary_data,
      aes(x = T2_concentration, y = T2_extent, z = T2_access,
          color = total_count, group = region)
    ) +
      Tlab("Access") +
      Llab("Concentration") +
      Rlab("Extent") +
      labs(title = title %||%
             "Foraging Opportunity Trajectories Over Time",
           color = "Bird Count")
  }
  
  p <- p +
    geom_path(linewidth = 0.5, alpha = 0.6) +
    geom_point(size = 2, alpha = 0.8) +
    geom_point(
      data   = ternary_data |>
        dplyr::group_by(region) |>
        dplyr::slice_min(year, n = 1),
      shape  = 21, size = 4,
      fill   = "white", stroke = 1.5
    ) +
    geom_point(
      data  = ternary_data |>
        dplyr::group_by(region) |>
        dplyr::slice_max(year, n = 1),
      shape = 24, size = 4, fill = "black"
    ) +
    facet_wrap(~region, ncol = 4) +
    scale_color_gradient2(
      low      = "#3498db",
      mid      = "#f39c12",
      high     = "#e74c3c",
      midpoint = midpoint,
      labels   = scales::comma
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      strip.text      = element_text(face = "bold"),
      plot.title      = element_text(face = "bold", size = 11)
    )
  
  return(p)
}

# =============================================================================
# PLOT 6: PIPER STYLE - Both triangles side by side
# =============================================================================

plot_piper_style <- function(ternary_data,
                             year_filter = NULL,
                             title       = NULL) {
  
  if (!is.null(year_filter)) {
    ternary_data <- ternary_data |>
      dplyr::filter(year %in% year_filter)
  }
  
  midpoint <- median(ternary_data$total_count, na.rm = TRUE)
  
  p1 <- ggtern(
    data = ternary_data,
    aes(x = T1_dry, y = T1_wet, z = T1_dynamic,
        color = total_count, size = total_count)
  ) +
    geom_point(alpha = 0.7) +
    Tlab("Dynamic") + Llab("Dry") + Rlab("Wet") +
    scale_color_gradient2(
      low = "#3498db", mid = "#f39c12", high = "#e74c3c",
      midpoint = midpoint, labels = scales::comma
    ) +
    scale_size_continuous(range = c(1, 8),
                          labels = scales::comma) +
    labs(title = "Water Regime",
         color = "Count", size = "Count") +
    theme_bw() +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold"))
  
  p2 <- ggtern(
    data = ternary_data,
    aes(x = T2_concentration, y = T2_extent, z = T2_access,
        color = total_count, size = total_count)
  ) +
    geom_point(alpha = 0.7) +
    Tlab("Access") + Llab("Concentration") + Rlab("Extent") +
    scale_color_gradient2(
      low = "#3498db", mid = "#f39c12", high = "#e74c3c",
      midpoint = midpoint, labels = scales::comma
    ) +
    scale_size_continuous(range = c(1, 8),
                          labels = scales::comma) +
    labs(title = "Foraging Opportunity",
         color = "Count", size = "Count") +
    theme_bw() +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold"))
  
  gridExtra::grid.arrange(
    p1, p2,
    ncol = 2,
    top  = grid::textGrob(
      title %||%
        "Piper-Style Ecological Classification\n(each point = 1 region x 1 year)",
      gp = grid::gpar(fontsize = 13, fontface = "bold")
    )
  )
}

# =============================================================================
# RUN ANALYSIS
# =============================================================================

cat("=== Preparing ternary data ===\n")

tern_data <- prepare_ternary_data(
  data           = data_sub,
  species_filter = species_to_use,
  year_min       = year_min,
  year_max       = year_max
)

out_dir <- "results/piper_plots"
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
  cat("Created:", out_dir, "\n")
}

cat("=== Generating plots ===\n\n")

# Plot 1: Water regime all years
cat("Plot 1: Water regime (all years)...\n")
p1 <- plot_spatial_ternary(tern_data, triangle = 1)
ggplot2::ggsave(file.path(out_dir, "water_regime_all_years.png"),
                p1, width = 8, height = 7, dpi = 300, bg = "white")

# Plot 2: Foraging all years
cat("Plot 2: Foraging (all years)...\n")
p2 <- plot_spatial_ternary(tern_data, triangle = 2)
ggplot2::ggsave(file.path(out_dir, "foraging_all_years.png"),
                p2, width = 8, height = 7, dpi = 300, bg = "white")

# Plot 3: Water by region
cat("Plot 3: Water regime by region...\n")
p3 <- plot_spatial_ternary(tern_data, triangle = 1,
                           facet_by = "region")
ggplot2::ggsave(file.path(out_dir, "water_regime_by_region.png"),
                p3, width = 16, height = 12, dpi = 300, bg = "white")

# Plot 4: Foraging by region
cat("Plot 4: Foraging by region...\n")
p4 <- plot_spatial_ternary(tern_data, triangle = 2,
                           facet_by = "region")
ggplot2::ggsave(file.path(out_dir, "foraging_by_region.png"),
                p4, width = 16, height = 12, dpi = 300, bg = "white")

# Plot 5: Piper combined
cat("Plot 5: Piper combined...\n")
png(file.path(out_dir, "piper_combined.png"),
    width = 14, height = 7, units = "in", res = 300, bg = "white")
plot_piper_style(tern_data)
dev.off()

# Plot 6: Water trajectories
cat("Plot 6: Water regime trajectories...\n")
p6 <- plot_temporal_ternary(tern_data, triangle = 1)
ggplot2::ggsave(file.path(out_dir, "water_regime_trajectories.png"),
                p6, width = 14, height = 10, dpi = 300, bg = "white")

# Plot 7: Foraging trajectories
cat("Plot 7: Foraging trajectories...\n")
p7 <- plot_temporal_ternary(tern_data, triangle = 2)
ggplot2::ggsave(file.path(out_dir, "foraging_trajectories.png"),
                p7, width = 14, height = 10, dpi = 300, bg = "white")

# Plot 8: Species water (independent scales)
cat("Plot 8: Species water regime (independent scales)...\n")
png(file.path(out_dir, "species_water_regime.png"),
    width = 16, height = 12, units = "in", res = 300, bg = "white")
plot_species_ternary(tern_data, data_sub, triangle = 1)
dev.off()

# Plot 9: Species foraging (independent scales)
cat("Plot 9: Species foraging (independent scales)...\n")
png(file.path(out_dir, "species_foraging.png"),
    width = 16, height = 12, units = "in", res = 300, bg = "white")
plot_species_ternary(tern_data, data_sub, triangle = 2)
dev.off()

# Plot 10: High count by species water
cat("Plot 10: High count by species (water)...\n")
png(file.path(out_dir, "high_count_species_water.png"),
    width = 16, height = 12, units = "in", res = 300, bg = "white")
plot_high_count_by_species(tern_data, data_sub, triangle = 1)
dev.off()

# Plot 11: High count by species foraging
cat("Plot 11: High count by species (foraging)...\n")
png(file.path(out_dir, "high_count_species_foraging.png"),
    width = 16, height = 12, units = "in", res = 300, bg = "white")
plot_high_count_by_species(tern_data, data_sub, triangle = 2)
dev.off()

# Individual species plots
cat("Individual species plots...\n")
for (sp in species_to_use) {
  cat(" ", sp, "\n")
  
  p <- plot_one_species_ternary(tern_data, data_sub,
                                species_name = sp, triangle = 1)
  ggplot2::ggsave(
    file.path(out_dir, sprintf("individual_%s_water.png", sp)),
    p, width = 7, height = 6, dpi = 300, bg = "white"
  )
  
  p <- plot_one_species_ternary(tern_data, data_sub,
                                species_name = sp, triangle = 2)
  ggplot2::ggsave(
    file.path(out_dir, sprintf("individual_%s_foraging.png", sp)),
    p, width = 7, height = 6, dpi = 300, bg = "white"
  )
}

# =============================================================================
# SUMMARY TABLE
# =============================================================================

cat("\n=== Region Classification Summary ===\n")

summary_table <- tern_data |>
  dplyr::group_by(region) |>
  dplyr::summarise(
    n_years         = dplyr::n(),
    mean_count      = round(mean(total_count, na.rm = TRUE), 0),
    max_count       = round(max(total_count,  na.rm = TRUE), 0),
    best_year       = year[which.max(total_count)],
    mean_dry        = round(mean(T1_dry,           na.rm = TRUE), 3),
    mean_wet        = round(mean(T1_wet,           na.rm = TRUE), 3),
    mean_dynamic    = round(mean(T1_dynamic,       na.rm = TRUE), 3),
    mean_conc       = round(mean(T2_concentration, na.rm = TRUE), 3),
    mean_extent     = round(mean(T2_extent,        na.rm = TRUE), 3),
    mean_access     = round(mean(T2_access,        na.rm = TRUE), 3),
    dominant_regime = dplyr::case_when(
      mean_dry  == pmax(mean_dry, mean_wet, mean_dynamic) ~ "Dry",
      mean_wet  == pmax(mean_dry, mean_wet, mean_dynamic) ~ "Wet",
      TRUE                                                  ~ "Dynamic"
    ),
    dominant_forage = dplyr::case_when(
      mean_conc   == pmax(mean_conc, mean_extent, mean_access) ~
        "Concentration",
      mean_extent == pmax(mean_conc, mean_extent, mean_access) ~
        "Extent",
      TRUE ~ "Access"
    ),
    .groups = "drop"
  ) |>
  dplyr::arrange(desc(mean_count))

print(summary_table)

cat("\n=== Complete ===\n")
cat("Plots saved to:", out_dir, "\n\n")
cat("Files generated:\n")
print(list.files(out_dir))