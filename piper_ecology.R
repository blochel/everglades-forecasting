# piper_ecology_ggplot.R
# ============================================================================
# Ecological Piper diagram using ggplot2
# Properly aligned diamond touching both triangles
#
# LEFT TRIANGLE  (Water Regime):   Dry / Wet / Dynamic
# RIGHT TRIANGLE (Foraging):       Concentration / Extent / Access
# CENTRAL DIAMOND:                 Combined ecological type
#
#        Dynamic+Access (top)
#             /\
#            /  \
#           /    \
#  Dry+    /      \ Wet+
#  Conc   /        \ Extent
#        /          \
#       /  DIAMOND   \
#      /              \
#     /________________\
#    /\    bottom      /\
#   /  \              /  \
#  / L  \            /  R \
# /  TRI \          /  TRI \
#/________\________/________\
# Dry    Wet    Conc       Extent
# ============================================================================

# =============================================================================
# PACKAGES
# =============================================================================

required_packages <- c("ggplot2", "dplyr", "scales",
                       "gridExtra", "grid",
                       "tsibble", "readr")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# Explicitly assign to avoid any masking
ggplot                 <- ggplot2::ggplot
aes                    <- ggplot2::aes
geom_point             <- ggplot2::geom_point
geom_text              <- ggplot2::geom_text
geom_polygon           <- ggplot2::geom_polygon
geom_line              <- ggplot2::geom_line
scale_colour_gradient2 <- ggplot2::scale_colour_gradient2
scale_size_identity    <- ggplot2::scale_size_identity
coord_fixed            <- ggplot2::coord_fixed
theme_void             <- ggplot2::theme_void
theme                  <- ggplot2::theme
element_text           <- ggplot2::element_text
margin                 <- ggplot2::margin
labs                   <- ggplot2::labs
ggsave                 <- ggplot2::ggsave

filter    <- dplyr::filter
mutate    <- dplyr::mutate
summarise <- dplyr::summarise
left_join <- dplyr::left_join
rename    <- dplyr::rename

# =============================================================================
# LOAD DATA
# =============================================================================

if (!exists("tern_data")) {
  
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
  
  if (!exists("data_sub")) {
    if (!dir.exists("SiteandMethods")) {
      library(wader)
      download_observations(".")
    }
    data_sub <- get_wading_bird_data(piper_config)
  }
  
  if (!exists("prepare_ternary_data")) source("piper_ecology.R")
}

# =============================================================================
# GEOMETRY
# All coordinates derived from W, GAP, H
# Key alignment: diamond vertices = triangle top corners + midpoint
# =============================================================================

H   <- sqrt(3) / 2
W   <- 1.0
GAP <- 0.2

# Key alignment points
L_BL  <- c(0,               0)     # left  tri bottom-left
L_BR  <- c(W,               0)     # left  tri bottom-right  = diamond left-bot
L_TOP <- c(W / 2,           H)     # left  tri top           = diamond left
R_BL  <- c(W + GAP,         0)     # right tri bottom-left   = diamond right-bot
R_BR  <- c(2 * W + GAP,     0)     # right tri bottom-right
R_TOP <- c(W + GAP + W / 2, H)     # right tri top           = diamond right
D_BOT <- c(W + GAP / 2,     0)     # diamond bottom (between inner corners)
D_TOP <- c(W + GAP / 2,     H * 2) # diamond top (projected above)

# =============================================================================
# COORDINATE TRANSFORMS
# =============================================================================

left_xy <- function(dry, wet, dynamic) {
  total <- dry + wet + dynamic
  b     <- wet     / total
  c     <- dynamic / total
  # x: starts at L_BL[1]=0, moves right with wet, center with dynamic
  x <- L_BL[1] + b * W + c * W / 2
  y <- c * H
  data.frame(x = x, y = y)
}

right_xy <- function(conc, extent, access) {
  total <- conc + extent + access
  b     <- extent / total
  c     <- access / total
  # x: starts at R_BL[1]
  x <- R_BL[1] + b * W + c * W / 2
  y <- c * H
  data.frame(x = x, y = y)
}

diamond_xy <- function(dry, wet, dynamic, conc, extent, access) {
  
  tw    <- dry  + wet     + dynamic
  tf    <- conc + extent  + access
  
  w_dry <- dry     / tw
  w_wet <- wet     / tw
  w_dyn <- dynamic / tw
  f_con <- conc    / tf
  f_ext <- extent  / tf
  f_acc <- access  / tf
  
  # Diamond is defined by 4 vertices:
  # L_TOP (left), R_TOP (right), D_BOT (bottom), D_TOP (top)
  #
  # A point in the diamond is parameterised by:
  #   s = position along left-right axis  (0=left, 1=right)
  #   t = position along bottom-top axis  (0=bottom, 1=top)
  #
  # s driven by: w_wet (left tri wet axis) + f_ext (right tri extent axis)
  # t driven by: w_dyn (left tri dynamic)  + f_acc (right tri access axis)
  
  s <- (w_wet + f_ext) / 2   # 0 to 1
  t <- (w_dyn + f_acc) / 2   # 0 to 1
  
  # Bilinear interpolation within diamond rhombus
  # Bottom half: interpolate between D_BOT and (L_TOP or R_TOP)
  # Top    half: interpolate between (L_TOP or R_TOP) and D_TOP
  
  # x: linear between L_TOP and R_TOP based on s
  x_mid <- L_TOP[1] + s * (R_TOP[1] - L_TOP[1])
  
  # y: linear between D_BOT and D_TOP based on t
  # but x shifts as we go up
  x <- D_BOT[1] + s * (R_TOP[1] - L_TOP[1]) * (1 - abs(t - 0.5) * 0.5)
  y <- D_BOT[2] + t * (D_TOP[2] - D_BOT[2])
  
  # Clamp within diamond bounds
  x <- pmax(L_TOP[1] + 0.02,
            pmin(R_TOP[1] - 0.02, x))
  y <- pmax(D_BOT[2] + 0.02,
            pmin(D_TOP[2] - 0.02, y))
  
  data.frame(x = x, y = y)
}

# =============================================================================
# STRUCTURE BUILDERS
# =============================================================================

make_triangle <- function(bl_x) {
  
  br_x  <- bl_x + W
  top_x <- bl_x + W / 2
  
  outline <- data.frame(
    x = c(bl_x, br_x, top_x, bl_x),
    y = c(0,    0,    H,     0)
  )
  
  grid_lines <- do.call(rbind, lapply(c(0.2, 0.4, 0.6, 0.8), function(f) {
    
    # Parallel to bottom
    x1 <- bl_x  + f * (top_x - bl_x)
    x2 <- br_x  + f * (top_x - br_x)
    l1 <- data.frame(x = c(x1, x2), y = c(f * H, f * H),
                     grp = paste0("h_", f, "_", bl_x))
    
    # Parallel to left side
    x1 <- bl_x  + f * (br_x  - bl_x)
    x2 <- top_x + f * (br_x  - top_x)
    y2 <- H     + f * (0     - H)
    l2 <- data.frame(x = c(x1, x2), y = c(0, y2),
                     grp = paste0("l_", f, "_", bl_x))
    
    # Parallel to right side
    x1 <- br_x  + f * (bl_x  - br_x)
    x2 <- top_x + f * (bl_x  - top_x)
    y2 <- H     + f * (0     - H)
    l3 <- data.frame(x = c(x1, x2), y = c(0, y2),
                     grp = paste0("r_", f, "_", bl_x))
    
    rbind(l1, l2, l3)
  }))
  
  list(outline = outline, grid = grid_lines)
}

make_diamond <- function() {
  
  # Vertices aligned exactly with triangle corners:
  # D_BOT  = midpoint between L_BR and R_BL (bottom)
  # L_TOP  = left  triangle top (left  vertex of diamond)
  # R_TOP  = right triangle top (right vertex of diamond)
  # D_TOP  = projected above center (top vertex)
  
  outline <- data.frame(
    x = c(L_TOP[1], D_BOT[1], R_TOP[1], D_TOP[1], L_TOP[1]),
    y = c(L_TOP[2], D_BOT[2], R_TOP[2], D_TOP[2], L_TOP[2])
  )
  
  grid_lines <- do.call(rbind, lapply(c(0.2, 0.4, 0.6, 0.8), function(f) {
    
    # Parallel to left side (L_TOP → D_BOT)
    x1 <- L_TOP[1] + f * (D_BOT[1]  - L_TOP[1])
    y1 <- L_TOP[2] + f * (D_BOT[2]  - L_TOP[2])
    x2 <- D_TOP[1] + f * (R_TOP[1]  - D_TOP[1])
    y2 <- D_TOP[2] + f * (R_TOP[2]  - D_TOP[2])
    l1 <- data.frame(x = c(x1, x2), y = c(y1, y2),
                     grp = paste0("dl_", f))
    
    # Parallel to right side (R_TOP → D_BOT)
    x1 <- R_TOP[1] + f * (D_BOT[1]  - R_TOP[1])
    y1 <- R_TOP[2] + f * (D_BOT[2]  - R_TOP[2])
    x2 <- D_TOP[1] + f * (L_TOP[1]  - D_TOP[1])
    y2 <- D_TOP[2] + f * (L_TOP[2]  - D_TOP[2])
    l2 <- data.frame(x = c(x1, x2), y = c(y1, y2),
                     grp = paste0("dr_", f))
    
    rbind(l1, l2)
  }))
  
  list(
    outline = outline,
    grid    = grid_lines
  )
}

# =============================================================================
# MAIN PLOT FUNCTION
# =============================================================================

plot_piper_gg <- function(tern_data,
                          data         = NULL,
                          species_name = NULL,
                          title        = NULL,
                          color_by     = "count",
                          point_size   = 3,
                          alpha        = 0.75) {
  
  # Get counts
  if (!is.null(species_name) && !is.null(data)) {
    sp_counts <- data |>
      as_tibble() |>
      dplyr::filter(species == species_name) |>
      dplyr::group_by(year, region) |>
      dplyr::summarise(count = sum(count, na.rm = TRUE),
                       .groups = "drop")
    plot_df  <- tern_data |>
      dplyr::left_join(sp_counts, by = c("year", "region"))
    sp_label <- toupper(species_name)
  } else {
    plot_df  <- tern_data |>
      dplyr::rename(count = total_count)
    sp_label <- "All Species"
  }
  
  plot_df <- plot_df |>
    dplyr::filter(!is.na(count), count >= 0)
  
  max_count <- max(plot_df$count, na.rm = TRUE)
  
  # Coordinates
  lc <- left_xy(plot_df$T1_dry,
                plot_df$T1_wet,
                plot_df$T1_dynamic)
  
  rc <- right_xy(plot_df$T2_concentration,
                 plot_df$T2_extent,
                 plot_df$T2_access)
  
  dc <- diamond_xy(plot_df$T1_dry,
                   plot_df$T1_wet,
                   plot_df$T1_dynamic,
                   plot_df$T2_concentration,
                   plot_df$T2_extent,
                   plot_df$T2_access)
  
  # All points
  pts <- rbind(
    data.frame(x = lc$x, y = lc$y,
               count  = plot_df$count,
               region = plot_df$region,
               year   = plot_df$year,
               panel  = "left"),
    data.frame(x = rc$x, y = rc$y,
               count  = plot_df$count,
               region = plot_df$region,
               year   = plot_df$year,
               panel  = "right"),
    data.frame(x = dc$x, y = dc$y,
               count  = plot_df$count,
               region = plot_df$region,
               year   = plot_df$year,
               panel  = "diamond")
  )
  
  # Color variable
  if (color_by == "region") {
    pts$cv   <- as.numeric(factor(pts$region))
    clabel   <- "Region"
  } else if (color_by == "year") {
    pts$cv   <- pts$year
    clabel   <- "Year"
  } else {
    pts$cv   <- pts$count
    clabel   <- "Bird Count"
  }
  
  # Size proportional to count
  pts$psz <- scales::rescale(
    pts$count,
    to   = c(0.8, point_size),
    from = c(0, max_count)
  )
  
  # Structure
  left_tri  <- make_triangle(bl_x = L_BL[1])
  right_tri <- make_triangle(bl_x = R_BL[1])
  diam      <- make_diamond()
  
  # ==========================================================================
  # LABELS
  # ==========================================================================
  
  # Triangle corner labels
  cl <- data.frame(
    x  = c(
      L_BL[1],  L_BR[1],  L_TOP[1],    # left  tri
      R_BL[1],  R_BR[1],  R_TOP[1],    # right tri
      L_TOP[1] - 0.14,                  # diamond left
      R_TOP[1] + 0.14,                  # diamond right
      D_TOP[1],                         # diamond top
      D_BOT[1]                          # diamond bottom
    ),
    y  = c(
      -0.09, -0.09, H + 0.08,
      -0.09, -0.09, H + 0.08,
      L_TOP[2],
      R_TOP[2],
      D_TOP[2] + 0.10,
      D_BOT[2] - 0.10
    ),
    lab = c(
      "Dry",   "Wet",    "Dynamic",
      "Concentration", "Extent", "Access",
      "Dry &\nConcentrated",
      "Wet &\nExtensive",
      "Dynamic &\nHigh Access",
      "Static &\nLow Forage"
    ),
    col = c(
      rep("grey20", 6),
      "#c0392b", "#2980b9", "#8e44ad", "#7f8c8d"
    ),
    sz  = c(rep(3.8, 6), rep(2.8, 4)),
    ff  = c(rep("bold", 6), rep("italic", 4)),
    stringsAsFactors = FALSE
  )
  
  # Panel titles
  ptl <- data.frame(
    x   = c(L_BL[1] + W / 2,   R_BL[1] + W / 2),
    y   = c(-0.20,               -0.20),
    lab = c("WATER REGIME",     "FORAGING OPPORTUNITY"),
    col = c("#2980b9",            "#e67e22"),
    stringsAsFactors = FALSE
  )
  
  # ==========================================================================
  # BUILD PLOT
  # ==========================================================================
  
  p <- ggplot() +
    
    # ---- FILLS ----
  geom_polygon(
    data = left_tri$outline,
    aes(x = x, y = y),
    fill = "#d6eaf8", colour = NA, alpha = 0.35,
    inherit.aes = FALSE
  ) +
    geom_polygon(
      data = right_tri$outline,
      aes(x = x, y = y),
      fill = "#fdebd0", colour = NA, alpha = 0.35,
      inherit.aes = FALSE
    ) +
    geom_polygon(
      data = diam$outline,
      aes(x = x, y = y),
      fill = "#e8daef", colour = NA, alpha = 0.35,
      inherit.aes = FALSE
    ) +
    
    # ---- GRID LINES ----
  geom_line(
    data = left_tri$grid,
    aes(x = x, y = y, group = grp),
    colour = "#aed6f1", linetype = "dashed",
    linewidth = 0.3, inherit.aes = FALSE
  ) +
    geom_line(
      data = right_tri$grid,
      aes(x = x, y = y, group = grp),
      colour = "#fad7a0", linetype = "dashed",
      linewidth = 0.3, inherit.aes = FALSE
    ) +
    geom_line(
      data = diam$grid,
      aes(x = x, y = y, group = grp),
      colour = "#d2b4de", linetype = "dashed",
      linewidth = 0.3, inherit.aes = FALSE
    ) +
    
    # ---- OUTLINES ----
  geom_polygon(
    data = left_tri$outline,
    aes(x = x, y = y),
    fill = NA, colour = "#2980b9",
    linewidth = 1.2, inherit.aes = FALSE
  ) +
    geom_polygon(
      data = right_tri$outline,
      aes(x = x, y = y),
      fill = NA, colour = "#e67e22",
      linewidth = 1.2, inherit.aes = FALSE
    ) +
    geom_polygon(
      data = diam$outline,
      aes(x = x, y = y),
      fill = NA, colour = "#8e44ad",
      linewidth = 1.2, inherit.aes = FALSE
    ) +
    
    # ---- DATA POINTS ----
  geom_point(
    data = pts,
    aes(x = x, y = y, colour = cv, size = psz),
    alpha = alpha,
    inherit.aes = FALSE
  ) +
    scale_colour_gradient2(
      low      = "#3498db",
      mid      = "#f7dc6f",
      high     = "#e74c3c",
      midpoint = median(pts$cv, na.rm = TRUE),
      name     = clabel,
      labels   = scales::comma
    ) +
    scale_size_identity(guide = "none") +
    
    # ---- CORNER LABELS ----
  geom_text(
    data = cl,
    aes(x = x, y = y, label = lab),
    colour   = cl$col,
    size     = cl$sz,
    fontface = cl$ff,
    inherit.aes = FALSE
  ) +
    
    # ---- PANEL TITLES ----
  geom_text(
    data = ptl,
    aes(x = x, y = y, label = lab),
    colour   = ptl$col,
    size     = 4.2,
    fontface = "bold",
    inherit.aes = FALSE
  ) +
    
    labs(
      title = title %||%
        sprintf("Ecological Piper Diagram - %s", sp_label)
    ) +
    
    coord_fixed() +
    theme_void() +
    theme(
      plot.title      = element_text(face   = "bold",
                                     size   = 14,
                                     hjust  = 0.5,
                                     margin = margin(b = 10)),
      legend.position = "right",
      legend.title    = element_text(face = "bold", size = 10),
      legend.text     = element_text(size = 9),
      plot.margin     = margin(20, 20, 30, 20)
    )
  
  return(p)
}

# =============================================================================
# RUN AND SAVE
# =============================================================================

out_dir <- "results/piper_plots"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

species_list <- c("gbhe", "greg", "rosp", "sneg", "whib", "wost")

cat("=== Generating Piper diagrams (ggplot2) ===\n\n")

# Per species
cat("Per species...\n")
for (sp in species_list) {
  cat(" ", sp, "...\n")
  p <- tryCatch(
    plot_piper_gg(tern_data    = tern_data,
                  data         = data_sub,
                  species_name = sp,
                  color_by     = "count"),
    error = function(e) {
      cat("  Error:", e$message, "\n")
      NULL
    }
  )
  if (!is.null(p)) {
    ggsave(
      file.path(out_dir, sprintf("piper_%s.png", sp)),
      plot = p, width = 12, height = 10,
      dpi = 300, bg = "white"
    )
    cat("  Saved\n")
  }
}

# All species by region
cat("All species (by region)...\n")
p <- tryCatch(
  plot_piper_gg(tern_data = tern_data,
                color_by  = "region"),
  error = function(e) {
    cat("Error:", e$message, "\n")
    NULL
  }
)
if (!is.null(p)) {
  ggsave(file.path(out_dir, "piper_all_by_region.png"),
         p, width = 12, height = 10, dpi = 300, bg = "white")
  cat("Saved\n")
}

# All species by year
cat("All species (by year)...\n")
p <- tryCatch(
  plot_piper_gg(tern_data = tern_data,
                color_by  = "year"),
  error = function(e) {
    cat("Error:", e$message, "\n")
    NULL
  }
)
if (!is.null(p)) {
  ggsave(file.path(out_dir, "piper_all_by_year.png"),
         p, width = 12, height = 10, dpi = 300, bg = "white")
  cat("Saved\n")
}

# All species grid
cat("All species grid...\n")
plot_list <- lapply(species_list, function(sp) {
  tryCatch(
    plot_piper_gg(tern_data    = tern_data,
                  data         = data_sub,
                  species_name = sp,
                  color_by     = "count",
                  point_size   = 2) +
      theme(legend.position = "none",
            plot.title = element_text(size = 9)),
    error = function(e) NULL
  )
})

plot_list <- Filter(Negate(is.null), plot_list)

if (length(plot_list) > 0) {
  grid_plot <- gridExtra::grid.arrange(
    grobs = plot_list,
    ncol  = 3,
    top   = grid::textGrob(
      "Ecological Piper Diagrams - All Species",
      gp = grid::gpar(fontsize = 14, fontface = "bold")
    )
  )
  ggsave(
    file.path(out_dir, "piper_all_species_grid.png"),
    plot   = grid_plot,
    width  = 24,
    height = 18,
    dpi    = 200,
    bg     = "white"
  )
  cat("Grid saved\n")
}

cat("\n=== Complete ===\n")
cat("Saved to:", out_dir, "\n\n")
print(list.files(out_dir, pattern = "\\.png$"))