# 12_refresh_main_figures_3_5_final.R
# Rebuild the final main-text Figures 3–5 in the paper output directory.

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tibble)
  library(purrr)
  library(rnaturalearth)
  library(scales)
  library(grid)
})

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------
project_dir <- getwd()

dir_models   <- file.path(project_dir, "models")
dir_fig_main <- file.path(project_dir, "outputs", "paper_revision1", "MAIN_TEXT", "FIGURES")
dir_tables   <- file.path(project_dir, "outputs", "tables")
dir_future   <- file.path(dir_models, "future_projections")

dir.create(dir_fig_main, recursive = TRUE, showWarnings = FALSE)

current_mean_file   <- file.path(dir_models, "current_ensemble_mean.tif")
future_summary_file <- file.path(dir_tables, "table_future_projection_summary.csv")

if (!file.exists(current_mean_file))   stop("Missing file: ", current_mean_file)
if (!file.exists(future_summary_file)) stop("Missing file: ", future_summary_file)
if (!dir.exists(dir_future))           stop("Missing directory: ", dir_future)

# ------------------------------------------------------------
# 2. Helpers
# ------------------------------------------------------------
message_section <- function(txt) {
  cat("\n", paste(rep("=", nchar(txt) + 4), collapse = ""), "\n", sep = "")
  cat("  ", txt, "\n", sep = "")
  cat(paste(rep("=", nchar(txt) + 4), collapse = ""), "\n\n", sep = "")
}

save_plot <- function(p, filename, width = 10, height = 7, dpi = 350) {
  ggsave(
    filename = filename,
    plot = p,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
}

theme_editorial_map <- function() {
  theme_minimal(base_size = 12) +
    theme(
      panel.grid.major = element_line(color = "grey82", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      axis.title = element_blank(),
      axis.text = element_text(color = "grey25"),
      legend.title = element_text(face = "bold"),
      legend.position = "right",
      strip.text = element_text(face = "bold"),
      strip.background = element_blank(),
      plot.title = element_text(face = "bold")
    )
}

build_europe_basemap <- function() {
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  world_4326 <- st_transform(world, 4326)
  bbox_eu <- st_as_sfc(st_bbox(c(xmin = -12, xmax = 45, ymin = 34, ymax = 72), crs = st_crs(4326)))
  st_crop(world_4326, bbox_eu)
}

raster_to_df_ll <- function(r, value_name = "value", panel_name = NULL) {
  df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(df)[3] <- value_name
  if (!is.null(panel_name)) df$panel <- panel_name
  as_tibble(df)
}

# double square rendering for visibility
add_square_points <- function(p, data, mapping,
                              size_under = 1.8, size_top = 0.95,
                              alpha_under = 0.45, alpha_top = 0.98,
                              show.legend.under = FALSE) {
  p +
    geom_point(
      data = data,
      mapping = mapping,
      shape = 15,
      size = size_under,
      alpha = alpha_under,
      show.legend = show.legend.under
    ) +
    geom_point(
      data = data,
      mapping = mapping,
      shape = 15,
      size = size_top,
      alpha = alpha_top
    )
}

coord_europe <- coord_sf(
  xlim = c(-12, 45),
  ylim = c(34, 72),
  expand = FALSE,
  datum = NA,
  default_crs = st_crs(4326)
)

# ------------------------------------------------------------
# 3. Read inputs
# ------------------------------------------------------------
message_section("Reading inputs")

r_current <- rast(current_mean_file)
future_summary <- read_csv(future_summary_file, show_col_types = FALSE)
world_eu <- build_europe_basemap()

scenario_grid <- tibble::tribble(
  ~ssp,      ~time_period, ~panel,
  "ssp126",  "2041-2070",  "SSP1-2.6\n2041–2070",
  "ssp126",  "2071-2100",  "SSP1-2.6\n2071–2100",
  "ssp585",  "2041-2070",  "SSP5-8.5\n2041–2070",
  "ssp585",  "2071-2100",  "SSP5-8.5\n2071–2100"
)

load_future_mean_rasters <- function(future_summary, dir_future, ssp_value, time_value) {
  rows <- future_summary |>
    filter(ssp == ssp_value, time_period == time_value)
  
  if (nrow(rows) == 0) {
    stop("No rows in future summary for: ", ssp_value, " / ", time_value)
  }
  
  files <- file.path(
    dir_future,
    paste(rows$gcm, rows$ssp, rows$time_period, sep = "_"),
    paste0("future_ensemble_mean_", paste(rows$gcm, rows$ssp, rows$time_period, sep = "_"), ".tif")
  )
  
  missing_files <- files[!file.exists(files)]
  if (length(missing_files) > 0) {
    stop("Missing future raster(s):\n", paste(missing_files, collapse = "\n"))
  }
  
  rast(files)
}

# ------------------------------------------------------------
# 4. Figure 3 - current suitability
# ------------------------------------------------------------
message_section("Building Figure 3")

df_current <- raster_to_df_ll(r_current, "suitability")
max_current <- max(df_current$suitability, na.rm = TRUE)

p3 <- ggplot() +
  geom_sf(data = world_eu, fill = "grey80", color = "white", linewidth = 0.25)

p3 <- add_square_points(
  p = p3,
  data = df_current,
  mapping = aes(x = x, y = y, color = suitability),
  size_under = 1.9,
  size_top = 1.0,
  alpha_under = 0.45,
  alpha_top = 0.98
) +
  scale_color_viridis_c(
    option = "plasma",
    limits = c(0, max_current),
    oob = squish,
    breaks = c(0, 0.10, 0.25, 0.50, 0.75),
    labels = number_format(accuracy = 0.01),
    guide = guide_colorbar(
      title = "Suitability",
      barheight = unit(60, "pt"),
      frame.colour = "black",
      ticks.colour = "black"
    )
  ) +
  coord_europe +
  theme_editorial_map()

save_plot(
  p3,
  file.path(dir_fig_main, "Figure_3_current_suitability_map.png"),
  width = 10.5,
  height = 7.6
)

# ------------------------------------------------------------
# 5. Compute grouped rasters for Figures 4 and 5
# ------------------------------------------------------------
message_section("Computing grouped future deltas and consensus")

group_delta_list <- list()
group_consensus_list <- list()

for (i in seq_len(nrow(scenario_grid))) {
  ssp_i   <- scenario_grid$ssp[i]
  time_i  <- scenario_grid$time_period[i]
  panel_i <- scenario_grid$panel[i]
  
  r_future_stack <- load_future_mean_rasters(
    future_summary = future_summary,
    dir_future = dir_future,
    ssp_value = ssp_i,
    time_value = time_i
  )
  
  r_delta_stack <- r_future_stack - r_current
  r_mean_delta <- app(r_delta_stack, mean, na.rm = TRUE)
  r_consensus <- app(r_delta_stack, fun = function(x) sum(sign(x), na.rm = TRUE))
  
  group_delta_list[[panel_i]] <- r_mean_delta
  group_consensus_list[[panel_i]] <- r_consensus
}

# ------------------------------------------------------------
# 6. Figure 4 - mean delta suitability
# ------------------------------------------------------------
message_section("Building Figure 4")

delta_df_list <- purrr::imap(group_delta_list, function(r, panel_name) {
  raster_to_df_ll(r, "delta", panel_name)
})

delta_df <- bind_rows(delta_df_list)

lim_delta <- quantile(abs(delta_df$delta), probs = 0.98, na.rm = TRUE)
lim_delta <- max(0.05, as.numeric(lim_delta))

p4 <- ggplot() +
  geom_sf(data = world_eu, fill = "grey80", color = "white", linewidth = 0.25)

p4 <- add_square_points(
  p = p4,
  data = delta_df,
  mapping = aes(x = x, y = y, color = delta),
  size_under = 1.9,
  size_top = 1.0,
  alpha_under = 0.45,
  alpha_top = 0.98
) +
  scale_color_gradient2(
    low = "#2C7BB6",
    mid = "#F7F7F7",
    high = "#D95F0E",
    midpoint = 0,
    limits = c(-lim_delta, lim_delta),
    oob = squish,
    guide = guide_colorbar(
      title = expression(Delta~"suitability"),
      barheight = unit(60, "pt"),
      frame.colour = "black",
      ticks.colour = "black"
    )
  ) +
  facet_wrap(~panel, ncol = 2) +
  coord_europe +
  theme_editorial_map()

save_plot(
  p4,
  file.path(dir_fig_main, "Figure_4_mean_delta_suitability_by_scenario_class.png"),
  width = 11.5,
  height = 8.8
)

# ------------------------------------------------------------
# 7. Figure 5 - directional consensus (DISCRETE, NON-ZERO ONLY)
# ------------------------------------------------------------
message_section("Building Figure 5")

cons_df_list <- purrr::imap(group_consensus_list, function(r, panel_name) {
  raster_to_df_ll(r, "consensus", panel_name)
})

cons_df <- bind_rows(cons_df_list) |>
  filter(!is.na(consensus), consensus != 0) |>
  mutate(
    consensus = as.integer(round(consensus)),
    consensus_f = factor(
      consensus,
      levels = c(-5, -4, -3, -2, -1, 1, 2, 3, 4, 5)
    )
  )

consensus_cols <- c(
  `-5` = "#2166AC",
  `-4` = "#4393C3",
  `-3` = "#92C5DE",
  `-2` = "#BFD9EA",
  `-1` = "#DCEAF4",
  `1`  = "#FEE0C5",
  `2`  = "#FDBB84",
  `3`  = "#FC8D59",
  `4`  = "#E66101",
  `5`  = "#B35806"
)

p5 <- ggplot() +
  geom_sf(data = world_eu, fill = "grey80", color = "white", linewidth = 0.25)

p5 <- add_square_points(
  p = p5,
  data = cons_df,
  mapping = aes(x = x, y = y, color = consensus_f),
  size_under = 1.9,
  size_top = 1.0,
  alpha_under = 0.45,
  alpha_top = 0.98
) +
  scale_color_manual(
    values = consensus_cols,
    drop = FALSE,
    guide = guide_legend(
      title = "Consensus",
      override.aes = list(shape = 15, size = 4, alpha = 1)
    )
  ) +
  facet_wrap(~panel, ncol = 2) +
  coord_europe +
  theme_editorial_map()

save_plot(
  p5,
  file.path(dir_fig_main, "Figure_5_directional_consensus_by_scenario_class.png"),
  width = 11.5,
  height = 8.8
)

message_section("Figures 3–5 rebuilt successfully")

cat("Saved:\n")
cat(file.path(dir_fig_main, "Figure_3_current_suitability_map.png"), "\n")
cat(file.path(dir_fig_main, "Figure_4_mean_delta_suitability_by_scenario_class.png"), "\n")
cat(file.path(dir_fig_main, "Figure_5_directional_consensus_by_scenario_class.png"), "\n")