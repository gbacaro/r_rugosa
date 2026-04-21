# 10_build_paper_outputs_revision1.R
# Rebuild the main-text and supplementary tables and figures used for the manuscript submission package.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(sf)
  library(terra)
  library(ggplot2)
  library(stringr)
  library(scales)
  library(patchwork)
})

# ------------------------------------------------------------
# 1. Project paths
# ------------------------------------------------------------
project_dir <- getwd()

# Input folders
input_clean   <- file.path(project_dir, "data_clean")
input_tables  <- file.path(project_dir, "outputs", "tables")
input_models  <- file.path(project_dir, "models")
input_climate <- file.path(project_dir, "data_raw", "climate_chelsa", "current")
input_future  <- file.path(project_dir, "models", "future_projections")

# Output folders
out_root <- file.path(project_dir, "outputs", "paper_revision1")
out_main_fig <- file.path(out_root, "MAIN_TEXT", "FIGURES")
out_main_tab <- file.path(out_root, "MAIN_TEXT", "TABLES")
out_sup_fig  <- file.path(out_root, "SUPPLEMENTARY", "FIGURES")
out_sup_tab  <- file.path(out_root, "SUPPLEMENTARY", "TABLES")

for (d in c(out_root, out_main_fig, out_main_tab, out_sup_fig, out_sup_tab)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ------------------------------------------------------------
# 2. Input files
# ------------------------------------------------------------
occ_all_file      <- file.path(input_clean,  "occurrences_native_invaded_all.csv")
occ_coast_file    <- file.path(input_clean,  "occurrences_native_invaded_coastal.csv")
bg_file           <- file.path(input_clean,  "background_points_europe_coastal.csv")

cv_file           <- file.path(input_tables, "table_model_performance_cv.csv")
cv_sum_file       <- file.path(input_tables, "table_model_performance_summary.csv")
future_sum_file   <- file.path(input_tables, "table_future_projection_summary.csv")
current_sum_file  <- file.path(input_tables, "table_current_prediction_summary.csv")
area_sum_file     <- file.path(input_tables, "table_suitable_area_above_threshold.csv")
main_pred_file    <- file.path(input_tables, "table_main_predictors.csv")
niche_sum_file    <- file.path(input_tables, "table_niche_summary.csv")
niche_vars_file   <- file.path(input_tables, "table_niche_variables.csv")

current_mean_file <- file.path(input_models, "current_ensemble_mean.tif")

clim_eu_file      <- file.path(input_climate, "current_europe_coast1km_selectedvars.tif")
clim_asia_file    <- file.path(input_climate, "current_asia_coast1km_selectedvars.tif")

required_files <- c(
  occ_all_file, occ_coast_file, bg_file,
  cv_file, cv_sum_file, future_sum_file, current_sum_file,
  area_sum_file, main_pred_file, current_mean_file,
  clim_eu_file, clim_asia_file
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required file(s):\n", paste(missing_files, collapse = "\n"))
}

# ------------------------------------------------------------
# 3. Styling helpers
# ------------------------------------------------------------
col_native   <- "#0072B2"   # Okabe-Ito blue
col_invaded  <- "#D55E00"   # Okabe-Ito vermillion
col_bg       <- "#56B4E9"   # light blue
col_land     <- "grey85"
col_border   <- "white"
col_na       <- "grey95"

pal_seq <- scale_fill_viridis_c(option = "C", na.value = col_na)
pal_div <- scale_fill_gradient2(
  low = col_native,
  mid = "grey97",
  high = col_invaded,
  midpoint = 0,
  na.value = col_na
)
pal_cons <- scale_fill_gradient2(
  low = col_native,
  mid = "grey97",
  high = col_invaded,
  midpoint = 0,
  na.value = col_na,
  limits = c(-5, 5)
)

theme_editorial <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.2, colour = "grey85"),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = base_size + 1),
      plot.subtitle = element_text(size = base_size),
      strip.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )
}

save_plot <- function(p, filename, width = 8, height = 5.5, dpi = 350) {
  ggsave(filename, plot = p, width = width, height = height, dpi = dpi, bg = "white")
}

message_section <- function(txt) {
  cat("\n", paste(rep("=", nchar(txt) + 4), collapse = ""), "\n", sep = "")
  cat("  ", txt, "\n", sep = "")
  cat(paste(rep("=", nchar(txt) + 4), collapse = ""), "\n\n", sep = "")
}

pretty_alg <- function(x) recode(x,
                                 GLM = "GLM",
                                 GAM = "GAM",
                                 RF = "Random forest",
                                 MAXNET = "Maxnet",
                                 .default = x
)

var_lookup <- tibble(
  variable = c("bio3","bio8","bio9","bio10","bio15","bio18","bio19","fcf","gsl","gst"),
  variable_name = c(
    "Isothermality",
    "Mean temperature of wettest quarter",
    "Mean temperature of driest quarter",
    "Mean temperature of warmest quarter",
    "Precipitation seasonality",
    "Precipitation of warmest quarter",
    "Precipitation of coldest quarter",
    "Frost change frequency",
    "Growing season length",
    "Growing season temperature"
  )
)

safe_read <- function(path) readr::read_csv(path, show_col_types = FALSE)

aggregate_for_plot <- function(r, max_cells = 200000) {
  n <- terra::ncell(r)
  if (!is.finite(n) || n <= max_cells) return(r)
  fact <- ceiling(sqrt(n / max_cells))
  terra::aggregate(r, fact = fact, fun = mean, na.rm = TRUE)
}

load_world <- function() {
  if (requireNamespace("rnaturalearth", quietly = TRUE)) {
    world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
    return(world)
  }
  stop("Package 'rnaturalearth' is required for the revised maps.")
}

# ------------------------------------------------------------
# 4. Read inputs
# ------------------------------------------------------------
message_section("Reading inputs")

occ_all   <- safe_read(occ_all_file)
occ_coast <- safe_read(occ_coast_file)
bg        <- safe_read(bg_file)
cv_tbl    <- safe_read(cv_file)
cv_sum    <- safe_read(cv_sum_file)
future_sum <- safe_read(future_sum_file)
current_sum <- safe_read(current_sum_file)
area_sum  <- safe_read(area_sum_file)
main_pred <- safe_read(main_pred_file)

niche_sum <- if (file.exists(niche_sum_file)) safe_read(niche_sum_file) else NULL
niche_vars <- if (file.exists(niche_vars_file)) safe_read(niche_vars_file) else NULL

current_mean <- terra::rast(current_mean_file)
clim_eu <- terra::rast(clim_eu_file)
clim_asia <- terra::rast(clim_asia_file)

world <- load_world()

# ------------------------------------------------------------
# 5. MAIN TABLE 1: occurrence / domain summary
# ------------------------------------------------------------
message_section("Building main Table 1")

range_year_summary <- occ_coast %>%
  filter(!is.na(year)) %>%
  group_by(range_status) %>%
  summarise(
    first_year = min(year, na.rm = TRUE),
    last_year = max(year, na.rm = TRUE),
    .groups = "drop"
  )

all_counts <- occ_all %>% count(range_status, name = "raw_records")
coast_counts <- occ_coast %>% count(range_status, name = "coastal_records_used")

main_table_1 <- tibble(
  dataset_component = c("Native range", "Invaded range", "European coastal background"),
  raw_records = c(
    all_counts$raw_records[match("native", all_counts$range_status)],
    all_counts$raw_records[match("invaded", all_counts$range_status)],
    NA_real_
  ),
  coastal_records_used = c(
    coast_counts$coastal_records_used[match("native", coast_counts$range_status)],
    coast_counts$coastal_records_used[match("invaded", coast_counts$range_status)],
    nrow(bg)
  ),
  temporal_coverage = c(
    paste0(range_year_summary$first_year[match("native", range_year_summary$range_status)], "–",
           range_year_summary$last_year[match("native", range_year_summary$range_status)]),
    paste0(range_year_summary$first_year[match("invaded", range_year_summary$range_status)], "–",
           range_year_summary$last_year[match("invaded", range_year_summary$range_status)]),
    "NA"
  ),
  notes = c(
    "Coastal records used for native-range niche comparison.",
    "Coastal records used for European SDM calibration and niche comparison.",
    "Background points sampled within the European coastal climatic domain."
  )
)

write_csv(main_table_1, file.path(out_main_tab, "Table_1_occurrence_and_domain_summary.csv"))

# ------------------------------------------------------------
# 6. MAIN TABLE 2: model performance summary
# ------------------------------------------------------------
message_section("Building main Table 2")

fmt_pm <- function(mean, sd, digits = 3) sprintf(paste0("%.", digits, "f ± %.", digits, "f"), mean, sd)

main_table_2 <- cv_sum %>%
  transmute(
    algorithm = pretty_alg(algorithm),
    `AUC (mean ± SD)` = fmt_pm(auc_mean, auc_sd, 3),
    `TSS (mean ± SD)` = fmt_pm(tss_mean, tss_sd, 3),
    `Sensitivity (mean ± SD)` = fmt_pm(sensitivity_mean, sensitivity_sd, 3),
    `Specificity (mean ± SD)` = fmt_pm(specificity_mean, specificity_sd, 3),
    `Boyce index (mean ± SD)` = fmt_pm(boyce_mean, boyce_sd, 3)
  )

write_csv(main_table_2, file.path(out_main_tab, "Table_2_model_performance_summary.csv"))

# ------------------------------------------------------------
# 7. MAIN TABLE 3: current + SSP x period summary
# ------------------------------------------------------------
message_section("Building main Table 3")

agg_future <- future_sum %>%
  group_by(ssp, time_period) %>%
  summarise(
    mean_suitability = mean(mean_suitability_mean, na.rm = TRUE),
    delta_suitability = mean(delta_mean_suitability, na.rm = TRUE),
    mean_uncertainty_sd = mean(uncertainty_sd_mean, na.rm = TRUE),
    mean_uncertainty_cv = mean(uncertainty_cv_mean, na.rm = TRUE),
    mean_agreement = mean(agreement_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    scenario = case_when(
      ssp == "ssp126" & time_period == "2041-2070" ~ "SSP1-2.6 (2041–2070)",
      ssp == "ssp126" & time_period == "2071-2100" ~ "SSP1-2.6 (2071–2100)",
      ssp == "ssp585" & time_period == "2041-2070" ~ "SSP5-8.5 (2041–2070)",
      ssp == "ssp585" & time_period == "2071-2100" ~ "SSP5-8.5 (2071–2100)",
      TRUE ~ paste(ssp, time_period)
    )
  ) %>%
  select(scenario, mean_suitability, delta_suitability, mean_uncertainty_sd, mean_uncertainty_cv, mean_agreement)

main_table_3 <- bind_rows(
  tibble(
    scenario = "Current",
    mean_suitability = current_sum$mean_suitability[1],
    delta_suitability = 0,
    mean_uncertainty_sd = current_sum$uncertainty_sd_mean[1],
    mean_uncertainty_cv = current_sum$uncertainty_cv_mean[1],
    mean_agreement = NA_real_
  ),
  agg_future
)

write_csv(main_table_3, file.path(out_main_tab, "Table_3_current_and_future_summary_by_scenario_class.csv"))

# ------------------------------------------------------------
# 8. SUPPLEMENTARY TABLES
# ------------------------------------------------------------
message_section("Building supplementary tables")

supp_table_s1 <- future_sum %>%
  transmute(
    gcm,
    ssp,
    time_period,
    mean_suitability = mean_suitability_mean,
    delta_suitability = delta_mean_suitability,
    uncertainty_sd = uncertainty_sd_mean,
    uncertainty_cv = uncertainty_cv_mean,
    mean_agreement = agreement_mean,
    max_agreement = agreement_max
  )
write_csv(supp_table_s1, file.path(out_sup_tab, "Table_S1_future_projection_summary_by_GCM.csv"))

niche_var_vec <- if (!is.null(niche_vars) && 'variable' %in% names(niche_vars)) niche_vars$variable else character()
supp_table_s2 <- var_lookup %>%
  mutate(
    used_in_niche_analysis = ifelse(variable %in% niche_var_vec, "Yes", "No"),
    retained_in_final_SDM = ifelse(variable %in% main_pred$variable, "Yes", "No")
  )
write_csv(supp_table_s2, file.path(out_sup_tab, "Table_S2_climatic_predictors_used.csv"))

supp_table_s3 <- cv_tbl %>%
  mutate(algorithm = pretty_alg(algorithm))
write_csv(supp_table_s3, file.path(out_sup_tab, "Table_S3_fold_level_model_performance.csv"))

supp_table_s4 <- area_sum %>%
  mutate(
    scenario = case_when(
      ssp == "ssp126" & time_period == "2041-2070" ~ "SSP1-2.6 (2041–2070)",
      ssp == "ssp126" & time_period == "2071-2100" ~ "SSP1-2.6 (2071–2100)",
      ssp == "ssp585" & time_period == "2041-2070" ~ "SSP5-8.5 (2041–2070)",
      ssp == "ssp585" & time_period == "2071-2100" ~ "SSP5-8.5 (2071–2100)",
      TRUE ~ paste(ssp, time_period)
    )
  ) %>%
  select(gcm, scenario, threshold, valid_domain_area_km2, current_suitable_area_km2,
         current_suitable_area_pct, future_suitable_area_km2, future_suitable_area_pct,
         delta_area_km2, delta_area_pct)
write_csv(supp_table_s4, file.path(out_sup_tab, "Table_S4_threshold_based_suitable_area_summary.csv"))

# ------------------------------------------------------------
# 9. FIGURE 1: occurrences, two-panel native vs invaded
# ------------------------------------------------------------
message_section("Building Figure 1")

world_4326 <- st_transform(world, 4326)

# Use cropping rather than intersection to avoid unnecessary geometry splitting
world_native <- st_crop(world_4326, xmin = 120, xmax = 156, ymin = 28, ymax = 56)
world_invaded <- st_crop(world_4326, xmin = -12, xmax = 46, ymin = 34, ymax = 72)

occ_native <- occ_coast %>%
  filter(range_status == "native")

occ_invaded <- occ_coast %>%
  filter(range_status == "invaded")

p_native <- ggplot() +
  geom_sf(data = world_native, fill = col_land, color = col_border, linewidth = 0.2) +
  geom_point(
    data = occ_native,
    aes(decimal_longitude, decimal_latitude, color = range_status),
    size = 0.7,
    alpha = 0.75
  ) +
  scale_color_manual(
    values = c(native = col_native, invaded = col_invaded),
    labels = c(native = "Native", invaded = "Invaded"),
    drop = FALSE
  ) +
  coord_sf(
    xlim = c(120, 156),
    ylim = c(28, 56),
    expand = FALSE
  ) +
  labs(title = "Native East Asia", x = NULL, y = NULL, color = NULL) +
  theme_editorial() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

p_invaded <- ggplot() +
  geom_sf(data = world_invaded, fill = col_land, color = col_border, linewidth = 0.2) +
  geom_point(
    data = occ_invaded,
    aes(decimal_longitude, decimal_latitude, color = range_status),
    size = 0.7,
    alpha = 0.75
  ) +
  scale_color_manual(
    values = c(native = col_native, invaded = col_invaded),
    labels = c(native = "Native", invaded = "Invaded"),
    drop = FALSE
  ) +
  coord_sf(
    xlim = c(-12, 46),
    ylim = c(34, 72),
    expand = FALSE
  ) +
  labs(title = "Invaded Europe", x = NULL, y = NULL, color = NULL) +
  theme_editorial() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

p_fig1 <- (p_native + p_invaded) +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

save_plot(
  p_fig1,
  file.path(out_main_fig, "Figure_1_occurrences_native_vs_invaded.png"),
  width = 10.5,
  height = 5.8
)

# ------------------------------------------------------------
# 10. FIGURE 2: PCA climatic space (rebuilt from climatic data)
# ------------------------------------------------------------
message_section("Building Figure 2")

# Use the variables from niche analysis if available, otherwise fallback
pca_vars <- if (!is.null(niche_vars) && "variable" %in% names(niche_vars)) niche_vars$variable else c("bio3","bio8","bio9","bio10","bio15","bio18","bio19","fcf","gsl","gst")
pca_vars <- intersect(pca_vars, intersect(names(clim_eu), names(clim_asia)))
if (length(pca_vars) < 3) stop("Too few niche-analysis variables available to rebuild PCA figure.")

sample_valid_values <- function(r, vars, n_max = 8000) {
  x <- as.data.frame(r[[vars]], na.rm = TRUE)
  if (nrow(x) > n_max) {
    set.seed(123)
    x <- dplyr::slice_sample(x, n = n_max)
  }
  x
}

bg_eu_vals   <- sample_valid_values(clim_eu, pca_vars, n_max = 8000)
bg_asia_vals <- sample_valid_values(clim_asia, pca_vars, n_max = 8000)

pca_bg <- bind_rows(
  mutate(bg_asia_vals, region = "Native background"),
  mutate(bg_eu_vals, region = "Invaded background")
)

pca_fit <- prcomp(pca_bg[, pca_vars], center = TRUE, scale. = TRUE)

extract_occ_env <- function(df, r, vars) {
  xy <- df[, c("decimal_longitude","decimal_latitude")]
  ex <- terra::extract(r[[vars]], xy)[, -1, drop = FALSE]
  bind_cols(df, as_tibble(ex)) %>%
    dplyr::filter(dplyr::if_all(dplyr::all_of(vars), ~ !is.na(.)))
}

occ_native_env  <- occ_coast %>% filter(range_status == "native")  %>% extract_occ_env(clim_asia, pca_vars) %>% mutate(group = "Native")
occ_invaded_env <- occ_coast %>% filter(range_status == "invaded") %>% extract_occ_env(clim_eu, pca_vars)   %>% mutate(group = "Invaded")
occ_pca <- bind_rows(occ_native_env, occ_invaded_env)

scores <- predict(pca_fit, newdata = occ_pca[, pca_vars])
occ_scores <- bind_cols(occ_pca %>% select(group), as.data.frame(scores)[, 1:2]) %>%
  rename(PC1 = PC1, PC2 = PC2)

# Convex hulls for cleaner shape envelopes
hull_df <- occ_scores %>%
  group_by(group) %>%
  slice(chull(PC1, PC2)) %>%
  arrange(group, PC1, PC2) %>%
  ungroup()

p_fig2 <- ggplot(occ_scores, aes(PC1, PC2, color = group)) +
  geom_point(alpha = 0.35, size = 0.8) +
  geom_polygon(
    data = hull_df,
    aes(x = PC1, y = PC2, group = group, fill = group),
    alpha = 0.10,
    color = NA,
    inherit.aes = FALSE
  ) +
  stat_ellipse(aes(group = group), linewidth = 1.0, level = 0.90, show.legend = FALSE) +
  scale_color_manual(values = c("Native" = col_native, "Invaded" = col_invaded)) +
  scale_fill_manual(values = c("Native" = col_native, "Invaded" = col_invaded)) +
  labs(x = "PC1", y = "PC2", color = NULL, fill = NULL) +
  theme_editorial() +
  theme(legend.position = "bottom")

save_plot(p_fig2, file.path(out_main_fig, "Figure_2_climatic_space_PCA.png"), width = 7.5, height = 6.2)

# ------------------------------------------------------------
# 11. Helpers for projected Europe maps
# ------------------------------------------------------------
make_europe_base <- function() {
  world %>%
    st_transform(3035) %>%
    st_crop(xmin = 1800000, xmax = 7200000, ymin = 1200000, ymax = 5600000)
}

raster_df_3035 <- function(r, max_cells = 180000) {
  rp <- terra::project(aggregate_for_plot(r, max_cells = max_cells), "EPSG:3035")
  as.data.frame(rp, xy = TRUE, na.rm = FALSE)
}

europe_base <- make_europe_base()

# ------------------------------------------------------------
# 12. FIGURE 3: current suitability
# ------------------------------------------------------------
message_section("Building Figure 3")

cur_df <- raster_df_3035(current_mean)
names(cur_df)[3] <- "suitability"

p_fig3 <- ggplot() +
  geom_sf(data = europe_base, fill = col_land, color = col_border, linewidth = 0.2) +
  geom_raster(data = cur_df, aes(x = x, y = y, fill = suitability)) +
  pal_seq +
  coord_sf(expand = FALSE) +
  labs(x = NULL, y = NULL, fill = "Suitability") +
  theme_editorial()

save_plot(p_fig3, file.path(out_main_fig, "Figure_3_current_suitability_map.png"), width = 8.6, height = 6.3)

# ------------------------------------------------------------
# 13. FIGURE 4: mean delta suitability by SSP x period
# ------------------------------------------------------------
message_section("Building Figure 4")

scenario_groups <- tribble(
  ~ssp, ~time_period, ~label,
  "ssp126", "2041-2070", "SSP1-2.6\n2041–2070",
  "ssp126", "2071-2100", "SSP1-2.6\n2071–2100",
  "ssp585", "2041-2070", "SSP5-8.5\n2041–2070",
  "ssp585", "2071-2100", "SSP5-8.5\n2071–2100"
)

build_group_mean_delta <- function(ssp_i, time_i) {
  rows <- future_sum %>% filter(ssp == ssp_i, time_period == time_i)
  if (nrow(rows) == 0) stop("No future scenarios found for ", ssp_i, " ", time_i)
  rasters <- lapply(seq_len(nrow(rows)), function(i) {
    scen_tag <- paste(rows$gcm[i], rows$ssp[i], rows$time_period[i], sep = "_")
    terra::rast(file.path(input_future, scen_tag, paste0("future_ensemble_mean_", scen_tag, ".tif")))
  })
  r_stack <- rast(rasters)
  r_mean <- app(r_stack, mean, na.rm = TRUE)
  r_mean - current_mean
}

delta_list <- lapply(seq_len(nrow(scenario_groups)), function(i) {
  g <- scenario_groups[i, ]
  r_delta <- build_group_mean_delta(g$ssp, g$time_period)
  df <- raster_df_3035(r_delta)
  names(df)[3] <- "delta"
  df$panel <- g$label
  df
})
delta_df <- bind_rows(delta_list)

p_fig4 <- ggplot() +
  geom_sf(data = europe_base, fill = col_land, color = col_border, linewidth = 0.2) +
  geom_raster(data = delta_df, aes(x = x, y = y, fill = delta)) +
  pal_div +
  facet_wrap(~panel, ncol = 2) +
  coord_sf(expand = FALSE) +
  labs(x = NULL, y = NULL, fill = expression(Delta~suitability)) +
  theme_editorial()

save_plot(p_fig4, file.path(out_main_fig, "Figure_4_mean_delta_suitability_by_scenario_class.png"), width = 10, height = 7.6)

# ------------------------------------------------------------
# 14. FIGURE 5: directional consensus across GCMs
# ------------------------------------------------------------
message_section("Building Figure 5")

build_group_consensus <- function(ssp_i, time_i) {
  rows <- future_sum %>% filter(ssp == ssp_i, time_period == time_i)
  rasters <- lapply(seq_len(nrow(rows)), function(i) {
    scen_tag <- paste(rows$gcm[i], rows$ssp[i], rows$time_period[i], sep = "_")
    fut <- terra::rast(file.path(input_future, scen_tag, paste0("future_ensemble_mean_", scen_tag, ".tif")))
    d <- fut - current_mean
    pos <- terra::ifel(d > 0, 1, terra::ifel(d < 0, -1, 0))
    pos
  })
  app(rast(rasters), sum, na.rm = TRUE)
}

cons_list <- lapply(seq_len(nrow(scenario_groups)), function(i) {
  g <- scenario_groups[i, ]
  r_cons <- build_group_consensus(g$ssp, g$time_period)
  df <- raster_df_3035(r_cons)
  names(df)[3] <- "consensus"
  df$panel <- g$label
  df
})
cons_df <- bind_rows(cons_list)

p_fig5 <- ggplot() +
  geom_sf(data = europe_base, fill = col_land, color = col_border, linewidth = 0.2) +
  geom_raster(data = cons_df, aes(x = x, y = y, fill = consensus)) +
  pal_cons +
  facet_wrap(~panel, ncol = 2) +
  coord_sf(expand = FALSE) +
  labs(x = NULL, y = NULL, fill = "Consensus") +
  theme_editorial()

save_plot(p_fig5, file.path(out_main_fig, "Figure_5_directional_consensus_by_scenario_class.png"), width = 10, height = 7.6)

# ------------------------------------------------------------
# 15. SUPPLEMENTARY FIGURES
# ------------------------------------------------------------
message_section("Building supplementary figures")

# Figure S1: cross-validated model performance (AUC)
p_s1 <- ggplot(cv_tbl, aes(x = pretty_alg(algorithm), y = auc)) +
  geom_boxplot(width = 0.65, fill = "grey88", color = "grey20") +
  geom_jitter(width = 0.08, alpha = 0.7, size = 1.2, color = col_native) +
  labs(x = NULL, y = "AUC") +
  theme_editorial()

save_plot(p_s1, file.path(out_sup_fig, "Figure_S1_cross_validated_model_performance.png"), width = 7, height = 5)

# Figure S2: background design
bg_sf <- st_as_sf(bg, coords = c("decimal_longitude", "decimal_latitude"), crs = 4326) %>% st_transform(3035)
pres_sf <- occ_coast %>% filter(range_status == "invaded") %>%
  st_as_sf(coords = c("decimal_longitude", "decimal_latitude"), crs = 4326) %>%
  st_transform(3035)

p_s2 <- ggplot() +
  geom_sf(data = europe_base, fill = col_land, color = col_border, linewidth = 0.2) +
  geom_sf(data = bg_sf, color = col_bg, alpha = 0.35, size = 0.25) +
  geom_sf(data = pres_sf, color = col_invaded, alpha = 0.7, size = 0.35) +
  coord_sf(expand = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_editorial()

save_plot(p_s2, file.path(out_sup_fig, "Figure_S2_background_design.png"), width = 8.6, height = 6.3)

# ------------------------------------------------------------
# 16. Session manifest
# ------------------------------------------------------------
message_section("Writing manifest")
manifest <- tibble::tibble(
  folder = c(
    out_main_fig, out_main_tab, out_sup_fig, out_sup_tab
  ),
  description = c(
    "Main text figures",
    "Main text tables",
    "Supplementary figures",
    "Supplementary tables"
  )
)
write_csv(manifest, file.path(out_root, "manifest_revision1_outputs.csv"))

message_section("Done")
cat("Outputs written to:\n", out_root, "\n")