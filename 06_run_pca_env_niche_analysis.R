# 06_run_pca_env_niche_analysis.R
# Compare native and invaded coastal climatic space using the PCA-env framework and export niche-summary objects and figures.

library(terra)
library(sf)
library(readr)
library(dplyr)
library(ggplot2)
library(tibble)
library(ecospat)

# ------------------------------------------------------------
# 1. Local project settings
# ------------------------------------------------------------
project_dir <- getwd()

dir_clean     <- file.path(project_dir, "data_clean")
dir_models    <- file.path(project_dir, "models")
dir_fig       <- file.path(project_dir, "outputs", "figures")
dir_tab       <- file.path(project_dir, "outputs", "tables")
dir_raw_clim  <- file.path(project_dir, "data_raw", "climate_chelsa")

dir.create(dir_models, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_tab, recursive = TRUE, showWarnings = FALSE)

# Variables actually downloaded in script 04
candidate_vars <- c("bio3", "bio8", "bio9", "bio10", "bio15", "bio18", "bio19", "fcf", "gsl", "gst")

# Background sample size
bg_sample_n <- 30000

# ------------------------------------------------------------
# 2. Helper functions
# ------------------------------------------------------------
message_section <- function(txt) {
  cat("\n", paste(rep("=", nchar(txt) + 4), collapse = ""), "\n", sep = "")
  cat("  ", txt, "\n", sep = "")
  cat(paste(rep("=", nchar(txt) + 4), collapse = ""), "\n\n", sep = "")
}

theme_editorial <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.25),
      axis.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 2),
      plot.subtitle = ggplot2::element_text(size = base_size),
      legend.position = "right",
      strip.text = ggplot2::element_text(face = "bold")
    )
}

save_editorial_plot <- function(p, filename, width = 8, height = 5, dpi = 320) {
  ggplot2::ggsave(
    filename = filename,
    plot = p,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
}

thin_by_cell_df <- function(df, raster_template, lon = "decimal_longitude", lat = "decimal_latitude") {
  xy <- as.matrix(df[, c(lon, lat)])
  cell_id <- terra::cellFromXY(raster_template, xy)
  df$cell_id_tmp <- cell_id
  df <- df[!is.na(df$cell_id_tmp), , drop = FALSE]
  df <- df |>
    dplyr::group_by(cell_id_tmp) |>
    dplyr::slice(1) |>
    dplyr::ungroup()
  df$cell_id_tmp <- NULL
  df
}

# ------------------------------------------------------------
# 3. Input files
# ------------------------------------------------------------
occ_file      <- file.path(dir_clean, "occurrences_native_invaded_coastal.csv")
bio_eu_file   <- file.path(dir_raw_clim, "current", "current_europe_coast1km_selectedvars.tif")
bio_asia_file <- file.path(dir_raw_clim, "current", "current_asia_coast1km_selectedvars.tif")

if (!file.exists(occ_file)) {
  stop("Missing file: ", occ_file, "\nRun script 03 first.")
}
if (!file.exists(bio_eu_file)) {
  stop("Missing file: ", bio_eu_file, "\nRun script 04 first.")
}
if (!file.exists(bio_asia_file)) {
  stop("Missing file: ", bio_asia_file, "\nRun script 04 first.")
}

# ------------------------------------------------------------
# 4. Read coastal occurrences and climate
# ------------------------------------------------------------
message_section("Read coastal occurrences and CHELSA coastal climate")

occ <- readr::read_csv(occ_file, show_col_types = FALSE)

bio_eu   <- terra::rast(bio_eu_file)
bio_asia <- terra::rast(bio_asia_file)

cat("Europe coastal climate layers:\n")
print(names(bio_eu))
cat("\nAsia coastal climate layers:\n")
print(names(bio_asia))

# ------------------------------------------------------------
# 5. Choose variable set
# ------------------------------------------------------------
message_section("Select variable set for niche analysis")

niche_vars <- intersect(candidate_vars, intersect(names(bio_eu), names(bio_asia)))

if (length(niche_vars) < 2) {
  stop("Too few common climate variables found between Europe and Asia coastal rasters.")
}

cat("Variables used for niche analysis:\n")
print(niche_vars)

readr::write_csv(
  tibble::tibble(variable = niche_vars),
  file.path(dir_tab, "table_niche_variables.csv")
)

# ------------------------------------------------------------
# 6. Coastal native and invaded subsets
# ------------------------------------------------------------
message_section("Prepare coastal native and invaded occurrence subsets")

if (!"range_status" %in% names(occ)) {
  stop("Column 'range_status' not found in coastal occurrences file. Check 05_define_raster_valid_coastal_occurrences_and_background.R output.")
}

native_occ <- occ |>
  dplyr::filter(range_status == "native") |>
  dplyr::select(decimal_longitude, decimal_latitude)

inv_occ <- occ |>
  dplyr::filter(range_status == "invaded") |>
  dplyr::select(decimal_longitude, decimal_latitude)

cat("Native coastal occurrences before thinning:", nrow(native_occ), "\n")
cat("Invaded European coastal occurrences before thinning:", nrow(inv_occ), "\n")

if (nrow(native_occ) < 10) stop("Too few native coastal records.")
if (nrow(inv_occ) < 10) stop("Too few invaded coastal records.")

native_occ <- thin_by_cell_df(native_occ, bio_asia[[1]])
inv_occ    <- thin_by_cell_df(inv_occ, bio_eu[[1]])

cat("Native coastal occurrences after thinning:", nrow(native_occ), "\n")
cat("Invaded European coastal occurrences after thinning:", nrow(inv_occ), "\n")

# ------------------------------------------------------------
# 7. Extract climate values
# ------------------------------------------------------------
message_section("Extract climate values")

native_vals <- terra::extract(
  bio_asia[[niche_vars]],
  native_occ[, c("decimal_longitude", "decimal_latitude")]
)[, -1, drop = FALSE]

inv_vals <- terra::extract(
  bio_eu[[niche_vars]],
  inv_occ[, c("decimal_longitude", "decimal_latitude")]
)[, -1, drop = FALSE]

native_dat <- cbind(native_occ, native_vals) |> stats::na.omit()
inv_dat    <- cbind(inv_occ, inv_vals) |> stats::na.omit()

cat("Native coastal occurrences used after extraction:", nrow(native_dat), "\n")
cat("Invaded European coastal occurrences used after extraction:", nrow(inv_dat), "\n")

if (nrow(native_dat) < 10) stop("Too few native coastal records remain after extraction.")
if (nrow(inv_dat) < 10) stop("Too few invaded coastal records remain after extraction.")

# ------------------------------------------------------------
# 8. Sample coastal environmental background
# ------------------------------------------------------------
message_section("Sample coastal environmental background")

set.seed(123)

bg_asia <- terra::spatSample(
  bio_asia[[niche_vars]],
  size = bg_sample_n,
  method = "random",
  xy = FALSE,
  na.rm = TRUE,
  as.df = TRUE
)

bg_eu <- terra::spatSample(
  bio_eu[[niche_vars]],
  size = bg_sample_n,
  method = "random",
  xy = FALSE,
  na.rm = TRUE,
  as.df = TRUE
)

cat("Asian coastal background size:", nrow(bg_asia), "\n")
cat("European coastal background size:", nrow(bg_eu), "\n")

# ------------------------------------------------------------
# 9. Run PCA-env
# ------------------------------------------------------------
message_section("Run PCA-env")

glob <- rbind(
  bg_asia,
  bg_eu,
  native_dat[, niche_vars, drop = FALSE],
  inv_dat[, niche_vars, drop = FALSE]
)

pca_env <- stats::prcomp(glob, center = TRUE, scale. = TRUE)

native_scores  <- stats::predict(pca_env, newdata = native_dat[, niche_vars, drop = FALSE])
inv_scores     <- stats::predict(pca_env, newdata = inv_dat[, niche_vars, drop = FALSE])
glob_scores    <- stats::predict(pca_env, newdata = glob[, niche_vars, drop = FALSE])
bg_asia_scores <- stats::predict(pca_env, newdata = bg_asia[, niche_vars, drop = FALSE])
bg_eu_scores   <- stats::predict(pca_env, newdata = bg_eu[, niche_vars, drop = FALSE])

cat("Explained variance (PC1, PC2):\n")
print(summary(pca_env)$importance[2, 1:2])

# ------------------------------------------------------------
# 10. ecospat grids and tests
# ------------------------------------------------------------
message_section("Build ecospat grids and run tests")

z_native <- ecospat::ecospat.grid.clim.dyn(
  glob  = as.matrix(glob_scores[, 1:2]),
  glob1 = as.matrix(bg_asia_scores[, 1:2]),
  sp    = as.matrix(native_scores[, 1:2]),
  R = 100
)

z_inv <- ecospat::ecospat.grid.clim.dyn(
  glob  = as.matrix(glob_scores[, 1:2]),
  glob1 = as.matrix(bg_eu_scores[, 1:2]),
  sp    = as.matrix(inv_scores[, 1:2]),
  R = 100
)

overlap <- ecospat::ecospat.niche.overlap(z_native, z_inv, cor = TRUE)
eq_test <- ecospat::ecospat.niche.equivalency.test(z_native, z_inv, rep = 100)
sim_test <- ecospat::ecospat.niche.similarity.test(z_native, z_inv, rep = 100, rand.type = 2)

eq_sim  <- if ("simulated.overlap" %in% names(eq_test)) eq_test$simulated.overlap else NA
sim_sim <- if ("simulated.overlap" %in% names(sim_test)) sim_test$simulated.overlap else NA

cat("Observed Schoener's D:", overlap$D, "\n")
if ("p.D" %in% names(eq_test)) cat("Equivalency test p.D:", eq_test$p.D, "\n")
if ("p.D" %in% names(sim_test)) cat("Similarity test p.D:", sim_test$p.D, "\n")

# ------------------------------------------------------------
# 11. Save results
# ------------------------------------------------------------
message_section("Save results")

niche_summary <- tibble::tibble(
  metric = c(
    "n_native_coastal_raw",
    "n_invaded_coastal_raw",
    "n_native_coastal_used",
    "n_invaded_coastal_used",
    "niche_overlap_D",
    "equivalency_pD",
    "similarity_pD"
  ),
  value = c(
    nrow(native_occ),
    nrow(inv_occ),
    nrow(native_dat),
    nrow(inv_dat),
    overlap$D,
    if ("p.D" %in% names(eq_test)) eq_test$p.D else NA_real_,
    if ("p.D" %in% names(sim_test)) sim_test$p.D else NA_real_
  )
)

readr::write_csv(
  niche_summary,
  file.path(dir_tab, "table_niche_summary.csv")
)

saveRDS(
  list(
    pca_env = pca_env,
    overlap = overlap,
    eq_test = eq_test,
    sim_test = sim_test,
    niche_vars = niche_vars,
    native_dat = native_dat,
    inv_dat = inv_dat,
    bg_asia = bg_asia,
    bg_eu = bg_eu,
    native_scores = native_scores,
    inv_scores = inv_scores
  ),
  file.path(dir_models, "niche_results.rds")
)

# ------------------------------------------------------------
# 12. Editorial figures
# ------------------------------------------------------------
message_section("Build editorial figures")

pca_df <- rbind(
  data.frame(PC1 = native_scores[, 1], PC2 = native_scores[, 2], group = "Native"),
  data.frame(PC1 = inv_scores[, 1],    PC2 = inv_scores[, 2],    group = "Invaded")
)

p1 <- ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2, color = group)) +
  ggplot2::geom_point(alpha = 0.30, size = 0.7) +
  ggplot2::stat_ellipse(linewidth = 0.9, alpha = 0.8) +
  ggplot2::scale_color_manual(values = c(Native = "#1A9641", Invaded = "#D7191C")) +
  ggplot2::labs(
    title = "Coastal climatic niche space of Rosa rugosa",
    subtitle = "Native East Asia versus invaded Europe",
    x = "PC1",
    y = "PC2",
    color = NULL
  ) +
  theme_editorial()

save_editorial_plot(
  p1,
  file.path(dir_fig, "fig05_pca_env_native_invaded.png"),
  width = 7.2,
  height = 6
)

if (!all(is.na(eq_sim))) {
  p2 <- ggplot2::ggplot(tibble::tibble(simulated_D = eq_sim), ggplot2::aes(x = simulated_D)) +
    ggplot2::geom_histogram(bins = 25) +
    ggplot2::geom_vline(xintercept = overlap$D, linewidth = 1) +
    ggplot2::labs(
      title = "Niche equivalency test",
      subtitle = paste("Observed D =", round(overlap$D, 3)),
      x = "Simulated D",
      y = "Frequency"
    ) +
    theme_editorial()
  
  save_editorial_plot(
    p2,
    file.path(dir_fig, "fig06_niche_equivalency_test.png"),
    width = 7.2,
    height = 5
  )
}

if (!all(is.na(sim_sim))) {
  p3 <- ggplot2::ggplot(tibble::tibble(simulated_D = sim_sim), ggplot2::aes(x = simulated_D)) +
    ggplot2::geom_histogram(bins = 25) +
    ggplot2::geom_vline(xintercept = overlap$D, linewidth = 1) +
    ggplot2::labs(
      title = "Niche similarity test",
      subtitle = paste("Observed D =", round(overlap$D, 3)),
      x = "Simulated D",
      y = "Frequency"
    ) +
    theme_editorial()
  
  save_editorial_plot(
    p3,
    file.path(dir_fig, "fig07_niche_similarity_test.png"),
    width = 7.2,
    height = 5
  )
}

message_section("Coastal niche analysis completed")
cat("Saved niche summary table at:\n", file.path(dir_tab, "table_niche_summary.csv"), "\n")
cat("Saved niche results object at:\n", file.path(dir_models, "niche_results.rds"), "\n")
cat("Saved figures in:\n", dir_fig, "\n")