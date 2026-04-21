# 11_refresh_figure_2_climatic_space.R
# Rebuild the final main-text Figure 2 in the paper output directory.

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(readr)
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

project_dir <- getwd()

dir_clean   <- file.path(project_dir, "data_clean")
dir_tab     <- file.path(project_dir, "outputs", "tables")
dir_paper   <- file.path(project_dir, "outputs", "paper_revision1")
out_main_fig <- file.path(dir_paper, "MAIN_TEXT", "FIGURES")
dir.create(out_main_fig, recursive = TRUE, showWarnings = FALSE)

occ_file       <- file.path(dir_clean, "occurrences_native_invaded_coastal.csv")
niche_varsfile <- file.path(dir_tab, "table_niche_variables.csv")
clim_eu_file   <- file.path(project_dir, "data_raw", "climate_chelsa", "current", "current_europe_coast1km_selectedvars.tif")
clim_asia_file <- file.path(project_dir, "data_raw", "climate_chelsa", "current", "current_asia_coast1km_selectedvars.tif")

req_files <- c(occ_file, niche_varsfile, clim_eu_file, clim_asia_file)
missing_files <- req_files[!file.exists(req_files)]
if (length(missing_files) > 0) {
  stop("Missing required files:\n", paste(missing_files, collapse = "\n"))
}

message_section <- function(txt) {
  cat("\n", paste(rep("=", nchar(txt) + 4), collapse = ""), "\n", sep = "")
  cat("  ", txt, "\n", sep = "")
  cat(paste(rep("=", nchar(txt) + 4), collapse = ""), "\n\n", sep = "")
}

save_plot <- function(p, filename, width, height, dpi = 400) {
  ggsave(
    filename = filename,
    plot = p,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
}

theme_editorial_scatter <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = "grey85"),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      axis.title       = element_text(face = "bold"),
      plot.title       = element_text(face = "bold", size = base_size + 2),
      legend.position  = "bottom",
      legend.title     = element_blank()
    )
}

sample_valid_values <- function(r, vars, n_max = 8000) {
  x <- as.data.frame(r[[vars]], na.rm = TRUE)
  if (nrow(x) > n_max) {
    set.seed(123)
    x <- dplyr::slice_sample(x, n = n_max)
  }
  x
}

extract_occ_env <- function(df, r, vars) {
  xy <- df[, c("decimal_longitude", "decimal_latitude")]
  ex <- terra::extract(r[[vars]], xy)[, -1, drop = FALSE]
  bind_cols(df, as_tibble(ex)) |>
    dplyr::filter(dplyr::if_all(dplyr::all_of(vars), ~ !is.na(.)))
}

message_section("Read data")

occ_coast <- readr::read_csv(occ_file, show_col_types = FALSE)
niche_vars <- readr::read_csv(niche_varsfile, show_col_types = FALSE)
clim_eu <- terra::rast(clim_eu_file)
clim_asia <- terra::rast(clim_asia_file)

pca_vars <- niche_vars$variable
pca_vars <- intersect(pca_vars, intersect(names(clim_eu), names(clim_asia)))
if (length(pca_vars) < 3) stop("Too few niche-analysis variables available to rebuild Figure 2.")

message_section("Fit PCA on coastal background climates")

bg_eu_vals   <- sample_valid_values(clim_eu, pca_vars, n_max = 8000)
bg_asia_vals <- sample_valid_values(clim_asia, pca_vars, n_max = 8000)

pca_bg <- bind_rows(
  mutate(bg_asia_vals, region = "Native background"),
  mutate(bg_eu_vals, region = "Invaded background")
)

pca_fit <- prcomp(pca_bg[, pca_vars], center = TRUE, scale. = TRUE)

occ_native_env <- occ_coast |>
  filter(range_status == "native") |>
  extract_occ_env(clim_asia, pca_vars) |>
  mutate(group = "Native")

occ_invaded_env <- occ_coast |>
  filter(range_status == "invaded") |>
  extract_occ_env(clim_eu, pca_vars) |>
  mutate(group = "Invaded")

occ_pca <- bind_rows(occ_native_env, occ_invaded_env)

scores <- predict(pca_fit, newdata = occ_pca[, pca_vars])
occ_scores <- bind_cols(
  occ_pca |> select(group),
  as.data.frame(scores)[, 1:2]
) |>
  rename(PC1 = PC1, PC2 = PC2)

col_native  <- "#1f78b4"
col_invaded <- "#d95f02"

p_fig2 <- ggplot(occ_scores, aes(PC1, PC2, color = group)) +
  geom_point(alpha = 0.35, size = 1.1) +
  stat_ellipse(aes(group = group), linewidth = 1.2, level = 0.90, show.legend = FALSE) +
  scale_color_manual(values = c("Native" = col_native, "Invaded" = col_invaded)) +
  labs(x = "PC1", y = "PC2") +
  theme_editorial_scatter()

save_plot(
  p_fig2,
  file.path(out_main_fig, "Figure_2_climatic_space_PCA.png"),
  width = 7.8,
  height = 6.2
)

message_section("Figure 2 saved")
