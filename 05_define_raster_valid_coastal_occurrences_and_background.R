# 05_define_raster_valid_coastal_occurrences_and_background.R
# Classify occurrences as native or invaded, snap them to raster-valid coastal cells, and sample the European coastal background.

library(terra)
library(readr)
library(dplyr)
library(sf)
library(ggplot2)
library(tibble)
library(rnaturalearth)
library(rnaturalearthdata)

project_dir <- getwd()

dir_clean <- file.path(project_dir, "data_clean")
dir_fig   <- file.path(project_dir, "outputs", "figures")
dir_tab   <- file.path(project_dir, "outputs", "tables")
dir_raw_clim <- file.path(project_dir, "data_raw", "climate_chelsa")

dir.create(dir_clean, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_tab, recursive = TRUE, showWarnings = FALSE)

eu_ext_vec   <- c(-15, 45, 30, 72)
asia_ext_vec <- c(118, 156, 28, 58)

candidate_vars <- c("bio3", "bio8", "bio9", "bio10", "bio15", "bio18", "bio19", "fcf", "gsl", "gst")
max_snap_distance_m <- 1500
background_n_europe <- 10000
background_exclusion_m <- 1000

native_country_names <- toupper(c(
  "RUSSIA", "RUSSIAN FEDERATION",
  "CHINA", "JAPAN",
  "SOUTH KOREA", "KOREA, REPUBLIC OF",
  "NORTH KOREA", "KOREA, DEMOCRATIC PEOPLE'S REPUBLIC OF"
))
native_country_codes <- c("RU", "CN", "JP", "KR", "KP")

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
      legend.position = "right"
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

ensure_col <- function(df, colname, default = NA_character_) {
  if (!colname %in% names(df)) df[[colname]] <- default
  df
}

bbox_to_sf <- function(ext_vec) {
  sf::st_as_sfc(sf::st_bbox(c(
    xmin = ext_vec[1], xmax = ext_vec[2],
    ymin = ext_vec[3], ymax = ext_vec[4]
  ), crs = 4326))
}

build_coastal_belt <- function(ext_vec, buffer_km) {
  bb <- bbox_to_sf(ext_vec)
  coast <- rnaturalearth::ne_coastline(scale = "medium", returnclass = "sf") |>
    suppressWarnings(sf::st_intersection(bb))
  land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") |>
    suppressWarnings(sf::st_intersection(bb))
  
  coast_m <- sf::st_transform(coast, 3857)
  land_m  <- sf::st_transform(land, 3857)
  coast_buf <- sf::st_buffer(coast_m, dist = buffer_km * 1000)
  
  belt <- suppressWarnings(sf::st_intersection(sf::st_union(land_m), sf::st_union(coast_buf)))
  belt <- sf::st_make_valid(belt)
  belt <- suppressWarnings(sf::st_collection_extract(belt, "POLYGON"))
  sf::st_transform(belt, 4326)
}

build_complete_mask <- function(r, vars) {
  vars <- intersect(vars, names(r))
  if (length(vars) == 0) stop("No candidate variables found in raster stack.")
  
  mask <- !is.na(r[[vars[1]]])
  if (length(vars) > 1) {
    for (v in vars[-1]) {
      mask <- mask & !is.na(r[[v]])
    }
  }
  terra::ifel(mask, 1, NA)
}

snap_occurrences_to_mask <- function(df, mask_rast, max_snap_m = 1500, region_label = "region") {
  if (nrow(df) == 0) {
    return(list(
      data = df,
      diagnostics = tibble::tibble(
        region = region_label,
        n_input = 0,
        n_exact_valid = 0,
        n_nearest_snapped = 0,
        n_retained = 0,
        n_unique_cells = 0,
        mean_snap_distance_m = NA_real_,
        max_snap_distance_m = NA_real_
      )
    ))
  }
  
  stopifnot(all(c("decimal_longitude", "decimal_latitude") %in% names(df)))
  
  df <- as.data.frame(df)
  df$row_id <- seq_len(nrow(df))
  df$orig_decimal_longitude <- df$decimal_longitude
  df$orig_decimal_latitude  <- df$decimal_latitude
  
  cell_guess <- terra::cellFromXY(mask_rast, as.matrix(df[, c("decimal_longitude", "decimal_latitude")]))
  mask_vals  <- terra::values(mask_rast, mat = FALSE)
  exact_ok   <- !is.na(cell_guess) & !is.na(mask_vals[cell_guess])
  
  out_list <- list()
  
  if (any(exact_ok)) {
    exact_df <- df[exact_ok, , drop = FALSE]
    exact_df$cell_id <- cell_guess[exact_ok]
    exact_xy <- terra::xyFromCell(mask_rast, exact_df$cell_id)
    exact_df$decimal_longitude <- exact_xy[, 1]
    exact_df$decimal_latitude  <- exact_xy[, 2]
    exact_df$snap_distance_m   <- 0
    exact_df$snap_type         <- "exact_valid_cell"
    out_list[["exact"]] <- exact_df
  }
  
  remaining <- df[!exact_ok, , drop = FALSE]
  
  if (nrow(remaining) > 0) {
    valid_df <- terra::as.data.frame(mask_rast, xy = TRUE, cells = TRUE, na.rm = TRUE)
    valid_df <- valid_df[, c("cell", "x", "y")]
    
    valid_sf <- sf::st_as_sf(valid_df, coords = c("x", "y"), crs = 4326)
    rem_sf   <- sf::st_as_sf(remaining, coords = c("decimal_longitude", "decimal_latitude"), crs = 4326, remove = FALSE)
    
    valid_m <- sf::st_transform(valid_sf, 3857)
    rem_m   <- sf::st_transform(rem_sf, 3857)
    
    idx <- sf::st_nearest_feature(rem_m, valid_m)
    d   <- as.numeric(sf::st_distance(rem_m, valid_m[idx, ], by_element = TRUE))
    
    snapped_df <- remaining
    snapped_df$cell_id <- valid_df$cell[idx]
    snapped_df$decimal_longitude <- valid_df$x[idx]
    snapped_df$decimal_latitude  <- valid_df$y[idx]
    snapped_df$snap_distance_m   <- d
    snapped_df$snap_type         <- "nearest_valid_cell"
    snapped_df <- snapped_df[is.finite(snapped_df$snap_distance_m) & snapped_df$snap_distance_m <= max_snap_m, , drop = FALSE]
    
    if (nrow(snapped_df) > 0) {
      out_list[["nearest"]] <- snapped_df
    }
  }
  
  if (length(out_list) == 0) {
    snapped_all <- df[0, , drop = FALSE]
  } else {
    snapped_all <- dplyr::bind_rows(out_list)
  }
  
  if (nrow(snapped_all) > 0) {
    snapped_all <- snapped_all |>
      dplyr::arrange(cell_id, snap_distance_m) |>
      dplyr::group_by(cell_id) |>
      dplyr::slice(1) |>
      dplyr::ungroup()
  }
  
  diag_tbl <- tibble::tibble(
    region = region_label,
    n_input = nrow(df),
    n_exact_valid = sum(exact_ok, na.rm = TRUE),
    n_nearest_snapped = if (nrow(snapped_all) > 0) sum(snapped_all$snap_type == "nearest_valid_cell", na.rm = TRUE) else 0,
    n_retained = nrow(snapped_all),
    n_unique_cells = dplyr::n_distinct(snapped_all$cell_id),
    mean_snap_distance_m = if (nrow(snapped_all) > 0) mean(snapped_all$snap_distance_m, na.rm = TRUE) else NA_real_,
    max_snap_distance_m = if (nrow(snapped_all) > 0) max(snapped_all$snap_distance_m, na.rm = TRUE) else NA_real_
  )
  
  snapped_all$row_id <- NULL
  
  list(data = snapped_all, diagnostics = diag_tbl)
}

sample_background_from_mask <- function(mask_rast, exclusion_df, size = 10000, exclusion_m = 1000) {
  if (nrow(exclusion_df) == 0) stop("No invaded occurrences available for background exclusion.")
  
  occ_vect <- terra::vect(exclusion_df, geom = c("decimal_longitude", "decimal_latitude"), crs = "EPSG:4326")
  occ_rast <- terra::rasterize(occ_vect, mask_rast, field = 1, background = NA)
  dist_r   <- terra::distance(occ_rast)
  
  bg_mask <- terra::ifel(mask_rast == 1 & dist_r > exclusion_m, 1, NA)
  
  bg_pts <- terra::spatSample(
    bg_mask,
    size = size,
    method = "random",
    xy = TRUE,
    na.rm = TRUE,
    as.df = TRUE
  )
  
  if (nrow(bg_pts) == 0) stop("Background sampling returned zero points.")
  
  bg_df <- tibble::tibble(
    pa = 0,
    decimal_longitude = bg_pts$x,
    decimal_latitude = bg_pts$y,
    cell_id = terra::cellFromXY(mask_rast, as.matrix(bg_pts[, c("x", "y")]))
  ) |>
    dplyr::distinct(cell_id, .keep_all = TRUE)
  
  list(points = bg_df, mask = bg_mask)
}

# ------------------------------------------------------------
# 3. Input files
# ------------------------------------------------------------
occ_file <- file.path(dir_clean, "occurrences_model_clean.csv")
bio_eu_file   <- file.path(dir_raw_clim, "current", "current_europe_coast1km_selectedvars.tif")
bio_asia_file <- file.path(dir_raw_clim, "current", "current_asia_coast1km_selectedvars.tif")

if (!file.exists(occ_file)) stop("Missing file: ", occ_file, "\nRun 02_clean_occurrences.R first.")
if (!file.exists(bio_eu_file)) stop("Missing file: ", bio_eu_file, "\nRun 04_download_chelsa_climate_coastal.R first.")
if (!file.exists(bio_asia_file)) stop("Missing file: ", bio_asia_file, "\nRun 04_download_chelsa_climate_coastal.R first.")

# ------------------------------------------------------------
# 4. Read cleaned occurrences and climate rasters
# ------------------------------------------------------------
message_section("Read cleaned occurrences and current coastal CHELSA rasters")

occ <- readr::read_csv(occ_file, show_col_types = FALSE)
bio_eu   <- terra::rast(bio_eu_file)
bio_asia <- terra::rast(bio_asia_file)

for (cl in c("country", "country_code", "country_use")) occ <- ensure_col(occ, cl, NA_character_)

cat("Europe climate layers available:\n")
print(names(bio_eu))
cat("\nAsia climate layers available:\n")
print(names(bio_asia))

# ------------------------------------------------------------
# 5. Classify native vs invaded
# ------------------------------------------------------------
message_section("Classify records as native vs invaded")

occ <- occ |>
  dplyr::mutate(
    country_name_std = toupper(trimws(dplyr::coalesce(country, country_use, ""))),
    country_code_std = toupper(trimws(dplyr::coalesce(country_code, ""))),
    in_native_bbox =
      decimal_longitude >= asia_ext_vec[1] & decimal_longitude <= asia_ext_vec[2] &
      decimal_latitude  >= asia_ext_vec[3] & decimal_latitude  <= asia_ext_vec[4],
    matches_native_country = country_name_std %in% native_country_names | country_code_std %in% native_country_codes,
    missing_country_info =
      (is.na(country_name_std) | country_name_std == "") &
      (is.na(country_code_std) | country_code_std == ""),
    range_status = dplyr::case_when(
      in_native_bbox & (matches_native_country | missing_country_info) ~ "native",
      TRUE ~ "invaded"
    )
  )

print(table(occ$range_status, useNA = "ifany"))

readr::write_csv(
  occ |>
    dplyr::filter(in_native_bbox) |>
    dplyr::count(country_name_std, country_code_std, sort = TRUE),
  file.path(dir_tab, "table_native_bbox_country_diagnostics.csv")
)

# ------------------------------------------------------------
# 6. Build valid coastal climate masks
# ------------------------------------------------------------
message_section("Build complete coastal climate masks")

valid_eu_mask   <- build_complete_mask(bio_eu, candidate_vars)
valid_asia_mask <- build_complete_mask(bio_asia, candidate_vars)

saveRDS(valid_eu_mask,   file.path(dir_clean, "europe_valid_coastal_cells_1km.rds"))
saveRDS(valid_asia_mask, file.path(dir_clean, "asia_valid_coastal_cells_1km.rds"))

cat("Valid European coastal cells:", terra::global(valid_eu_mask, "sum", na.rm = TRUE)[1, 1], "\n")
cat("Valid East Asian coastal cells:", terra::global(valid_asia_mask, "sum", na.rm = TRUE)[1, 1], "\n")

# ------------------------------------------------------------
# 7. Snap occurrences to valid coastal cells
# ------------------------------------------------------------
message_section("Link occurrences to valid coastal climate cells")

native_occ_raw <- occ |>
  dplyr::filter(range_status == "native") |>
  dplyr::select(-country_name_std, -country_code_std, -in_native_bbox, -matches_native_country, -missing_country_info)

inv_occ_raw <- occ |>
  dplyr::filter(range_status == "invaded") |>
  dplyr::filter(
    decimal_longitude >= eu_ext_vec[1], decimal_longitude <= eu_ext_vec[2],
    decimal_latitude  >= eu_ext_vec[3], decimal_latitude  <= eu_ext_vec[4]
  ) |>
  dplyr::select(-country_name_std, -country_code_std, -in_native_bbox, -matches_native_country, -missing_country_info)

native_snap <- snap_occurrences_to_mask(native_occ_raw, valid_asia_mask, max_snap_m = max_snap_distance_m, region_label = "native_asia")
inv_snap    <- snap_occurrences_to_mask(inv_occ_raw,    valid_eu_mask,   max_snap_m = max_snap_distance_m, region_label = "invaded_europe")

native_occ_coast <- native_snap$data
inv_occ_coast    <- inv_snap$data

diag_tbl <- dplyr::bind_rows(native_snap$diagnostics, inv_snap$diagnostics)
readr::write_csv(diag_tbl, file.path(dir_tab, "table_coastal_cellvalid_diagnostics.csv"))

cat("Native records retained on valid coastal cells:", nrow(native_occ_coast), "\n")
cat("Invaded European records retained on valid coastal cells:", nrow(inv_occ_coast), "\n")

occ_coastal <- dplyr::bind_rows(native_occ_coast, inv_occ_coast)

readr::write_csv(occ,         file.path(dir_clean, "occurrences_native_invaded_all.csv"))
readr::write_csv(occ_coastal, file.path(dir_clean, "occurrences_native_invaded_coastal.csv"))

# ------------------------------------------------------------
# 8. Generate background directly from valid coastal cells
# ------------------------------------------------------------
message_section("Generate European coastal background from valid cells")

bg_obj <- sample_background_from_mask(
  mask_rast = valid_eu_mask,
  exclusion_df = inv_occ_coast[, c("decimal_longitude", "decimal_latitude")],
  size = background_n_europe,
  exclusion_m = background_exclusion_m
)

bg_df <- bg_obj$points
bg_mask <- bg_obj$mask

readr::write_csv(bg_df, file.path(dir_clean, "background_points_europe_coastal.csv"))
saveRDS(bg_mask, file.path(dir_clean, "background_area_europe_coastal_1km.rds"))

cat("European coastal background points sampled from valid cells:", nrow(bg_df), "\n")

# ------------------------------------------------------------
# 9. Editorial figures
# ------------------------------------------------------------
message_section("Build editorial figures")

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
eu_belt_1km   <- build_coastal_belt(eu_ext_vec, 1)
asia_belt_1km <- build_coastal_belt(asia_ext_vec, 1)

p1 <- ggplot() +
  geom_sf(data = world, fill = "grey96", color = "grey75", linewidth = 0.2) +
  geom_sf(data = eu_belt_1km, fill = NA, color = "#2C7FB8", linewidth = 0.5) +
  geom_sf(data = asia_belt_1km, fill = NA, color = "#1A9641", linewidth = 0.5) +
  geom_point(data = occ_coastal, aes(x = decimal_longitude, y = decimal_latitude, color = range_status), alpha = 0.45, size = 0.45) +
  coord_sf(xlim = c(-20, 180), ylim = c(15, 75), expand = FALSE) +
  scale_color_manual(values = c(native = "#1A9641", invaded = "#D7191C")) +
  labs(
    title = "Coastal occurrences of Rosa rugosa",
    subtitle = "Occurrences linked to valid 1 km coastal climate cells in Europe and East Asia",
    x = "Longitude", y = "Latitude", color = NULL
  ) +
  theme_editorial()

save_editorial_plot(
  p1,
  file.path(dir_fig, "fig03_coastal_occurrences_native_invaded.png"),
  width = 10.5,
  height = 5.8
)

p2 <- ggplot() +
  geom_point(data = bg_df, aes(x = decimal_longitude, y = decimal_latitude), color = "#2C7FB8", alpha = 0.25, size = 0.15) +
  geom_point(data = inv_occ_coast, aes(x = decimal_longitude, y = decimal_latitude), color = "#D7191C", alpha = 0.5, size = 0.3) +
  coord_cartesian(xlim = c(-15, 45), ylim = c(30, 72), expand = FALSE) +
  labs(
    title = "European coastal background design",
    subtitle = "Valid coastal climate cells sampled with 1 km exclusion around invaded occurrences",
    x = "Longitude", y = "Latitude"
  ) +
  theme_editorial()

save_editorial_plot(
  p2,
  file.path(dir_fig, "fig04_coastal_background_design.png"),
  width = 8.8,
  height = 6
)

message_section("Cell-valid coastal occurrence and background preparation completed")
cat("Saved updated coastal occurrence files, background points, and diagnostics.\n")