# 02_clean_occurrences.R
# Clean GBIF occurrences, retain a flagged archive, and export the modelling dataset used downstream.

source("00_setup_project.R")
s <- readRDS("project_settings.rds")

message_section("Read GBIF raw data")
raw_path <- file.path(s$dir_raw_gbif, "rosa_rugosa_gbif_raw.csv")
if (!file.exists(raw_path)) stop("Run 01_download_gbif_occurrences.R first.")

# Read everything as character first: safer for heterogeneous GBIF tables
occ <- readr::read_csv(
  raw_path,
  col_types = readr::cols(.default = readr::col_character()),
  show_col_types = FALSE
) |>
  janitor::clean_names()

message_section("Inspect available columns")
print(names(occ))

message_section("Keep relevant columns")
keep_cols <- c(
  "gbif_id", "species", "scientific_name", "accepted_scientific_name",
  "decimal_longitude", "decimal_latitude", "country", "country_code",
  "year", "month", "day", "basis_of_record", "dataset_name",
  "coordinate_uncertainty_in_meters", "issues", "institution_code",
  "locality", "occurrence_remarks", "habitat", "establishment_means"
)
keep_cols <- keep_cols[keep_cols %in% names(occ)]
occ <- occ[, keep_cols, drop = FALSE]

message_section("Create safe helper fields")

# Ensure expected text columns always exist
ensure_col <- function(df, colname, default = NA_character_) {
  if (!colname %in% names(df)) df[[colname]] <- default
  df
}

text_cols_needed <- c(
  "species", "scientific_name", "accepted_scientific_name",
  "country", "country_code", "locality", "occurrence_remarks",
  "habitat", "establishment_means", "institution_code"
)

for (cl in text_cols_needed) {
  occ <- ensure_col(occ, cl, NA_character_)
}

# Species label used in downstream steps
occ <- occ |>
  dplyr::mutate(
    species_use = dplyr::coalesce(species, accepted_scientific_name, scientific_name),
    country_use = dplyr::coalesce(country, country_code)
  )

message_section("Convert numeric fields safely")
to_numeric_if_present <- function(df, colname) {
  if (colname %in% names(df)) {
    df[[colname]] <- suppressWarnings(as.numeric(df[[colname]]))
  }
  df
}

occ <- to_numeric_if_present(occ, "decimal_longitude")
occ <- to_numeric_if_present(occ, "decimal_latitude")
occ <- to_numeric_if_present(occ, "year")
occ <- to_numeric_if_present(occ, "month")
occ <- to_numeric_if_present(occ, "day")
occ <- to_numeric_if_present(occ, "coordinate_uncertainty_in_meters")

message_section("Basic filters")
occ <- occ |>
  dplyr::filter(!is.na(decimal_longitude), !is.na(decimal_latitude)) |>
  dplyr::filter(decimal_longitude >= -180, decimal_longitude <= 180) |>
  dplyr::filter(decimal_latitude >= -90, decimal_latitude <= 90)

message_section("Create duplicate keys")
occ <- occ |>
  dplyr::mutate(
    coord_key = paste(round(decimal_longitude, 5), round(decimal_latitude, 5), sep = "_"),
    text_key = paste(species_use, country_use, year, coord_key, sep = "_")
  )

message_section("CoordinateCleaner tests")

# CoordinateCleaner needs a species field name
occ$species_cc <- occ$species_use

cc <- CoordinateCleaner::clean_coordinates(
  x = occ,
  lon = "decimal_longitude",
  lat = "decimal_latitude",
  species = "species_cc",
  tests = c("capitals", "centroids", "equal", "gbif", "institutions", "seas", "zeros"),
  value = "spatialvalid"
)

occ <- cc

print(summary(occ))

flag_cols <- c("cc_cap", "cc_cen", "cc_equ", "cc_gbf", "cc_inst", "cc_sea", "cc_zero")
flag_cols <- flag_cols[flag_cols %in% names(occ)]

if (length(flag_cols) == 0) {
  occ$cc_pass <- TRUE
} else {
  occ$cc_pass <- apply(occ[, flag_cols, drop = FALSE], 1, function(x) all(x))
}

# Additional diagnostic plot from CoordinateCleaner
png(
  filename = file.path(s$dir_fig, "diag05_coordinatecleaner_flags.png"),
  width = 1800, height = 1200, res = 220
)
plot(occ, lon = "decimal_longitude", lat = "decimal_latitude")
dev.off()

message_section("Flag likely cultivated / ornamental records")

cult_terms <- c(
  "garden", "botanic", "botanical", "arboret", "nursery", "park",
  "cultiv", "ornamental", "planted", "hortic", "campus"
)

collapse_text_vec <- function(...) {
  x <- list(...)
  x <- lapply(x, function(v) {
    v[is.na(v)] <- ""
    v
  })
  out <- do.call(paste, c(x, sep = " | "))
  stringr::str_to_lower(out)
}

occ <- occ |>
  dplyr::mutate(
    text_blob = collapse_text_vec(
      locality,
      occurrence_remarks,
      habitat,
      establishment_means,
      institution_code
    ),
    likely_cultivated = stringr::str_detect(
      text_blob,
      paste(cult_terms, collapse = "|")
    )
  )

message_section("Flag duplicate coordinates")
occ <- occ |>
  dplyr::group_by(coord_key) |>
  dplyr::mutate(duplicate_coord = dplyr::n() > 1) |>
  dplyr::ungroup()

message_section("Build cleaned datasets")

occ_flagged <- occ

occ_model <- occ |>
  dplyr::filter(cc_pass) |>
  dplyr::filter(!likely_cultivated) |>
  dplyr::distinct(coord_key, .keep_all = TRUE)

if ("coordinate_uncertainty_in_meters" %in% names(occ_model)) {
  occ_model <- occ_model |>
    dplyr::filter(
      is.na(coordinate_uncertainty_in_meters) |
        coordinate_uncertainty_in_meters <= 10000
    )
}

readr::write_csv(occ_flagged, file.path(s$dir_clean, "occurrences_flagged.csv"))
readr::write_csv(occ_model,   file.path(s$dir_clean, "occurrences_model_clean.csv"))

message_section("Diagnostics")

summary_tbl <- tibble::tibble(
  step = c(
    "Raw after numeric/geographic filter",
    "Pass CoordinateCleaner",
    "Not likely cultivated",
    "Distinct coordinates",
    "Final model dataset"
  ),
  n = c(
    nrow(occ),
    sum(occ$cc_pass, na.rm = TRUE),
    sum(occ$cc_pass & !occ$likely_cultivated, na.rm = TRUE),
    nrow(dplyr::distinct(
      occ |>
        dplyr::filter(cc_pass, !likely_cultivated),
      coord_key,
      .keep_all = TRUE
    )),
    nrow(occ_model)
  )
)

readr::write_csv(summary_tbl, file.path(s$dir_tab, "table_cleaning_summary.csv"))

p1 <- ggplot2::ggplot(summary_tbl, ggplot2::aes(x = step, y = n)) +
  ggplot2::geom_col() +
  ggplot2::coord_flip() +
  ggplot2::labs(
    title = "Occurrence filtering workflow",
    x = NULL,
    y = "Number of records"
  ) +
  theme_editorial()

save_editorial_plot(
  p1,
  file.path(s$dir_fig, "diag03_cleaning_flow.png"),
  width = 8.2,
  height = 4.8
)

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

p2a <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = world, fill = "grey95", color = "grey70", linewidth = 0.2) +
  ggplot2::geom_point(
    data = occ,
    ggplot2::aes(x = decimal_longitude, y = decimal_latitude),
    alpha = 0.22,
    size = 0.32
  ) +
  ggplot2::coord_sf(xlim = c(-20, 180), ylim = c(15, 75), expand = FALSE) +
  ggplot2::labs(
    title = "Records after basic geographic filtering",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_editorial()

p2b <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = world, fill = "grey95", color = "grey70", linewidth = 0.2) +
  ggplot2::geom_point(
    data = occ_model,
    ggplot2::aes(x = decimal_longitude, y = decimal_latitude),
    alpha = 0.30,
    size = 0.32
  ) +
  ggplot2::coord_sf(xlim = c(-20, 180), ylim = c(15, 75), expand = FALSE) +
  ggplot2::labs(
    title = "Cleaned modelling dataset",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_editorial()

p2 <- p2a + p2b + patchwork::plot_layout(ncol = 2)

save_editorial_plot(
  p2,
  file.path(s$dir_fig, "fig02_raw_vs_clean_map.png"),
  width = 11,
  height = 5.5
)

if ("year" %in% names(occ_model)) {
  p3 <- ggplot2::ggplot(occ_model, ggplot2::aes(x = year)) +
    ggplot2::geom_histogram(bins = 35) +
    ggplot2::labs(
      title = "Temporal profile of cleaned records",
      x = "Year",
      y = "Number of records"
    ) +
    theme_editorial()
  
  save_editorial_plot(
    p3,
    file.path(s$dir_fig, "diag04_clean_years.png"),
    width = 8,
    height = 4.5
  )
}

message_section("Cleaning completed")
cat("Flagged dataset saved at:\n", file.path(s$dir_clean, "occurrences_flagged.csv"), "\n")
cat("Modelling dataset saved at:\n", file.path(s$dir_clean, "occurrences_model_clean.csv"), "\n")
cat("Cleaning summary table saved at:\n", file.path(s$dir_tab, "table_cleaning_summary.csv"), "\n")