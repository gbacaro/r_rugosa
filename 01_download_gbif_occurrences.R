# 01_download_gbif_occurrences.R
# Download the GBIF occurrence archive for Rosa rugosa and save the raw table and metadata.

source("00_setup_project.R")
s <- readRDS("project_settings.rds")

message_section("Resolve taxon key from GBIF backbone")
bb <- rgbif::name_backbone(name = s$species_name)
print(bb)

if (is.null(bb$usageKey) || length(bb$usageKey) == 0) {
  stop("GBIF taxon key not found. Please inspect the backbone match manually.")
}
taxon_key <- as.character(bb$usageKey)[1]

# Output paths
raw_csv_path  <- file.path(s$dir_raw_gbif, "rosa_rugosa_gbif_raw.csv")
meta_json_path <- file.path(s$dir_raw_gbif, "gbif_download_metadata.json")

# ------------------------------------------------------------
# Helper: robust extraction of GBIF download key
# ------------------------------------------------------------
extract_download_key <- function(x) {
  if (is.character(x) && length(x) >= 1) {
    return(as.character(x)[1])
  }
  if (is.list(x) && "key" %in% names(x)) {
    return(as.character(x$key)[1])
  }
  y <- tryCatch(as.character(x)[1], error = function(e) NA_character_)
  return(y)
}

# ------------------------------------------------------------
# If raw CSV already exists, do not repeat the download
# ------------------------------------------------------------
if (file.exists(raw_csv_path)) {
  
  message_section("Raw GBIF CSV already exists")
  cat("Existing file found:\n", raw_csv_path, "\n")
  cat("The download step will be skipped.\n\n")
  
  raw_occ <- readr::read_csv(raw_csv_path, show_col_types = FALSE)
  
  # try to recover metadata if already present
  download_key <- NA_character_
  gbif_doi <- NA_character_
  
  if (file.exists(meta_json_path)) {
    old_meta <- tryCatch(jsonlite::read_json(meta_json_path, simplifyVector = TRUE),
                         error = function(e) NULL)
    if (!is.null(old_meta)) {
      if ("gbif_download_key" %in% names(old_meta)) download_key <- old_meta$gbif_download_key
      if ("gbif_download_doi" %in% names(old_meta)) gbif_doi <- old_meta$gbif_download_doi
    }
  }
  
} else {
  
  message_section("Check GBIF credentials")
  
  gbif_user  <- Sys.getenv("GBIF_USER")
  gbif_pwd   <- Sys.getenv("GBIF_PWD")
  gbif_email <- Sys.getenv("GBIF_EMAIL")
  
  if (gbif_user == "" || gbif_pwd == "" || gbif_email == "") {
    stop(
      paste(
        "GBIF credentials not found.",
        "Please set GBIF_USER, GBIF_PWD and GBIF_EMAIL in your .Renviron",
        "or with Sys.setenv() before running this script.",
        sep = "\n"
      )
    )
  }
  
  message_section("Submit GBIF download")
  
  preds <- rgbif::pred_in(
    "basisOfRecord",
    c("HUMAN_OBSERVATION", "OBSERVATION", "PRESERVED_SPECIMEN", "LIVING_SPECIMEN")
  )
  
  dl <- rgbif::occ_download(
    rgbif::pred("taxonKey", taxon_key),
    rgbif::pred("hasCoordinate", TRUE),
    preds,
    format = "SIMPLE_CSV",
    user = gbif_user,
    pwd = gbif_pwd,
    email = gbif_email
  )
  
  print(dl)
  
  download_key <- extract_download_key(dl)
  
  if (is.na(download_key) || download_key == "") {
    stop("Could not extract the GBIF download key from the object returned by occ_download().")
  }
  
  message_section("Wait for completion")
  rgbif::occ_download_wait(download_key)
  
  message_section("Retrieve the download")
  dl_file <- rgbif::occ_download_get(download_key, overwrite = TRUE, path = s$dir_raw_gbif)
  raw_occ <- rgbif::occ_download_import(dl_file)
  
  readr::write_csv(raw_occ, raw_csv_path)
  
  # Try to retrieve DOI / citation
  gbif_doi <- tryCatch({
    cit <- rgbif::gbif_citation(download_key)
    if (is.list(cit) && "download" %in% names(cit) && "doi" %in% names(cit$download)) {
      as.character(cit$download$doi)
    } else if (is.list(cit) && "doi" %in% names(cit)) {
      as.character(cit$doi)
    } else {
      NA_character_
    }
  }, error = function(e) NA_character_)
  
  meta <- list(
    species_name = s$species_name,
    taxon_key = taxon_key,
    gbif_download_key = download_key,
    gbif_download_doi = gbif_doi,
    n_records = nrow(raw_occ),
    download_date = as.character(Sys.time())
  )
  
  jsonlite::write_json(
    meta,
    path = meta_json_path,
    auto_unbox = TRUE,
    pretty = TRUE
  )
}

# ------------------------------------------------------------
# Final metadata check / rewrite
# ------------------------------------------------------------
if (!exists("download_key")) download_key <- NA_character_
if (!exists("gbif_doi")) gbif_doi <- NA_character_

meta <- list(
  species_name = s$species_name,
  taxon_key = taxon_key,
  gbif_download_key = download_key,
  gbif_download_doi = gbif_doi,
  n_records = nrow(raw_occ),
  download_date = as.character(Sys.time())
)

jsonlite::write_json(
  meta,
  path = meta_json_path,
  auto_unbox = TRUE,
  pretty = TRUE
)

message_section("Basic diagnostics")

plot_df <- raw_occ |>
  dplyr::mutate(
    year = suppressWarnings(as.integer(year)),
    basisOfRecord = as.character(basisOfRecord)
  )

# Plot 1: basisOfRecord
p1 <- ggplot2::ggplot(plot_df, ggplot2::aes(x = forcats::fct_infreq(basisOfRecord))) +
  ggplot2::geom_bar() +
  ggplot2::coord_flip() +
  ggplot2::labs(
    title = "GBIF records by basis of record",
    x = NULL,
    y = "Number of records"
  ) +
  theme_editorial()

save_editorial_plot(
  p1,
  file.path(s$dir_fig, "diag01_basis_of_record.png"),
  width = 8,
  height = 4.5
)

# Plot 2: year histogram
p2 <- ggplot2::ggplot(plot_df, ggplot2::aes(x = year)) +
  ggplot2::geom_histogram(bins = 40) +
  ggplot2::labs(
    title = "Temporal distribution of GBIF records",
    x = "Year",
    y = "Number of records"
  ) +
  theme_editorial()

save_editorial_plot(
  p2,
  file.path(s$dir_fig, "diag02_record_years_raw.png"),
  width = 8,
  height = 4.5
)

# Plot 3: raw occurrence map
# NOTE: GBIF imported columns usually retain camelCase names
lon_col <- if ("decimalLongitude" %in% names(plot_df)) "decimalLongitude" else "decimal_longitude"
lat_col <- if ("decimalLatitude" %in% names(plot_df)) "decimalLatitude" else "decimal_latitude"

plot_map <- plot_df |>
  dplyr::filter(!is.na(.data[[lon_col]]), !is.na(.data[[lat_col]]))

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

p3 <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = world, fill = "grey95", color = "grey70", linewidth = 0.2) +
  ggplot2::geom_point(
    data = plot_map,
    ggplot2::aes(x = .data[[lon_col]], y = .data[[lat_col]]),
    alpha = 0.25,
    size = 0.4
  ) +
  ggplot2::coord_sf(xlim = c(-20, 180), ylim = c(15, 75), expand = FALSE) +
  ggplot2::labs(
    title = "Raw GBIF occurrences of Rosa rugosa",
    subtitle = ifelse(is.na(download_key), "Existing local file reused", paste("Download key:", download_key)),
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_editorial()

save_editorial_plot(
  p3,
  file.path(s$dir_fig, "fig01_raw_occurrences_map.png"),
  width = 10,
  height = 5.6
)

message_section("Done")
cat("Raw file saved at:\n", raw_csv_path, "\n")
cat("Metadata saved at:\n", meta_json_path, "\n")
cat("Download key:", download_key, "\n")
cat("Download DOI:", gbif_doi, "\n")