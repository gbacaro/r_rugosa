# ============================================================
# ROSA RUGOSA - REVISION 2
# 01_occurrence_audit_REV2.R
#
# Reconstruct occurrence filtering from the REV1 flagged archive
# using the criteria stated in the final REV1 manuscript, then
# compare the resulting records with the Ribotta/Sara dataset.
#
# IMPORTANT:
# - REV1 files are read only.
# - "likely_cultivated" is audited but NOT used as an exclusion.
# - coordinate uncertainty limit is 5,000 m.
# - this script does NOT yet refit models.
# ============================================================

rm(list = ls())
gc()

# ------------------------------------------------------------
# 1. Paths and packages
# ------------------------------------------------------------

paths_file <- file.path(getwd(), "02_scripts", "REV2_paths.R")
if (!file.exists(paths_file)) {
  stop("Missing 02_scripts/REV2_paths.R. Run 00_setup_REV2.R first.")
}
source(paths_file)

required_pkgs <- c("sf", "dplyr", "readr", "tibble")
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs) > 0) {
  stop("Missing required R packages: ", paste(missing_pkgs, collapse = ", "))
}

dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. Load REV1 workspace into a private environment
# ------------------------------------------------------------

if (!file.exists(rev1_rdata_general)) {
  stop("REV1 RData not found: ", rev1_rdata_general)
}

rev1 <- new.env(parent = baseenv())
load(rev1_rdata_general, envir = rev1)

get_rev1 <- function(nm, required = TRUE) {
  if (!exists(nm, envir = rev1, inherits = FALSE)) {
    if (required) stop("Required REV1 object not found: ", nm)
    return(NULL)
  }
  get(nm, envir = rev1, inherits = FALSE)
}

flagged  <- as.data.frame(get_rev1("flagged"))
old_clean <- as.data.frame(get_rev1("clean"))
old_coastal <- as.data.frame(get_rev1("coastal"))
main_occ <- as.data.frame(get_rev1("main_occ"))
final_invaded <- as.data.frame(get_rev1("occ"))

# ------------------------------------------------------------
# 3. Basic field checks
# ------------------------------------------------------------

required_flagged <- c(
  "decimal_longitude",
  "decimal_latitude",
  "basis_of_record",
  "cc_pass"
)

missing_flagged <- setdiff(required_flagged, names(flagged))
if (length(missing_flagged) > 0) {
  stop(
    "The flagged REV1 archive lacks required fields: ",
    paste(missing_flagged, collapse = ", ")
  )
}

if (!"coordinate_uncertainty_in_meters" %in% names(flagged)) {
  flagged$coordinate_uncertainty_in_meters <- NA_real_
}

if (!"likely_cultivated" %in% names(flagged)) {
  flagged$likely_cultivated <- NA
}

if (!"gbif_id" %in% names(flagged)) {
  stop("The flagged REV1 archive lacks gbif_id; record-level reconciliation is not possible.")
}

# Normalize core types
flagged$decimal_longitude <- suppressWarnings(as.numeric(flagged$decimal_longitude))
flagged$decimal_latitude <- suppressWarnings(as.numeric(flagged$decimal_latitude))
flagged$coordinate_uncertainty_in_meters <- suppressWarnings(
  as.numeric(flagged$coordinate_uncertainty_in_meters)
)
flagged$basis_of_record <- toupper(trimws(as.character(flagged$basis_of_record)))
flagged$gbif_id <- as.character(flagged$gbif_id)

for (nm in c("old_clean", "old_coastal", "main_occ", "final_invaded")) {
  x <- get(nm)
  if ("gbif_id" %in% names(x)) x$gbif_id <- as.character(x$gbif_id)
  assign(nm, x)
}

# ------------------------------------------------------------
# 4. Rebuild the pre-coastal occurrence archive according to
#    the criteria explicitly stated in the final REV1 Methods
# ------------------------------------------------------------

allowed_basis <- c(
  "HUMAN_OBSERVATION",
  "OBSERVATION",
  "PRESERVED_SPECIMEN"
)

flagged_audit <- flagged |>
  dplyr::mutate(
    pass_valid_coordinates =
      !is.na(decimal_longitude) &
      !is.na(decimal_latitude) &
      decimal_longitude >= -180 &
      decimal_longitude <= 180 &
      decimal_latitude >= -90 &
      decimal_latitude <= 90,

    pass_basis =
      basis_of_record %in% allowed_basis,

    pass_coordinate_cleaner =
      !is.na(cc_pass) & cc_pass,

    pass_uncertainty_5km =
      is.na(coordinate_uncertainty_in_meters) |
      coordinate_uncertainty_in_meters <= 5000,

    pass_declared_record_filters =
      pass_valid_coordinates &
      pass_basis &
      pass_coordinate_cleaner &
      pass_uncertainty_5km,

    coord_key_rev2 = paste(
      round(decimal_longitude, 5),
      round(decimal_latitude, 5),
      sep = "_"
    )
  )

declared_pre_dedup <- flagged_audit |>
  dplyr::filter(pass_declared_record_filters)

# Duplicate handling: retain one record per identical rounded coordinate.
# Missing uncertainty is not penalized; records with known smaller uncertainty
# are preferred only when duplicate coordinates need to be collapsed.
declared_precoast <- declared_pre_dedup |>
  dplyr::mutate(
    uncertainty_sort = dplyr::if_else(
      is.na(coordinate_uncertainty_in_meters),
      Inf,
      coordinate_uncertainty_in_meters
    )
  ) |>
  dplyr::arrange(coord_key_rev2, uncertainty_sort, gbif_id) |>
  dplyr::distinct(coord_key_rev2, .keep_all = TRUE) |>
  dplyr::select(-uncertainty_sort)

readr::write_csv(
  declared_precoast,
  file.path(work_dir, "occurrences_rebuilt_precoastal_declared_filters.csv")
)

# ------------------------------------------------------------
# 5. Audit the hidden cultivated filter and the 10-km -> 5-km
#    discrepancy in the earlier cleaned pipeline
# ------------------------------------------------------------

old_clean_ids <- if ("gbif_id" %in% names(old_clean)) unique(old_clean$gbif_id) else character()
old_coastal_ids <- if ("gbif_id" %in% names(old_coastal)) unique(old_coastal$gbif_id) else character()
main_occ_ids <- if ("gbif_id" %in% names(main_occ)) unique(main_occ$gbif_id) else character()
final_invaded_ids <- if ("gbif_id" %in% names(final_invaded)) unique(final_invaded$gbif_id) else character()

flagged_audit <- flagged_audit |>
  dplyr::mutate(
    in_old_clean = gbif_id %in% old_clean_ids,
    in_old_coastal = gbif_id %in% old_coastal_ids,
    in_corrected_all_range = gbif_id %in% main_occ_ids,
    in_final_invaded = gbif_id %in% final_invaded_ids
  )

recovered_due_to_rebuild <- flagged_audit |>
  dplyr::filter(pass_declared_record_filters, !in_old_clean)

cultivated_recovery <- recovered_due_to_rebuild |>
  dplyr::filter(!is.na(likely_cultivated) & likely_cultivated)

old_clean_over_5km <- old_clean
if ("coordinate_uncertainty_in_meters" %in% names(old_clean_over_5km)) {
  old_clean_over_5km$coordinate_uncertainty_in_meters <- suppressWarnings(
    as.numeric(old_clean_over_5km$coordinate_uncertainty_in_meters)
  )
  old_clean_over_5km <- old_clean_over_5km |>
    dplyr::filter(
      !is.na(coordinate_uncertainty_in_meters),
      coordinate_uncertainty_in_meters > 5000
    )
} else {
  old_clean_over_5km <- old_clean_over_5km[0, , drop = FALSE]
}

readr::write_csv(
  recovered_due_to_rebuild,
  file.path(audit_dir, "01A_records_passing_declared_filters_but_absent_old_clean.csv")
)

readr::write_csv(
  cultivated_recovery,
  file.path(audit_dir, "01B_recovered_records_flagged_likely_cultivated.csv")
)

readr::write_csv(
  old_clean_over_5km,
  file.path(audit_dir, "01C_old_clean_records_uncertainty_gt5km.csv")
)

# ------------------------------------------------------------
# 6. Country summaries using country_code, never country names
# ------------------------------------------------------------

normalize_country_code <- function(df) {
  if ("country_code" %in% names(df)) {
    toupper(trimws(as.character(df$country_code)))
  } else {
    rep(NA_character_, nrow(df))
  }
}

country_summary <- function(df, dataset_name) {
  cc <- normalize_country_code(df)
  tibble::tibble(country_code = cc) |>
    dplyr::filter(!is.na(country_code), nzchar(country_code)) |>
    dplyr::count(country_code, name = "n", sort = TRUE) |>
    dplyr::mutate(dataset = dataset_name, .before = 1)
}

country_tables <- dplyr::bind_rows(
  country_summary(flagged, "REV1_flagged"),
  country_summary(old_clean, "REV1_old_clean"),
  country_summary(old_coastal, "REV1_old_coastal"),
  country_summary(main_occ, "REV1_corrected_all_range"),
  country_summary(final_invaded, "REV1_final_invaded"),
  country_summary(declared_precoast, "REV2_rebuilt_precoastal")
)

readr::write_csv(
  country_tables,
  file.path(audit_dir, "01D_country_counts_all_stages.csv")
)

target_codes <- c(
  "GB", "FR", "IT", "ES", "PT", "SI", "HR", "BA", "ME", "AL",
  "GR", "TR", "CY", "MT", "BG", "RO"
)

target_country_summary <- tidyr::complete(
  country_tables |>
    dplyr::filter(country_code %in% target_codes),
  dataset,
  country_code = target_codes,
  fill = list(n = 0)
) |>
  dplyr::arrange(dataset, country_code)

readr::write_csv(
  target_country_summary,
  file.path(audit_dir, "01E_target_country_counts.csv")
)

# ------------------------------------------------------------
# 7. Read Sara/Ribotta shapefile and validate actual geometry
# ------------------------------------------------------------

sara_shp <- list.files(
  sara_dir,
  pattern = "\\.shp$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

if (length(sara_shp) != 1) {
  stop("Expected exactly one Sara/Ribotta shapefile; found ", length(sara_shp))
}

sara_raw <- sf::st_read(sara_shp, quiet = TRUE)
sara_raw$sara_feature_id <- seq_len(nrow(sara_raw))

if (is.na(sf::st_crs(sara_raw))) {
  stop("Sara/Ribotta shapefile has no declared CRS. No CRS will be guessed.")
}

# MULTIPOINT is allowed, but comparison is point-based.
sara_pts <- suppressWarnings(sf::st_cast(sara_raw, "POINT"))

geom_xy <- sf::st_coordinates(sara_pts)

declared_lonlat <- isTRUE(sf::st_is_longlat(sara_pts))
geometry_lonlat_valid <- all(
  is.finite(geom_xy[, "X"]) &
  is.finite(geom_xy[, "Y"]) &
  geom_xy[, "X"] >= -180 &
  geom_xy[, "X"] <= 180 &
  geom_xy[, "Y"] >= -90 &
  geom_xy[, "Y"] <= 90
)

crs_diag <- data.frame(
  declared_crs = sf::st_crs(sara_pts)$input,
  st_is_longlat = declared_lonlat,
  geometry_x_min = min(geom_xy[, "X"], na.rm = TRUE),
  geometry_x_max = max(geom_xy[, "X"], na.rm = TRUE),
  geometry_y_min = min(geom_xy[, "Y"], na.rm = TRUE),
  geometry_y_max = max(geom_xy[, "Y"], na.rm = TRUE),
  geometry_compatible_with_lonlat = geometry_lonlat_valid,
  n_features_source = nrow(sara_raw),
  n_point_geometries = nrow(sara_pts),
  stringsAsFactors = FALSE
)

readr::write_csv(
  crs_diag,
  file.path(audit_dir, "01F_sara_geometry_crs_diagnostic.csv")
)

if (declared_lonlat && !geometry_lonlat_valid) {
  stop(
    "Sara shapefile declares longitude/latitude but its GEOMETRY coordinates are outside valid lon/lat bounds. ",
    "See 03_audit/01F_sara_geometry_crs_diagnostic.csv. CRS was not guessed or changed."
  )
}

# Transform the actual geometry to WGS84 for comparison.
sara_wgs84 <- sf::st_transform(sara_pts, 4326)

# Add country code from Natural Earth if available.
if (requireNamespace("rnaturalearth", quietly = TRUE)) {
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") |>
    sf::st_transform(4326)

  code_field <- if ("iso_a2_eh" %in% names(world)) {
    "iso_a2_eh"
  } else if ("iso_a2" %in% names(world)) {
    "iso_a2"
  } else {
    NA_character_
  }

  if (!is.na(code_field)) {
    country_join <- sf::st_join(
      sara_wgs84,
      world[, code_field, drop = FALSE],
      left = TRUE,
      join = sf::st_within
    )
    sara_wgs84$sara_country_code <- toupper(
      as.character(country_join[[code_field]])
    )
  } else {
    sara_wgs84$sara_country_code <- NA_character_
  }
} else {
  sara_wgs84$sara_country_code <- NA_character_
}

# ------------------------------------------------------------
# 8. Nearest-neighbour reconciliation with REV1 records
# ------------------------------------------------------------

make_points_sf <- function(df, label) {
  if (!all(c("decimal_longitude", "decimal_latitude") %in% names(df))) {
    stop(label, " lacks decimal coordinate fields.")
  }

  x <- df |>
    dplyr::filter(
      !is.na(decimal_longitude),
      !is.na(decimal_latitude),
      decimal_longitude >= -30,
      decimal_longitude <= 60,
      decimal_latitude >= 20,
      decimal_latitude <= 75
    )

  sf::st_as_sf(
    x,
    coords = c("decimal_longitude", "decimal_latitude"),
    crs = 4326,
    remove = FALSE
  )
}

flagged_eu_sf <- make_points_sf(flagged_audit, "flagged")
final_eu_sf <- make_points_sf(final_invaded, "final_invaded")

sara_m <- sf::st_transform(sara_wgs84, 3035)
flagged_m <- sf::st_transform(flagged_eu_sf, 3035)
final_m <- sf::st_transform(final_eu_sf, 3035)

nearest_flagged_idx <- sf::st_nearest_feature(sara_m, flagged_m)
nearest_flagged_dist <- as.numeric(
  sf::st_distance(
    sara_m,
    flagged_m[nearest_flagged_idx, ],
    by_element = TRUE
  )
)

nearest_final_idx <- sf::st_nearest_feature(sara_m, final_m)
nearest_final_dist <- as.numeric(
  sf::st_distance(
    sara_m,
    final_m[nearest_final_idx, ],
    by_element = TRUE
  )
)

nearest_flagged <- sf::st_drop_geometry(
  flagged_eu_sf[nearest_flagged_idx, , drop = FALSE]
)

nearest_final <- sf::st_drop_geometry(
  final_eu_sf[nearest_final_idx, , drop = FALSE]
)

sara_xy <- sf::st_coordinates(sara_wgs84)
sara_attr <- sf::st_drop_geometry(sara_wgs84)

audit_tbl <- tibble::tibble(
  sara_feature_id = sara_attr$sara_feature_id,
  sara_longitude = sara_xy[, "X"],
  sara_latitude = sara_xy[, "Y"],
  sara_country_code = sara_attr$sara_country_code,
  nearest_flagged_distance_m = nearest_flagged_dist,
  nearest_flagged_gbif_id = as.character(nearest_flagged$gbif_id),
  nearest_flagged_country_code = if ("country_code" %in% names(nearest_flagged)) {
    as.character(nearest_flagged$country_code)
  } else {
    NA_character_
  },
  nearest_flagged_basis = if ("basis_of_record" %in% names(nearest_flagged)) {
    as.character(nearest_flagged$basis_of_record)
  } else {
    NA_character_
  },
  nearest_flagged_uncertainty_m =
    if ("coordinate_uncertainty_in_meters" %in% names(nearest_flagged)) {
      suppressWarnings(as.numeric(nearest_flagged$coordinate_uncertainty_in_meters))
    } else {
      NA_real_
    },
  nearest_flagged_cc_pass =
    if ("cc_pass" %in% names(nearest_flagged)) as.logical(nearest_flagged$cc_pass) else NA,
  nearest_flagged_likely_cultivated =
    if ("likely_cultivated" %in% names(nearest_flagged)) {
      as.logical(nearest_flagged$likely_cultivated)
    } else {
      NA
    },
  nearest_flagged_passes_declared_filters =
    as.character(nearest_flagged$gbif_id) %in% declared_precoast$gbif_id,
  nearest_flagged_in_old_clean =
    as.character(nearest_flagged$gbif_id) %in% old_clean_ids,
  nearest_flagged_in_old_coastal =
    as.character(nearest_flagged$gbif_id) %in% old_coastal_ids,
  nearest_flagged_in_corrected_all_range =
    as.character(nearest_flagged$gbif_id) %in% main_occ_ids,
  nearest_flagged_in_final_invaded =
    as.character(nearest_flagged$gbif_id) %in% final_invaded_ids,
  flagged_match_within_1km = nearest_flagged_dist <= 1000,
  flagged_match_within_5km = nearest_flagged_dist <= 5000,
  nearest_final_distance_m = nearest_final_dist,
  nearest_final_gbif_id = as.character(nearest_final$gbif_id),
  represented_in_final_within_1km = nearest_final_dist <= 1000,
  represented_in_final_within_5km = nearest_final_dist <= 5000
)

readr::write_csv(
  audit_tbl,
  file.path(audit_dir, "01G_sara_vs_REV1_record_level_audit.csv")
)

# Sara country summary
sara_country_summary <- audit_tbl |>
  dplyr::count(sara_country_code, name = "n_sara", sort = TRUE)

readr::write_csv(
  sara_country_summary,
  file.path(audit_dir, "01H_sara_country_counts.csv")
)

# ------------------------------------------------------------
# 9. Summary evidence table
# ------------------------------------------------------------

summary_tbl <- tibble::tibble(
  metric = c(
    "REV1 flagged records",
    "REV1 old cleaned records",
    "REV1 old coastal records",
    "REV1 corrected all-range coastal records",
    "REV1 final invaded coastal records",
    "Records passing final declared pre-coastal filters",
    "Records passing final declared filters but absent from old clean",
    "Of recovered records: flagged likely_cultivated",
    "Old clean records with uncertainty > 5 km",
    "Sara source features",
    "Sara point geometries after cast",
    "Sara points with flagged record within 1 km",
    "Sara points with flagged record within 5 km",
    "Sara points represented in final invaded set within 1 km",
    "Sara points represented in final invaded set within 5 km"
  ),
  value = c(
    nrow(flagged),
    nrow(old_clean),
    nrow(old_coastal),
    nrow(main_occ),
    nrow(final_invaded),
    nrow(declared_precoast),
    nrow(recovered_due_to_rebuild),
    nrow(cultivated_recovery),
    nrow(old_clean_over_5km),
    nrow(sara_raw),
    nrow(sara_wgs84),
    sum(audit_tbl$flagged_match_within_1km, na.rm = TRUE),
    sum(audit_tbl$flagged_match_within_5km, na.rm = TRUE),
    sum(audit_tbl$represented_in_final_within_1km, na.rm = TRUE),
    sum(audit_tbl$represented_in_final_within_5km, na.rm = TRUE)
  )
)

readr::write_csv(
  summary_tbl,
  file.path(audit_dir, "01_occurrence_audit_summary.csv")
)

# ------------------------------------------------------------
# 10. Console report
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("REV2 OCCURRENCE AUDIT COMPLETE\n")
cat("============================================================\n\n")
print(summary_tbl, n = Inf)
cat("\nSara geometry diagnostic:\n")
print(crs_diag)
cat("\nTarget-country counts written to:\n")
cat(file.path(audit_dir, "01E_target_country_counts.csv"), "\n")
cat("\nRecord-level Sara reconciliation written to:\n")
cat(file.path(audit_dir, "01G_sara_vs_REV1_record_level_audit.csv"), "\n")
cat("\nNo REV1 files were modified. No models were fitted.\n")
cat("============================================================\n")
