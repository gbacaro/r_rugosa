# ============================================================
# ROSA RUGOSA - REVISION 2
# 03_reconcile_Ribotta_with_REV1.R
#
# IMPORTANT
# ---------
# This script supersedes the previous REV2 scripts 03, 04 and 05
# for all analyses involving the European coastal domain.
#
# It uses the SAME raster-valid coastal mask logic as the official
# REV1 pipeline:
#   current_europe_coast1km_selectedvars.tif
#   -> complete-data mask across the 10 current predictors
#   -> direct cell assignment
#   -> nearest valid cell only when distance <= 1500 m
#
# No occurrence is added to the modelling dataset.
# No model is fitted.
# No REV1 file is modified.
# ============================================================

rm(list = ls())
gc()

required_pkgs <- c(
  "terra", "sf", "dplyr", "readr", "tibble"
)

missing_pkgs <- required_pkgs[
  !vapply(
    required_pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing packages: ",
    paste(missing_pkgs, collapse = ", ")
  )
}

library(terra)
library(sf)
library(dplyr)
library(readr)
library(tibble)

sf::sf_use_s2(TRUE)

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

rev2_dir <- normalizePath(
  getwd(),
  winslash = "/",
  mustWork = TRUE
)

paths_file <- file.path(
  rev2_dir,
  "02_scripts",
  "REV2_paths.R"
)

if (file.exists(paths_file)) {
  source(paths_file)
}

if (!exists("rev1_dir")) {
  rev1_dir <- normalizePath(
    "C:/Users/works/Desktop/PAPER ROSA RUGOSA/Rosa rugosa bis/REVISIONE_1_ROSA",
    winslash = "/",
    mustWork = TRUE
  )
}

if (!exists("input_dir")) input_dir <- file.path(rev2_dir, "01_input")
if (!exists("audit_dir")) audit_dir <- file.path(rev2_dir, "03_audit")
if (!exists("work_dir"))  work_dir  <- file.path(rev2_dir, "04_work")

dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1B. CLEAN OBSOLETE REV2 AUDIT ARTIFACTS
# ------------------------------------------------------------

obsolete_scripts <- file.path(
  rev2_dir,
  "02_scripts",
  c(
    "03_reconcile_Ribotta_with_REV1_CORRECTED.R",
    "04_trace_missing_Ribotta_cells.R",
    "05_audit_Ribotta_geography_by_country.R"
  )
)

obsolete_audit_patterns <- c(
  "^03[A-I]_",
  "^03COR_",
  "^03_MISSING_",
  "^04[A-F]_",
  "^05[A-I]_",
  "^REV2_AUDIT_SUMMARY_STAGE1\\.csv$"
)

obsolete_audit_files <- list.files(
  audit_dir,
  full.names = TRUE
)

if (length(obsolete_audit_files) > 0) {
  obsolete_audit_files <- obsolete_audit_files[
    Reduce(
      `|`,
      lapply(
        obsolete_audit_patterns,
        function(p) grepl(p, basename(obsolete_audit_files))
      )
    )
  ]
}

obsolete_work_files <- file.path(
  work_dir,
  c(
    "03_Ribotta_REV1_reconciliation.rds",
    "03COR_Ribotta_REV1_reconciliation.rds",
    "04_missing_Ribotta_cells_trace.rds",
    "05_Ribotta_geography_audit.rds"
  )
)

obsolete_files <- unique(
  c(
    obsolete_scripts,
    obsolete_audit_files,
    obsolete_work_files
  )
)

obsolete_files <- obsolete_files[
  file.exists(obsolete_files)
]

if (length(obsolete_files) > 0) {
  unlink(
    obsolete_files,
    recursive = TRUE,
    force = TRUE
  )
}

cat(
  "\nRemoved obsolete REV2 audit artifacts:",
  length(obsolete_files),
  "\n"
)

rdata_revision <- file.path(
  rev1_dir,
  "dati_per_revisione.RData"
)

if (!file.exists(rdata_revision)) {
  stop(
    "Cannot find REV1 RData: ",
    rdata_revision
  )
}

# ------------------------------------------------------------
# 2. Locate the OFFICIAL REV1 European coastal mask
#
# Priority:
# 1) saved mask from the official REV1 pipeline:
#    europe_valid_coastal_cells_1km.rds
# 2) original current CHELSA raster:
#    current_europe_coast1km_selectedvars.tif
#
# Search REV1 plus its parent project directories because large
# raster assets may have been stored outside REVISIONE_1_ROSA.
# ------------------------------------------------------------

search_roots <- unique(
  normalizePath(
    c(
      rev1_dir,
      dirname(rev1_dir),
      dirname(dirname(rev1_dir))
    ),
    winslash = "/",
    mustWork = FALSE
  )
)

search_roots <- search_roots[
  dir.exists(search_roots)
]

cat("\nSearching official REV1 spatial assets under:\n")
cat(
  paste0(" - ", search_roots, collapse = "\n"),
  "\n"
)

find_exact_file <- function(roots, filename) {

  out <- unlist(
    lapply(
      roots,
      function(root) {
        list.files(
          root,
          pattern = paste0(
            "^",
            gsub("\\.", "\\\\.", filename),
            "$"
          ),
          recursive = TRUE,
          full.names = TRUE,
          ignore.case = TRUE
        )
      }
    ),
    use.names = FALSE
  )

  unique(out[file.exists(out)])
}

mask_candidates <- find_exact_file(
  search_roots,
  "europe_valid_coastal_cells_1km.rds"
)

climate_candidates <- find_exact_file(
  search_roots,
  "current_europe_coast1km_selectedvars.tif"
)

candidate_vars <- c(
  "bio3", "bio8", "bio9", "bio10",
  "bio15", "bio18", "bio19",
  "fcf", "gsl", "gst"
)

mask_source_type <- NA_character_
mask_source_file <- NA_character_
climate_file <- NA_character_

if (length(mask_candidates) > 0) {

  # Prefer the canonical data_clean location when multiple copies exist.
  preferred_mask <- grep(
    "data_clean",
    mask_candidates,
    ignore.case = TRUE
  )

  if (length(preferred_mask) > 0) {
    mask_source_file <- mask_candidates[preferred_mask[1]]
  } else {
    mask_source_file <- mask_candidates[1]
  }

  cat(
    "\nUsing saved REV1 valid coastal mask:\n",
    mask_source_file,
    "\n",
    sep = ""
  )

  valid_eu_mask <- readRDS(
    mask_source_file
  )

  if (!inherits(valid_eu_mask, "SpatRaster")) {
    stop(
      "The saved object europe_valid_coastal_cells_1km.rds ",
      "is not a terra SpatRaster."
    )
  }

  if (terra::nlyr(valid_eu_mask) != 1) {
    stop(
      "The saved European valid coastal mask should contain one layer, ",
      "but contains ",
      terra::nlyr(valid_eu_mask),
      "."
    )
  }

  mask_source_type <- "saved_REV1_mask_RDS"

} else if (length(climate_candidates) > 0) {

  preferred_tif <- grep(
    "data_raw[/\\\\]climate_chelsa[/\\\\]current",
    climate_candidates,
    ignore.case = TRUE
  )

  if (length(preferred_tif) > 0) {
    climate_file <- climate_candidates[preferred_tif[1]]
  } else {
    climate_file <- climate_candidates[1]
  }

  cat(
    "\nSaved mask not found. Rebuilding from official REV1 raster:\n",
    climate_file,
    "\n",
    sep = ""
  )

  bio_eu <- terra::rast(
    climate_file
  )

  available_vars <- intersect(
    candidate_vars,
    names(bio_eu)
  )

  if (length(available_vars) == 0) {
    stop(
      "None of the expected current predictor layers were found. ",
      "Raster layers are: ",
      paste(names(bio_eu), collapse = ", ")
    )
  }

  build_complete_mask <- function(r, vars) {

    vars <- intersect(
      vars,
      names(r)
    )

    mask <- !is.na(
      r[[vars[1]]]
    )

    if (length(vars) > 1) {

      for (v in vars[-1]) {
        mask <- mask &
          !is.na(
            r[[v]]
          )
      }
    }

    terra::ifel(
      mask,
      1,
      NA
    )
  }

  valid_eu_mask <- build_complete_mask(
    bio_eu,
    candidate_vars
  )

  mask_source_type <- "rebuilt_from_REV1_current_TIF"
  mask_source_file <- climate_file

} else {

  # Write a diagnostic before stopping, so the missing dependency is explicit.
  asset_diag <- data.frame(
    searched_root = search_roots,
    stringsAsFactors = FALSE
  )

  readr::write_csv(
    asset_diag,
    file.path(
      audit_dir,
      "03_MISSING_SPATIAL_ASSET_search_roots.csv"
    )
  )

  stop(
    "Neither europe_valid_coastal_cells_1km.rds nor ",
    "current_europe_coast1km_selectedvars.tif could be found in the ",
    "REV1 directory or its parent project directories. ",
    "The search roots have been written to ",
    file.path(
      audit_dir,
      "03_MISSING_SPATIAL_ASSET_search_roots.csv"
    ),
    ". Do not reconstruct the domain from eu_df."
  )
}

# ------------------------------------------------------------
# 3. Validate the official mask
# ------------------------------------------------------------

n_valid_cells <- as.integer(
  terra::global(
    valid_eu_mask,
    "sum",
    na.rm = TRUE
  )[1, 1]
)

expected_domain_cells <- 79579L

domain_count_matches_REV1 <-
  isTRUE(
    n_valid_cells == expected_domain_cells
  )

cat(
  "\nSpatial source type:",
  mask_source_type,
  "\n"
)

cat(
  "Spatial source file:",
  mask_source_file,
  "\n"
)

cat(
  "Reconstructed raster-valid European coastal cells:",
  n_valid_cells,
  "\n"
)

cat(
  "Matches REV1 expected 79,579 cells:",
  domain_count_matches_REV1,
  "\n"
)

if (!domain_count_matches_REV1) {

  stop(
    "STOP: the recovered REV1 coastal mask contains ",
    n_valid_cells,
    " valid cells instead of the expected 79,579. ",
    "No reconciliation output will be written until the spatial ",
    "domain is reproduced exactly."
  )
}

# ------------------------------------------------------------
# 4. Load REV1 occurrence stages
# ------------------------------------------------------------

rev1 <- new.env(
  parent = baseenv()
)

load(
  rdata_revision,
  envir = rev1
)

get_rev1 <- function(
  name,
  required = TRUE
) {

  if (
    !exists(
      name,
      envir = rev1,
      inherits = FALSE
    )
  ) {

    if (required) {
      stop(
        "Required REV1 object not found: ",
        name
      )
    }

    return(NULL)
  }

  get(
    name,
    envir = rev1,
    inherits = FALSE
  )
}

flagged  <- get_rev1("flagged")
clean    <- get_rev1("clean")
coastal  <- get_rev1("coastal")
main_occ <- get_rev1("main_occ")
occ      <- get_rev1("occ")
world    <- get_rev1("world", required = FALSE)

if (!"cell_id" %in% names(coastal)) {
  stop("REV1 object 'coastal' does not contain cell_id.")
}

if (!"cell_id" %in% names(main_occ)) {
  stop("REV1 object 'main_occ' does not contain cell_id.")
}

if (!"cell_id" %in% names(occ)) {
  stop("REV1 object 'occ' does not contain cell_id.")
}

# ------------------------------------------------------------
# 5. Locate and read Sara/Ribotta shapefile
# ------------------------------------------------------------

shp_files <- list.files(
  rev2_dir,
  pattern = "\\.shp$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(shp_files) == 0) {

  zip_files <- list.files(
    rev2_dir,
    pattern = "\\.zip$",
    recursive = FALSE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (length(zip_files) == 0) {
    stop(
      "No Sara/Ribotta shapefile or ZIP file was found."
    )
  }

  sara_dir <- file.path(
    input_dir,
    "sara_ribotta"
  )

  dir.create(
    sara_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  unzip(
    zip_files[1],
    exdir = sara_dir,
    overwrite = TRUE
  )

  shp_files <- list.files(
    sara_dir,
    pattern = "\\.shp$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
}

if (length(shp_files) == 0) {
  stop(
    "No shapefile found for Sara/Ribotta data."
  )
}

# Select the shapefile with the greatest number of features.
shp_info <- lapply(
  shp_files,
  function(f) {

    x <- try(
      sf::st_read(
        f,
        quiet = TRUE
      ),
      silent = TRUE
    )

    if (inherits(x, "try-error")) {
      return(NULL)
    }

    tibble(
      file = f,
      n_features = nrow(x)
    )
  }
)

shp_info <- bind_rows(
  shp_info
)

if (nrow(shp_info) == 0) {
  stop(
    "No readable Sara/Ribotta shapefile found."
  )
}

sara_file <- shp_info$file[
  which.max(
    shp_info$n_features
  )
]

cat(
  "\nSara/Ribotta shapefile:\n",
  sara_file,
  "\n",
  sep = ""
)

sara_raw <- sf::st_read(
  sara_file,
  quiet = TRUE
)

if (is.na(sf::st_crs(sara_raw))) {
  stop(
    "Sara/Ribotta shapefile has no declared CRS."
  )
}

sara_raw <- sf::st_transform(
  sara_raw,
  4326
)

sara <- suppressWarnings(
  sf::st_cast(
    sara_raw,
    "POINT"
  )
)

sara$sara_id <- seq_len(
  nrow(sara)
)

sara_xy <- sf::st_coordinates(
  sara
)

sara$sara_lon <- sara_xy[, "X"]
sara$sara_lat <- sara_xy[, "Y"]

# ------------------------------------------------------------
# 6. Assign Sara points to the official REV1 valid coastal mask
#    using the SAME rule as the REV1 pipeline
# ------------------------------------------------------------

assign_to_valid_mask <- function(
  points_sf,
  mask_rast,
  max_snap_m = 1500
) {

  xy <- sf::st_coordinates(
    points_sf
  )

  cell_guess <- terra::cellFromXY(
    mask_rast,
    xy
  )

  mask_vals <- terra::values(
    mask_rast,
    mat = FALSE
  )

  exact_ok <-
    !is.na(cell_guess) &
    !is.na(mask_vals[cell_guess])

  out <- tibble(
    point_row = seq_len(
      nrow(points_sf)
    ),

    exact_cell_guess =
      cell_guess,

    exact_valid =
      exact_ok,

    assigned_cell_id =
      NA_integer_,

    assigned_lon =
      NA_real_,

    assigned_lat =
      NA_real_,

    snap_distance_m =
      NA_real_,

    snap_type =
      NA_character_
  )

  # Exact valid cells
  if (any(exact_ok)) {

    exact_cells <- cell_guess[
      exact_ok
    ]

    exact_xy <- terra::xyFromCell(
      mask_rast,
      exact_cells
    )

    out$assigned_cell_id[
      exact_ok
    ] <- exact_cells

    out$assigned_lon[
      exact_ok
    ] <- exact_xy[, 1]

    out$assigned_lat[
      exact_ok
    ] <- exact_xy[, 2]

    out$snap_distance_m[
      exact_ok
    ] <- 0

    out$snap_type[
      exact_ok
    ] <- "exact_valid_cell"
  }

  # Nearest valid cell for points outside the valid mask
  remaining <- which(
    !exact_ok
  )

  if (length(remaining) > 0) {

    valid_df <- terra::as.data.frame(
      mask_rast,
      xy = TRUE,
      cells = TRUE,
      na.rm = TRUE
    )

    valid_df <- valid_df[
      ,
      c(
        "cell",
        "x",
        "y"
      )
    ]

    valid_sf <- sf::st_as_sf(
      valid_df,
      coords = c("x", "y"),
      crs = 4326
    )

    rem_sf <- points_sf[
      remaining,
    ]

    # Use the same projected CRS used in the REV1 snapping function.
    valid_m <- sf::st_transform(
      valid_sf,
      3857
    )

    rem_m <- sf::st_transform(
      rem_sf,
      3857
    )

    idx <- sf::st_nearest_feature(
      rem_m,
      valid_m
    )

    d <- as.numeric(
      sf::st_distance(
        rem_m,
        valid_m[idx, ],
        by_element = TRUE
      )
    )

    keep <- is.finite(d) &
      d <= max_snap_m

    if (any(keep)) {

      retained_rows <- remaining[
        keep
      ]

      valid_idx <- idx[
        keep
      ]

      out$assigned_cell_id[
        retained_rows
      ] <- valid_df$cell[
        valid_idx
      ]

      out$assigned_lon[
        retained_rows
      ] <- valid_df$x[
        valid_idx
      ]

      out$assigned_lat[
        retained_rows
      ] <- valid_df$y[
        valid_idx
      ]

      out$snap_distance_m[
        retained_rows
      ] <- d[
        keep
      ]

      out$snap_type[
        retained_rows
      ] <- "nearest_valid_cell"
    }

    # Keep actual nearest distance even when > 1500 m for audit.
    not_kept <- !keep

    if (any(not_kept)) {

      rejected_rows <- remaining[
        not_kept
      ]

      out$snap_distance_m[
        rejected_rows
      ] <- d[
        not_kept
      ]

      out$snap_type[
        rejected_rows
      ] <- "outside_1500m"
    }
  }

  out$retained_by_REV1_spatial_rule <-
    !is.na(
      out$assigned_cell_id
    )

  out
}

sara_assignment <- assign_to_valid_mask(
  sara,
  valid_eu_mask,
  max_snap_m = 1500
)

sara_df <- bind_cols(
  sf::st_drop_geometry(
    sara
  ),
  sara_assignment %>%
    select(
      -point_row
    )
)

# ------------------------------------------------------------
# 7. Robust country assignment
#
# Country names from Natural Earth are treated as the primary source.
# ISO codes are derived from the country name when the native ISO field
# is missing or coded as -99. Points falling just outside a polygon
# because of coastline generalisation are assigned to the nearest
# country only when the distance is <= 20 km.
# ------------------------------------------------------------

country_name_to_iso2 <- function(x) {

  y <- toupper(trimws(as.character(x)))

  map <- c(
    "UNITED KINGDOM" = "GB",
    "FRANCE" = "FR",
    "ITALY" = "IT",
    "SPAIN" = "ES",
    "PORTUGAL" = "PT",
    "GREECE" = "GR",
    "CROATIA" = "HR",
    "SLOVENIA" = "SI",
    "MONTENEGRO" = "ME",
    "ALBANIA" = "AL",
    "BOSNIA AND HERZEGOVINA" = "BA",
    "TURKEY" = "TR",
    "TÜRKIYE" = "TR",
    "BULGARIA" = "BG",
    "ROMANIA" = "RO",
    "CYPRUS" = "CY",
    "MALTA" = "MT",
    "GERMANY" = "DE",
    "DENMARK" = "DK",
    "NORWAY" = "NO",
    "SWEDEN" = "SE",
    "FINLAND" = "FI",
    "POLAND" = "PL",
    "NETHERLANDS" = "NL",
    "BELGIUM" = "BE",
    "IRELAND" = "IE",
    "ESTONIA" = "EE",
    "LATVIA" = "LV",
    "LITHUANIA" = "LT",
    "ICELAND" = "IS",
    "UKRAINE" = "UA",
    "GEORGIA" = "GE",
    "RUSSIA" = "RU",
    "RUSSIAN FEDERATION" = "RU"
  )

  unname(
    map[
      match(
        y,
        names(map)
      )
    ]
  )
}

sara_df$country_code <- NA_character_
sara_df$country_name <- NA_character_
sara_df$country_assignment_method <- NA_character_
sara_df$country_polygon_distance_m <- NA_real_

if (
  !is.null(world) &&
  inherits(world, "sf")
) {

  world_4326 <- sf::st_transform(
    world,
    4326
  )

  world_names <- names(
    world_4326
  )

  name_candidates <- c(
    "admin",
    "ADMIN",
    "name",
    "NAME",
    "sovereignt",
    "SOVEREIGNT"
  )

  name_col <- name_candidates[
    name_candidates %in%
      world_names
  ][1]

  if (is.na(name_col)) {
    stop(
      "Could not identify a country-name field in the REV1 world polygons."
    )
  }

  # One-to-one polygon assignment without row multiplication.
  within_list <- sf::st_within(
    sara,
    world_4326
  )

  inside_idx <- vapply(
    within_list,
    function(z) {
      if (length(z) == 0) NA_integer_ else z[1]
    },
    integer(1)
  )

  inside <- !is.na(
    inside_idx
  )

  if (any(inside)) {

    sara_df$country_name[
      inside
    ] <- as.character(
      world_4326[[name_col]][
        inside_idx[inside]
      ]
    )

    sara_df$country_assignment_method[
      inside
    ] <- "polygon_within"

    sara_df$country_polygon_distance_m[
      inside
    ] <- 0
  }

  # Coastal points may fall just outside medium-resolution polygons.
  outside <- which(
    !inside
  )

  if (length(outside) > 0) {

    sara_out_m <- sf::st_transform(
      sara[
        outside,
      ],
      3035
    )

    world_m <- sf::st_transform(
      world_4326,
      3035
    )

    nearest_country_idx <- sf::st_nearest_feature(
      sara_out_m,
      world_m
    )

    nearest_country_dist <- as.numeric(
      sf::st_distance(
        sara_out_m,
        world_m[
          nearest_country_idx,
        ],
        by_element = TRUE
      )
    )

    near_enough <- is.finite(
      nearest_country_dist
    ) &
      nearest_country_dist <= 20000

    if (any(near_enough)) {

      rows <- outside[
        near_enough
      ]

      poly_idx <- nearest_country_idx[
        near_enough
      ]

      sara_df$country_name[
        rows
      ] <- as.character(
        world_4326[[name_col]][
          poly_idx
        ]
      )

      sara_df$country_assignment_method[
        rows
      ] <- "nearest_polygon_le20km"

      sara_df$country_polygon_distance_m[
        rows
      ] <- nearest_country_dist[
        near_enough
      ]
    }
  }

  sara_df$country_code <-
    country_name_to_iso2(
      sara_df$country_name
    )
}

# ------------------------------------------------------------
# 8. Compare assigned cells directly with REV1 occurrence cells
# ------------------------------------------------------------

old_coastal_cells <- unique(
  coastal$cell_id
)

corrected_allrange_cells <- unique(
  main_occ$cell_id
)

final_invaded_cells <- unique(
  occ$cell_id
)

sara_df$cell_in_old_coastal <-
  sara_df$assigned_cell_id %in%
  old_coastal_cells

sara_df$cell_in_corrected_allrange <-
  sara_df$assigned_cell_id %in%
  corrected_allrange_cells

sara_df$cell_in_final_invaded <-
  sara_df$assigned_cell_id %in%
  final_invaded_cells

# ------------------------------------------------------------
# 9. Nearest raw GBIF distance
#    Used only for provenance / independent-field-data audit.
# ------------------------------------------------------------

lon_candidates <- c(
  "decimal_longitude",
  "decimalLongitude",
  "longitude",
  "lon"
)

lat_candidates <- c(
  "decimal_latitude",
  "decimalLatitude",
  "latitude",
  "lat"
)

lon_col <- lon_candidates[
  lon_candidates %in%
    names(flagged)
][1]

lat_col <- lat_candidates[
  lat_candidates %in%
    names(flagged)
][1]

if (
  is.na(lon_col) ||
  is.na(lat_col)
) {

  stop(
    "Could not identify longitude/latitude fields in REV1 flagged data."
  )
}

flag_lon <- suppressWarnings(
  as.numeric(
    flagged[[lon_col]]
  )
)

flag_lat <- suppressWarnings(
  as.numeric(
    flagged[[lat_col]]
  )
)

flag_ok <-
  is.finite(flag_lon) &
  is.finite(flag_lat) &
  abs(flag_lon) <= 180 &
  abs(flag_lat) <= 90

flagged_spatial <- flagged[
  flag_ok,
  ,
  drop = FALSE
]

flagged_spatial$.flagged_row_original <-
  which(flag_ok)

flagged_spatial$.lon <-
  flag_lon[
    flag_ok
  ]

flagged_spatial$.lat <-
  flag_lat[
    flag_ok
  ]

flagged_sf <- sf::st_as_sf(
  flagged_spatial,
  coords = c(
    ".lon",
    ".lat"
  ),
  crs = 4326,
  remove = FALSE
)

sara_m <- sf::st_transform(
  sara,
  3035
)

flagged_m <- sf::st_transform(
  flagged_sf,
  3035
)

nearest_idx <- sf::st_nearest_feature(
  sara_m,
  flagged_m
)

nearest_dist <- as.numeric(
  sf::st_distance(
    sara_m,
    flagged_m[
      nearest_idx,
    ],
    by_element = TRUE
  )
)

sara_df$nearest_flagged_distance_m <-
  nearest_dist

sara_df$nearest_flagged_original_row <-
  flagged_spatial$.flagged_row_original[
    nearest_idx
  ]


# Country reported by the nearest GBIF record, used as an independent
# provenance check rather than as the primary geographic assignment.
flag_country_candidates <- c(
  "country_code",
  "countryCode",
  "country"
)

flag_country_col <- flag_country_candidates[
  flag_country_candidates %in%
    names(flagged_spatial)
][1]

sara_df$country_code_nearest_GBIF <- NA_character_

if (!is.na(flag_country_col)) {

  raw_gbif_country <- as.character(
    flagged_spatial[[flag_country_col]][
      nearest_idx
    ]
  )

  raw_upper <- toupper(
    trimws(
      raw_gbif_country
    )
  )

  already_iso2 <- grepl(
    "^[A-Z]{2}$",
    raw_upper
  )

  sara_df$country_code_nearest_GBIF[
    already_iso2
  ] <- raw_upper[
    already_iso2
  ]

  needs_name_map <- !already_iso2 &
    !is.na(raw_upper) &
    raw_upper != ""

  if (any(needs_name_map)) {

    sara_df$country_code_nearest_GBIF[
      needs_name_map
    ] <- country_name_to_iso2(
      raw_upper[
        needs_name_map
      ]
    )
  }
}

sara_df$country_polygon_vs_GBIF_agree <-
  ifelse(
    is.na(sara_df$country_code) |
      is.na(sara_df$country_code_nearest_GBIF),
    NA,
    sara_df$country_code ==
      sara_df$country_code_nearest_GBIF
  )

# ------------------------------------------------------------
# 10. Point-level status
# ------------------------------------------------------------

sara_df$gbif_match_100m <-
  sara_df$nearest_flagged_distance_m <= 100

sara_df$gbif_match_1km <-
  sara_df$nearest_flagged_distance_m <= 1000

sara_df$gbif_match_5km <-
  sara_df$nearest_flagged_distance_m <= 5000

sara_df$potential_independent_field_record <-
  sara_df$retained_by_REV1_spatial_rule &
  !sara_df$gbif_match_1km

# ------------------------------------------------------------
# 11. Cell-level reconciliation
# ------------------------------------------------------------

sara_valid <- sara_df %>%
  filter(
    retained_by_REV1_spatial_rule,
    !is.na(
      assigned_cell_id
    )
  )

cell_audit <- sara_valid %>%
  group_by(
    assigned_cell_id
  ) %>%
  summarise(
    n_sara_points = n(),

    country_code =
      first(
        country_code[
          !is.na(country_code)
        ],
        default = NA_character_
      ),

    country_name =
      first(
        country_name[
          !is.na(country_name)
        ],
        default = NA_character_
      ),

    min_GBIF_distance_m =
      min(
        nearest_flagged_distance_m,
        na.rm = TRUE
      ),

    cell_in_old_coastal =
      any(
        cell_in_old_coastal
      ),

    cell_in_corrected_allrange =
      any(
        cell_in_corrected_allrange
      ),

    cell_in_final_invaded =
      any(
        cell_in_final_invaded
      ),

    .groups = "drop"
  )

cell_audit$min_GBIF_distance_m[
  !is.finite(
    cell_audit$min_GBIF_distance_m
  )
] <- NA_real_

cell_audit$sara_cell_absent_final <-
  !cell_audit$cell_in_final_invaded

cell_audit$potential_independent_cell <-
  cell_audit$sara_cell_absent_final &
  (
    is.na(
      cell_audit$min_GBIF_distance_m
    ) |
    cell_audit$min_GBIF_distance_m > 1000
  )

# ------------------------------------------------------------
# 12. Country summaries
# ------------------------------------------------------------

point_country_summary <- sara_df %>%
  mutate(
    country_code =
      ifelse(
        is.na(country_code) |
          country_code == "",
        "UNRESOLVED",
        country_code
      ),

    country_name =
      ifelse(
        is.na(country_name) |
          country_name == "",
        country_code,
        country_name
      )
  ) %>%
  group_by(
    country_code,
    country_name
  ) %>%
  summarise(
    n_sara_points = n(),

    exact_valid_points =
      sum(
        exact_valid,
        na.rm = TRUE
      ),

    snapped_within_1500m =
      sum(
        snap_type ==
          "nearest_valid_cell",
        na.rm = TRUE
      ),

    retained_by_REV1_spatial_rule =
      sum(
        retained_by_REV1_spatial_rule,
        na.rm = TRUE
      ),

    outside_1500m =
      sum(
        !retained_by_REV1_spatial_rule,
        na.rm = TRUE
      ),

    GBIF_within_100m =
      sum(
        gbif_match_100m,
        na.rm = TRUE
      ),

    GBIF_within_1km =
      sum(
        gbif_match_1km,
        na.rm = TRUE
      ),

    country_assigned_by_polygon =
      sum(
        country_assignment_method == "polygon_within",
        na.rm = TRUE
      ),

    country_assigned_by_nearest_polygon =
      sum(
        country_assignment_method == "nearest_polygon_le20km",
        na.rm = TRUE
      ),

    polygon_GBIF_country_agreement =
      sum(
        country_polygon_vs_GBIF_agree %in% TRUE,
        na.rm = TRUE
      ),

    polygon_GBIF_country_disagreement =
      sum(
        country_polygon_vs_GBIF_agree %in% FALSE,
        na.rm = TRUE
      ),

    final_REV1_cell =
      sum(
        retained_by_REV1_spatial_rule &
          cell_in_final_invaded,
        na.rm = TRUE
      ),

    .groups = "drop"
  ) %>%
  arrange(
    desc(
      n_sara_points
    )
  )

cell_country_summary <- cell_audit %>%
  mutate(
    country_code =
      ifelse(
        is.na(country_code) |
          country_code == "",
        "UNRESOLVED",
        country_code
      ),

    country_name =
      ifelse(
        is.na(country_name) |
          country_name == "",
        country_code,
        country_name
      )
  ) %>%
  group_by(
    country_code,
    country_name
  ) %>%
  summarise(
    n_unique_Sara_valid_cells =
      n(),

    n_cells_in_final_REV1 =
      sum(
        cell_in_final_invaded,
        na.rm = TRUE
      ),

    n_cells_absent_final_REV1 =
      sum(
        sara_cell_absent_final,
        na.rm = TRUE
      ),

    n_potential_independent_cells =
      sum(
        potential_independent_cell,
        na.rm = TRUE
      ),

    .groups = "drop"
  ) %>%
  arrange(
    desc(
      n_unique_Sara_valid_cells
    )
  )

reviewer_codes <- c(
  "GB", "FR", "IT",
  "ES", "PT",
  "SI", "HR", "BA", "ME", "AL", "GR",
  "TR", "BG", "RO",
  "CY", "MT"
)

reviewer_points_summary <-
  point_country_summary %>%
  filter(
    country_code %in%
      reviewer_codes
  )

reviewer_cells_summary <-
  cell_country_summary %>%
  filter(
    country_code %in%
      reviewer_codes
  )

# ------------------------------------------------------------
# 13. Critical reviewer-country points
# ------------------------------------------------------------

critical_points <- sara_df %>%
  filter(
    country_code %in%
      reviewer_codes
  ) %>%
  arrange(
    country_code,
    retained_by_REV1_spatial_rule,
    snap_distance_m
  )

independent_candidates <- sara_df %>%
  filter(
    potential_independent_field_record
  ) %>%
  arrange(
    country_code,
    desc(
      nearest_flagged_distance_m
    )
  )

independent_reviewer_candidates <-
  independent_candidates %>%
  filter(
    country_code %in%
      reviewer_codes
  )

# ------------------------------------------------------------
# 14. Global validation summary
# ------------------------------------------------------------

summary_tbl <- tibble(
  metric = c(
    "Official reconstructed European valid coastal cells",
    "Expected REV1 European valid coastal cells",
    "Domain cell count matches REV1",
    "Sara/Ribotta point geometries",
    "Sara points directly in valid coastal cell",
    "Sara points snapped to nearest valid cell <=1500 m",
    "Sara points retained by REV1 spatial rule",
    "Sara points farther than 1500 m from valid coastal mask",
    "Sara unique valid coastal cells",
    "Sara valid cells represented in final REV1",
    "Sara valid cells absent from final REV1",
    "Sara valid cells potentially independent from GBIF",
    "Sara points with GBIF within 100 m",
    "Sara points with GBIF within 1 km",
    "Sara points with GBIF within 5 km"
  ),

  value = c(
    as.character(n_valid_cells),
    as.character(expected_domain_cells),
    as.character(domain_count_matches_REV1),
    as.character(nrow(sara_df)),
    as.character(
      sum(
        sara_df$exact_valid,
        na.rm = TRUE
      )
    ),
    as.character(
      sum(
        sara_df$snap_type ==
          "nearest_valid_cell",
        na.rm = TRUE
      )
    ),
    as.character(
      sum(
        sara_df$retained_by_REV1_spatial_rule,
        na.rm = TRUE
      )
    ),
    as.character(
      sum(
        !sara_df$retained_by_REV1_spatial_rule,
        na.rm = TRUE
      )
    ),
    as.character(
      nrow(
        cell_audit
      )
    ),
    as.character(
      sum(
        cell_audit$cell_in_final_invaded,
        na.rm = TRUE
      )
    ),
    as.character(
      sum(
        cell_audit$sara_cell_absent_final,
        na.rm = TRUE
      )
    ),
    as.character(
      sum(
        cell_audit$potential_independent_cell,
        na.rm = TRUE
      )
    ),
    as.character(
      sum(
        sara_df$gbif_match_100m,
        na.rm = TRUE
      )
    ),
    as.character(
      sum(
        sara_df$gbif_match_1km,
        na.rm = TRUE
      )
    ),
    as.character(
      sum(
        sara_df$gbif_match_5km,
        na.rm = TRUE
      )
    )
  )
)

# ------------------------------------------------------------
# 15. Save corrected outputs
# ------------------------------------------------------------

write_csv(
  summary_tbl,
  file.path(
    audit_dir,
    "03COR_A_reconciliation_summary.csv"
  )
)

write_csv(
  sara_df,
  file.path(
    audit_dir,
    "03COR_B_Sara_point_reconciliation.csv"
  )
)

write_csv(
  cell_audit,
  file.path(
    audit_dir,
    "03COR_C_Sara_cell_reconciliation.csv"
  )
)

write_csv(
  point_country_summary,
  file.path(
    audit_dir,
    "03COR_D_point_summary_by_country.csv"
  )
)

write_csv(
  cell_country_summary,
  file.path(
    audit_dir,
    "03COR_E_cell_summary_by_country.csv"
  )
)

write_csv(
  reviewer_points_summary,
  file.path(
    audit_dir,
    "03COR_F_reviewer_point_summary.csv"
  )
)


country_assignment_diagnostic <- sara_df %>%
  count(
    country_code,
    country_name,
    country_assignment_method,
    country_code_nearest_GBIF,
    country_polygon_vs_GBIF_agree,
    name = "n_points"
  ) %>%
  arrange(
    country_code,
    desc(n_points)
  )

write_csv(
  country_assignment_diagnostic,
  file.path(
    audit_dir,
    "03COR_F2_country_assignment_diagnostic.csv"
  )
)

write_csv(
  reviewer_cells_summary,
  file.path(
    audit_dir,
    "03COR_G_reviewer_cell_summary.csv"
  )
)

write_csv(
  critical_points,
  file.path(
    audit_dir,
    "03COR_H_all_reviewer_country_points.csv"
  )
)

write_csv(
  independent_reviewer_candidates,
  file.path(
    audit_dir,
    "03COR_I_candidate_independent_reviewer_points.csv"
  )
)

saveRDS(
  list(
    summary = summary_tbl,
    sara_points = sara_df,
    cell_audit = cell_audit,
    point_country_summary = point_country_summary,
    cell_country_summary = cell_country_summary,
    reviewer_points_summary = reviewer_points_summary,
    reviewer_cells_summary = reviewer_cells_summary,
    independent_candidates = independent_candidates,
    mask_source_type = mask_source_type,
    mask_source_file = mask_source_file,
    climate_file = climate_file,
    domain_cells = n_valid_cells
  ),
  file.path(
    work_dir,
    "03COR_Ribotta_REV1_reconciliation.rds"
  )
)

# ------------------------------------------------------------
# 16. WRITE CANONICAL REV2 WORKFLOW MANIFEST
# ------------------------------------------------------------

canonical_workflow <- c(
  "ROSA RUGOSA REV2 - CANONICAL WORKFLOW",
  "",
  "VALID SCRIPTS",
  "00_setup_REV2.R",
  "01_occurrence_audit_REV2.R",
  "02_model_input_audit_REV2.R",
  "03_reconcile_Ribotta_with_REV1.R",
  "",
  "SUPERSEDED / INVALIDATED",
  "03_reconcile_Ribotta_with_REV1_CORRECTED.R",
  "04_trace_missing_Ribotta_cells.R",
  "05_audit_Ribotta_geography_by_country.R",
  "",
  "IMPORTANT",
  "Previous 03/04/05 audit outputs based on eu_df as a geographic grid were deleted.",
  "The canonical step 03 uses the official REV1 raster-valid coastal mask recovered from the saved RDS when available, otherwise rebuilt from the current CHELSA TIFF.",
  "Country assignment uses Natural Earth country names with a <=20 km nearest-polygon fallback for coastal points and is cross-checked against the nearest GBIF record.",
  "Only outputs prefixed 03COR_ and the RDS 03COR_Ribotta_REV1_reconciliation.rds are valid for the Ribotta reconciliation stage."
)

writeLines(
  canonical_workflow,
  file.path(
    audit_dir,
    "REV2_CANONICAL_WORKFLOW.txt"
  )
)

# ------------------------------------------------------------
# 17. Console report
# ------------------------------------------------------------

cat(
  "\n============================================================\n"
)

cat(
  "REV2 - CORRECTED RIBOTTA / REV1 RECONCILIATION\n"
)

cat(
  "============================================================\n\n"
)

print(
  as.data.frame(
    summary_tbl
  ),
  row.names = FALSE
)

cat(
  "\nReviewer-focused POINT summary:\n"
)

print(
  as.data.frame(
    reviewer_points_summary
  ),
  row.names = FALSE
)

cat(
  "\nReviewer-focused CELL summary:\n"
)

print(
  as.data.frame(
    reviewer_cells_summary
  ),
  row.names = FALSE
)

cat(
  "\nPotentially independent reviewer-country points:\n"
)

if (
  nrow(
    independent_reviewer_candidates
  ) == 0
) {

  cat(
    "None.\n"
  )

} else {

  print(
    as.data.frame(
      independent_reviewer_candidates
    ),
    row.names = FALSE
  )
}

cat(
  "\nIMPORTANT: outputs from previous REV2 scripts 03, 04 and 05 ",
  "that treated eu_df as the geographic coastal grid must be ignored.\n"
)

cat(
  "\nNo occurrence has been added. No model has been fitted.\n"
)

cat(
  "============================================================\n"
)
