# ============================================================
# ROSA RUGOSA - REVISION 2
# 05_rebuild_REV2_occurrences_FINAL7.R
#
# Purpose
# -------
# Rebuild the European invaded occurrence set on a corrected
# 1-km coastal domain using:
#   - the occurrence-quality filters declared in REV1;
#   - the unmasked European current-climate raster stack;
#   - a 1,000-m terrestrial coastal belt built in EPSG:3035;
#   - raster-cell centres as the coastal-domain inclusion rule;
#   - ONLY the seven predictors retained in the final SDM;
#   - <=1,500 m occurrence-to-cell snapping in EPSG:3035;
#   - one occurrence per raster cell.
#
# The previous EPSG:3857 coastal mask is not used.
# Sara/Ribotta data are audited but not added to calibration here.
# No model is fitted and no REV1 file is modified.
# ============================================================

rm(list = ls())
gc()

required_pkgs <- c(
  "terra", "sf", "dplyr", "readr", "tibble", "tidyr",
  "rnaturalearth", "rnaturalearthdata"
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
library(tidyr)
library(rnaturalearth)
library(rnaturalearthdata)

sf::sf_use_s2(TRUE)

# ------------------------------------------------------------
# 1. Paths and cleanup
# ------------------------------------------------------------

paths_file <- file.path(
  getwd(),
  "02_scripts",
  "REV2_paths.R"
)

if (!file.exists(paths_file)) {
  stop(
    "Missing 02_scripts/REV2_paths.R. Run 00_setup_REV2.R first."
  )
}

source(paths_file)

dir.create(
  audit_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  work_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

mask_file <- file.path(
  work_dir,
  "05_REV2_FINAL7_TRUE1KM_mask.tif"
)

belt_file <- file.path(
  work_dir,
  "05_REV2_europe_coastal_belt_1km_EPSG3035.rds"
)

invalidated_outputs <- c(
  list.files(
    audit_dir,
    pattern = "^(05|06|07)[A-Z]_",
    full.names = TRUE
  ),
  file.path(
    work_dir,
    "05_REV2_occurrences_FINAL7_VALID.rds"
  ),
  file.path(
    work_dir,
    "occurrences_REV2_FINAL7_GBIF.csv"
  ),
  mask_file,
  belt_file,
  file.path(
    work_dir,
    "06_REV2_occurrence_dataset_FROZEN.rds"
  ),
  file.path(
    work_dir,
    "occurrences_REV2_FROZEN.csv"
  ),
  file.path(
    work_dir,
    "background_REV2_FINAL7.csv"
  ),
  file.path(
    work_dir,
    "model_input_REV2_FINAL7.csv"
  ),
  file.path(
    work_dir,
    "07_REV2_background_FINAL7_FROZEN.rds"
  )
)

invalidated_outputs <- unique(
  invalidated_outputs[
    file.exists(
      invalidated_outputs
    )
  ]
)

if (length(invalidated_outputs) > 0) {
  unlink(
    invalidated_outputs,
    force = TRUE
  )
}

# ------------------------------------------------------------
# 2. Inputs
# ------------------------------------------------------------

precoast_file <- file.path(
  work_dir,
  "occurrences_rebuilt_precoastal_declared_filters.csv"
)

if (!file.exists(precoast_file)) {
  stop(
    "Missing canonical pre-coastal occurrence file: ",
    precoast_file,
    "\nRun 01_occurrence_audit_REV2.R first."
  )
}

step03_file <- file.path(
  work_dir,
  "03_Ribotta_REV1_reconciliation_VALID.rds"
)

if (!file.exists(step03_file)) {
  stop(
    "Missing canonical step-03 reconciliation RDS: ",
    step03_file
  )
}

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

find_exact_file <- function(
  roots,
  filename
) {

  out <- unlist(
    lapply(
      roots,
      function(root) {
        list.files(
          root,
          pattern = paste0(
            "^",
            gsub(
              "\\.",
              "\\\\.",
              filename
            ),
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

  unique(
    out[
      file.exists(out)
    ]
  )
}

bbox_candidates <- find_exact_file(
  search_roots,
  "current_europe_bbox_selectedvars.tif"
)

if (length(bbox_candidates) == 0) {
  stop(
    "Could not find current_europe_bbox_selectedvars.tif. ",
    "The corrected REV2 coastal domain must be rebuilt from the ",
    "unmasked European raster."
  )
}

preferred <- grep(
  "data_raw[/\\\\]climate_chelsa[/\\\\]current",
  bbox_candidates,
  ignore.case = TRUE
)

if (length(preferred) > 0) {
  bbox_file <- bbox_candidates[
    preferred[1]
  ]
} else {
  bbox_file <- bbox_candidates[1]
}

clim_bbox <- terra::rast(
  bbox_file
)

final7_vars <- c(
  "bio3",
  "bio8",
  "bio9",
  "bio15",
  "bio18",
  "gsl",
  "gst"
)

missing_final7 <- setdiff(
  final7_vars,
  names(clim_bbox)
)

if (length(missing_final7) > 0) {
  stop(
    "European bbox raster is missing FINAL7 predictors: ",
    paste(
      missing_final7,
      collapse = ", "
    )
  )
}

# ------------------------------------------------------------
# 3. Correct 1-km European coastal domain
# ------------------------------------------------------------

eu_ext_vec <- c(
  -15,
  45,
  30,
  72
)

distance_crs <- 3035
coastal_buffer_m <- 1000

bb <- sf::st_as_sfc(
  sf::st_bbox(
    c(
      xmin = eu_ext_vec[1],
      xmax = eu_ext_vec[2],
      ymin = eu_ext_vec[3],
      ymax = eu_ext_vec[4]
    ),
    crs = sf::st_crs(4326)
  )
)

coast <- rnaturalearth::ne_coastline(
  scale = "medium",
  returnclass = "sf"
)

land <- rnaturalearth::ne_countries(
  scale = "medium",
  returnclass = "sf"
)

coast <- suppressWarnings(
  sf::st_crop(
    coast,
    sf::st_bbox(bb)
  )
)

land <- suppressWarnings(
  sf::st_crop(
    land,
    sf::st_bbox(bb)
  )
)

coast_m <- sf::st_transform(
  coast,
  distance_crs
)

land_m <- sf::st_transform(
  sf::st_make_valid(
    land
  ),
  distance_crs
)

coast_union_m <- sf::st_union(
  coast_m
)

land_union_m <- sf::st_union(
  land_m
)

coast_buffer_m <- sf::st_buffer(
  coast_union_m,
  dist = coastal_buffer_m
)

belt_geom_m <- suppressWarnings(
  sf::st_intersection(
    land_union_m,
    coast_buffer_m
  )
)

belt_geom_m <- sf::st_make_valid(
  belt_geom_m
)

if (length(belt_geom_m) == 0) {
  stop(
    "Corrected 1-km coastal belt is empty."
  )
}

belt_m <- sf::st_sf(
  belt_id = seq_along(
    belt_geom_m
  ),
  geometry = belt_geom_m
)

belt_ll <- sf::st_transform(
  belt_m,
  4326
)

saveRDS(
  belt_ll,
  belt_file
)

# ------------------------------------------------------------
# 4. FINAL7 raster mask: cell centres inside corrected belt
# ------------------------------------------------------------

complete_mask <- function(
  r,
  vars
) {

  x <- !is.na(
    r[[vars[1]]]
  )

  if (length(vars) > 1) {
    for (v in vars[-1]) {
      x <- x &
        !is.na(
          r[[v]]
        )
    }
  }

  terra::ifel(
    x,
    1,
    NA
  )
}

final7_complete_bbox <- complete_mask(
  clim_bbox,
  final7_vars
)

belt_vect <- terra::vect(
  belt_ll
)

final7_complete_crop <- terra::crop(
  final7_complete_bbox,
  belt_vect,
  snap = "out"
)

belt_center_mask <- terra::rasterize(
  belt_vect,
  final7_complete_crop,
  field = 1,
  background = NA,
  touches = FALSE
)

final7_mask <- terra::ifel(
  final7_complete_crop == 1 &
    !is.na(
      belt_center_mask
    ),
  1,
  NA
)

terra::writeRaster(
  final7_mask,
  mask_file,
  overwrite = TRUE,
  wopt = list(
    gdal = c(
      "COMPRESS=LZW",
      "TILED=YES"
    )
  )
)

n_final7_cells <- as.integer(
  terra::global(
    final7_mask,
    "sum",
    na.rm = TRUE
  )[1, 1]
)

if (!is.finite(n_final7_cells) ||
    n_final7_cells < 1L) {
  stop(
    "Corrected FINAL7 mask contains no valid cells."
  )
}

valid_df <- terra::as.data.frame(
  final7_mask,
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

names(valid_df) <- c(
  "cell_id",
  "decimal_longitude",
  "decimal_latitude"
)

valid_df$cell_id <- as.integer(
  valid_df$cell_id
)

valid_sf <- sf::st_as_sf(
  valid_df,
  coords = c(
    "decimal_longitude",
    "decimal_latitude"
  ),
  crs = 4326,
  remove = FALSE
)

valid_m <- sf::st_transform(
  valid_sf,
  distance_crs
)

belt_area_km2 <- as.numeric(
  sum(
    sf::st_area(
      belt_m
    )
  )
) / 1e6

raster_area_km2 <- as.numeric(
  terra::global(
    terra::cellSize(
      final7_mask,
      unit = "km",
      mask = TRUE
    ),
    "sum",
    na.rm = TRUE
  )[1, 1]
)

domain_summary <- tibble(
  metric = c(
    "Coastal buffer distance m",
    "Distance CRS EPSG",
    "Raster inclusion rule",
    "FINAL7 valid coastal cells",
    "Vector coastal-belt area km2",
    "FINAL7 raster-domain area km2"
  ),
  value = c(
    as.character(
      coastal_buffer_m
    ),
    as.character(
      distance_crs
    ),
    "cell_center_inside_terrestrial_1km_belt",
    as.character(
      n_final7_cells
    ),
    format(
      belt_area_km2,
      digits = 12,
      scientific = FALSE
    ),
    format(
      raster_area_km2,
      digits = 12,
      scientific = FALSE
    )
  )
)

# ------------------------------------------------------------
# 5. Read declared-filter GBIF occurrence archive
# ------------------------------------------------------------

precoast <- readr::read_csv(
  precoast_file,
  show_col_types = FALSE
)

required_occ_cols <- c(
  "decimal_longitude",
  "decimal_latitude",
  "gbif_id"
)

missing_occ_cols <- setdiff(
  required_occ_cols,
  names(precoast)
)

if (length(missing_occ_cols) > 0) {
  stop(
    "Pre-coastal occurrence file lacks: ",
    paste(
      missing_occ_cols,
      collapse = ", "
    )
  )
}

precoast$decimal_longitude <- suppressWarnings(
  as.numeric(
    precoast$decimal_longitude
  )
)

precoast$decimal_latitude <- suppressWarnings(
  as.numeric(
    precoast$decimal_latitude
  )
)

precoast$gbif_id <- as.character(
  precoast$gbif_id
)

eu_precoast <- precoast %>%
  filter(
    is.finite(
      decimal_longitude
    ),
    is.finite(
      decimal_latitude
    ),
    decimal_longitude >= eu_ext_vec[1],
    decimal_longitude <= eu_ext_vec[2],
    decimal_latitude >= eu_ext_vec[3],
    decimal_latitude <= eu_ext_vec[4]
  )

# ------------------------------------------------------------
# 6. Occurrence-to-mask assignment in EPSG:3035
# ------------------------------------------------------------

assign_occurrences_to_mask <- function(
  df,
  mask_rast,
  valid_df,
  valid_m,
  bbox_template,
  max_snap_m = 1500,
  deduplicate = TRUE
) {

  if (nrow(df) == 0) {
    return(
      df
    )
  }

  df <- as.data.frame(
    df
  )

  df$.source_row <- seq_len(
    nrow(df)
  )

  df$orig_decimal_longitude <-
    df$decimal_longitude

  df$orig_decimal_latitude <-
    df$decimal_latitude

  xy <- as.matrix(
    df[
      ,
      c(
        "decimal_longitude",
        "decimal_latitude"
      )
    ]
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
    !is.na(
      mask_vals[
        cell_guess
      ]
    )

  out_list <- list()

  if (any(exact_ok)) {

    exact <- df[
      exact_ok,
      ,
      drop = FALSE
    ]

    exact$cell_id <-
      as.integer(
        cell_guess[
          exact_ok
        ]
      )

    exact_xy <- terra::xyFromCell(
      mask_rast,
      exact$cell_id
    )

    exact$decimal_longitude <-
      exact_xy[
        ,
        1
      ]

    exact$decimal_latitude <-
      exact_xy[
        ,
        2
      ]

    exact$snap_distance_m <- 0
    exact$snap_type <- "exact_valid_cell"

    out_list[["exact"]] <- exact
  }

  remaining <- df[
    !exact_ok,
    ,
    drop = FALSE
  ]

  if (nrow(remaining) > 0) {

    rem_sf <- sf::st_as_sf(
      remaining,
      coords = c(
        "decimal_longitude",
        "decimal_latitude"
      ),
      crs = 4326,
      remove = FALSE
    )

    rem_m <- sf::st_transform(
      rem_sf,
      distance_crs
    )

    idx <- sf::st_nearest_feature(
      rem_m,
      valid_m
    )

    d <- as.numeric(
      sf::st_distance(
        rem_m,
        valid_m[
          idx,
        ],
        by_element = TRUE
      )
    )

    keep <- is.finite(d) &
      d <= max_snap_m

    if (any(keep)) {

      snapped <- remaining[
        keep,
        ,
        drop = FALSE
      ]

      nearest <- idx[
        keep
      ]

      snapped$cell_id <-
        valid_df$cell_id[
          nearest
        ]

      snapped$decimal_longitude <-
        valid_df$decimal_longitude[
          nearest
        ]

      snapped$decimal_latitude <-
        valid_df$decimal_latitude[
          nearest
        ]

      snapped$snap_distance_m <-
        d[
          keep
        ]

      snapped$snap_type <-
        "nearest_valid_cell"

      out_list[["snapped"]] <- snapped
    }
  }

  if (length(out_list) == 0) {
    return(
      df[
        0,
        ,
        drop = FALSE
      ]
    )
  }

  assigned <- bind_rows(
    out_list
  )

  assigned$bbox_cell_id <- as.integer(
    terra::cellFromXY(
      bbox_template,
      as.matrix(
        assigned[
          ,
          c(
            "decimal_longitude",
            "decimal_latitude"
          )
        ]
      )
    )
  )

  if (deduplicate) {

    assigned <- assigned %>%
      arrange(
        cell_id,
        snap_distance_m,
        gbif_id
      ) %>%
      group_by(
        cell_id
      ) %>%
      slice(1) %>%
      ungroup()
  }

  assigned$.source_row <- NULL

  assigned
}

gbif_final7_all <- assign_occurrences_to_mask(
  eu_precoast,
  final7_mask,
  valid_df,
  valid_m,
  clim_bbox[[1]],
  max_snap_m = 1500,
  deduplicate = FALSE
)

gbif_final7 <- gbif_final7_all %>%
  arrange(
    cell_id,
    snap_distance_m,
    gbif_id
  ) %>%
  group_by(
    cell_id
  ) %>%
  slice(1) %>%
  ungroup()

gbif_final7$pa <- 1L

# ------------------------------------------------------------
# 7. Compare with REV1 on the common CHELSA base grid
# ------------------------------------------------------------

if (!file.exists(rev1_rdata_general)) {
  stop(
    "REV1 RData not found: ",
    rev1_rdata_general
  )
}

rev1 <- new.env(
  parent = baseenv()
)

load(
  rev1_rdata_general,
  envir = rev1
)

if (!exists(
  "occ",
  envir = rev1,
  inherits = FALSE
)) {
  stop(
    "REV1 final invaded occurrence object 'occ' not found."
  )
}

rev1_occ <- as.data.frame(
  get(
    "occ",
    envir = rev1,
    inherits = FALSE
  )
)

required_rev1_cols <- c(
  "gbif_id",
  "decimal_longitude",
  "decimal_latitude"
)

missing_rev1_cols <- setdiff(
  required_rev1_cols,
  names(rev1_occ)
)

if (length(missing_rev1_cols) > 0) {
  stop(
    "REV1 occurrence object lacks: ",
    paste(
      missing_rev1_cols,
      collapse = ", "
    )
  )
}

rev1_occ$gbif_id <- as.character(
  rev1_occ$gbif_id
)

rev1_occ$bbox_cell_id <- as.integer(
  terra::cellFromXY(
    clim_bbox[[1]],
    as.matrix(
      rev1_occ[
        ,
        c(
          "decimal_longitude",
          "decimal_latitude"
        )
      ]
    )
  )
)

rev1_grid_cells <- unique(
  rev1_occ$bbox_cell_id[
    !is.na(
      rev1_occ$bbox_cell_id
    )
  ]
)

rev2_grid_cells <- unique(
  gbif_final7$bbox_cell_id[
    !is.na(
      gbif_final7$bbox_cell_id
    )
  ]
)

cell_compare_summary <- tibble(
  metric = c(
    "REV1 invaded cells on common CHELSA grid",
    "REV2 corrected FINAL7 GBIF cells",
    "Common CHELSA-grid cells shared",
    "New REV2 CHELSA-grid cells absent REV1",
    "REV1 CHELSA-grid cells absent REV2"
  ),
  value = c(
    length(
      rev1_grid_cells
    ),
    length(
      rev2_grid_cells
    ),
    length(
      intersect(
        rev1_grid_cells,
        rev2_grid_cells
      )
    ),
    length(
      setdiff(
        rev2_grid_cells,
        rev1_grid_cells
      )
    ),
    length(
      setdiff(
        rev1_grid_cells,
        rev2_grid_cells
      )
    )
  )
)

absent_rev1_grid_cells <- setdiff(
  rev1_grid_cells,
  rev2_grid_cells
)

rev1_absent_trace <- rev1_occ %>%
  filter(
    bbox_cell_id %in%
      absent_rev1_grid_cells
  ) %>%
  mutate(
    in_declared_precoast =
      gbif_id %in%
      precoast$gbif_id,
    in_europe_precoast =
      gbif_id %in%
      eu_precoast$gbif_id
  )

rev2_all_map <- gbif_final7_all %>%
  transmute(
    gbif_id = as.character(
      gbif_id
    ),
    rev2_bbox_cell_id_same_record =
      bbox_cell_id,
    rev2_mask_cell_id_same_record =
      cell_id,
    rev2_snap_distance_m_same_record =
      snap_distance_m,
    rev2_snap_type_same_record =
      snap_type
  )

rev1_absent_trace <- rev1_absent_trace %>%
  left_join(
    rev2_all_map,
    by = "gbif_id"
  ) %>%
  mutate(
    same_record_reassigned =
      !is.na(
        rev2_bbox_cell_id_same_record
      ) &
      rev2_bbox_cell_id_same_record !=
        bbox_cell_id,
    same_record_same_grid_cell =
      !is.na(
        rev2_bbox_cell_id_same_record
      ) &
      rev2_bbox_cell_id_same_record ==
        bbox_cell_id
  )

rev1_absent_trace$replacement_gbif_id_same_coord <-
  NA_character_

rev1_absent_trace$replacement_rev2_bbox_cell_id_same_coord <-
  NA_integer_

if (
  "coord_key" %in% names(rev1_absent_trace) &&
  "coord_key" %in% names(precoast)
) {

  precoast_coord_map <- precoast %>%
    filter(
      !is.na(
        coord_key
      ),
      coord_key != ""
    ) %>%
    transmute(
      coord_key,
      replacement_gbif_id_same_coord =
        as.character(
          gbif_id
        )
    ) %>%
    distinct(
      coord_key,
      .keep_all = TRUE
    )

  replacement_map <- rev1_absent_trace %>%
    select(
      gbif_id,
      coord_key
    ) %>%
    left_join(
      precoast_coord_map,
      by = "coord_key"
    ) %>%
    left_join(
      rev2_all_map %>%
        rename(
          replacement_gbif_id_same_coord =
            gbif_id,
          replacement_rev2_bbox_cell_id_same_coord =
            rev2_bbox_cell_id_same_record
        ) %>%
        select(
          replacement_gbif_id_same_coord,
          replacement_rev2_bbox_cell_id_same_coord
        ),
      by = "replacement_gbif_id_same_coord"
    )

  rev1_absent_trace <-
    rev1_absent_trace %>%
    select(
      -replacement_gbif_id_same_coord,
      -replacement_rev2_bbox_cell_id_same_coord
    ) %>%
    left_join(
      replacement_map,
      by = c(
        "gbif_id",
        "coord_key"
      )
    )
}

if (exists(
  "flagged",
  envir = rev1,
  inherits = FALSE
)) {

  flagged_rev1 <- as.data.frame(
    get(
      "flagged",
      envir = rev1,
      inherits = FALSE
    )
  )

  if ("gbif_id" %in%
      names(flagged_rev1)) {

    flagged_rev1$gbif_id <- as.character(
      flagged_rev1$gbif_id
    )

    if (!"coordinate_uncertainty_in_meters" %in%
        names(flagged_rev1)) {
      flagged_rev1$coordinate_uncertainty_in_meters <-
        NA_real_
    }

    if (!"basis_of_record" %in%
        names(flagged_rev1)) {
      flagged_rev1$basis_of_record <-
        NA_character_
    }

    if (!"cc_pass" %in%
        names(flagged_rev1)) {
      flagged_rev1$cc_pass <-
        NA
    }

    allowed_basis <- c(
      "HUMAN_OBSERVATION",
      "OBSERVATION",
      "PRESERVED_SPECIMEN"
    )

    flagged_status <- flagged_rev1 %>%
      mutate(
        basis_norm =
          toupper(
            trimws(
              as.character(
                basis_of_record
              )
            )
          ),
        uncertainty_num =
          suppressWarnings(
            as.numeric(
              coordinate_uncertainty_in_meters
            )
          ),
        pass_basis =
          basis_norm %in%
          allowed_basis,
        pass_cc =
          !is.na(
            cc_pass
          ) &
          as.logical(
            cc_pass
          ),
        pass_uncertainty_5km =
          is.na(
            uncertainty_num
          ) |
          uncertainty_num <= 5000
      ) %>%
      select(
        gbif_id,
        basis_norm,
        uncertainty_num,
        pass_basis,
        pass_cc,
        pass_uncertainty_5km
      )

    rev1_absent_trace <- rev1_absent_trace %>%
      left_join(
        flagged_status,
        by = "gbif_id"
      ) %>%
      mutate(
        likely_reason =
          case_when(
            is.na(
              pass_basis
            ) ~
              "record_not_found_in_flagged_archive",
            !pass_basis ~
              "fails_declared_basisOfRecord_filter",
            !pass_cc ~
              "fails_declared_CoordinateCleaner_filter",
            !pass_uncertainty_5km ~
              "fails_declared_uncertainty_5km_filter",
            !in_declared_precoast &
              !is.na(
                replacement_gbif_id_same_coord
              ) ~
              "old_REV1_record_removed_by_coordinate_dedup_but_same_coordinate_retained",
            !in_declared_precoast ~
              "passes_declared_filters_but_removed_by_other_precoast_step",
            !in_europe_precoast ~
              "outside_Europe_bbox",
            same_record_reassigned ~
              "same_GBIF_record_reassigned_under_corrected_true1km_domain",
            same_record_same_grid_cell ~
              "same_GBIF_record_same_grid_cell_but_not_selected_as_representative",
            TRUE ~
              "unresolved_after_corrected_domain_trace"
          )
      )
  }
}

if (!"likely_reason" %in%
    names(rev1_absent_trace)) {
  rev1_absent_trace$likely_reason <-
    "unresolved_no_flagged_status"
}

rev1_absent_reason_summary <- rev1_absent_trace %>%
  count(
    likely_reason,
    name = "n_REV1_cells",
    sort = TRUE
  )

# ------------------------------------------------------------
# 8. GBIF country summaries
# ------------------------------------------------------------

if (!"country_code" %in%
    names(gbif_final7)) {
  gbif_final7$country_code <-
    NA_character_
}

gbif_final7$country_code <-
  toupper(
    trimws(
      as.character(
        gbif_final7$country_code
      )
    )
  )

country_summary <- gbif_final7 %>%
  filter(
    !is.na(
      country_code
    ),
    nzchar(
      country_code
    )
  ) %>%
  count(
    country_code,
    name = "n_REV2_FINAL7_cells",
    sort = TRUE
  )

target_codes <- c(
  "GB", "FR", "IT",
  "ES", "PT",
  "SI", "HR", "BA", "ME", "AL", "GR",
  "TR", "BG", "RO",
  "CY", "MT",
  "DE", "DK", "NO", "SE", "FI"
)

target_country_summary <- tidyr::complete(
  country_summary %>%
    filter(
      country_code %in%
        target_codes
    ),
  country_code = target_codes,
  fill = list(
    n_REV2_FINAL7_cells = 0
  )
) %>%
  arrange(
    country_code
  )

# ------------------------------------------------------------
# 9. Reconcile Sara/Ribotta against corrected GBIF set
# ------------------------------------------------------------

step03 <- readRDS(
  step03_file
)

sara_points <- step03$sara_points

required_sara <- c(
  "sara_id",
  "sara_lon",
  "sara_lat",
  "country_code",
  "country_name"
)

missing_sara <- setdiff(
  required_sara,
  names(sara_points)
)

if (length(missing_sara) > 0) {
  stop(
    "Step-03 Sara table lacks: ",
    paste(
      missing_sara,
      collapse = ", "
    )
  )
}

sara_input <- sara_points %>%
  select(
    sara_id,
    sara_lon,
    sara_lat,
    country_code,
    country_name
  ) %>%
  mutate(
    decimal_longitude =
      sara_lon,
    decimal_latitude =
      sara_lat,
    gbif_id =
      paste0(
        "SARA_",
        sara_id
      ),
    .source_row =
      row_number()
  )

assign_points_no_dedup <- function(
  df,
  mask_rast,
  valid_df,
  valid_m,
  max_snap_m = 1500
) {

  df <- as.data.frame(
    df
  )

  xy <- as.matrix(
    df[
      ,
      c(
        "decimal_longitude",
        "decimal_latitude"
      )
    ]
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
    !is.na(
      cell_guess
    ) &
    !is.na(
      mask_vals[
        cell_guess
      ]
    )

  df$assigned_cell_id <-
    NA_integer_

  df$assigned_longitude <-
    NA_real_

  df$assigned_latitude <-
    NA_real_

  df$snap_distance_m <-
    NA_real_

  df$snap_type <-
    "outside_1500m"

  if (any(exact_ok)) {

    exact_cells <- as.integer(
      cell_guess[
        exact_ok
      ]
    )

    exact_xy <- terra::xyFromCell(
      mask_rast,
      exact_cells
    )

    df$assigned_cell_id[
      exact_ok
    ] <- exact_cells

    df$assigned_longitude[
      exact_ok
    ] <- exact_xy[
      ,
      1
    ]

    df$assigned_latitude[
      exact_ok
    ] <- exact_xy[
      ,
      2
    ]

    df$snap_distance_m[
      exact_ok
    ] <- 0

    df$snap_type[
      exact_ok
    ] <- "exact_valid_cell"
  }

  rem <- which(
    !exact_ok
  )

  if (length(rem) > 0) {

    rem_sf <- sf::st_as_sf(
      df[
        rem,
        ,
        drop = FALSE
      ],
      coords = c(
        "decimal_longitude",
        "decimal_latitude"
      ),
      crs = 4326,
      remove = FALSE
    )

    rem_m <- sf::st_transform(
      rem_sf,
      distance_crs
    )

    idx <- sf::st_nearest_feature(
      rem_m,
      valid_m
    )

    d <- as.numeric(
      sf::st_distance(
        rem_m,
        valid_m[
          idx,
        ],
        by_element = TRUE
      )
    )

    keep <- is.finite(d) &
      d <= max_snap_m

    if (any(keep)) {

      rows <- rem[
        keep
      ]

      nearest <- idx[
        keep
      ]

      df$assigned_cell_id[
        rows
      ] <- valid_df$cell_id[
        nearest
      ]

      df$assigned_longitude[
        rows
      ] <- valid_df$decimal_longitude[
        nearest
      ]

      df$assigned_latitude[
        rows
      ] <- valid_df$decimal_latitude[
        nearest
      ]

      df$snap_distance_m[
        rows
      ] <- d[
        keep
      ]

      df$snap_type[
        rows
      ] <- "nearest_valid_cell"
    }

    if (any(!keep)) {

      rows <- rem[
        !keep
      ]

      df$snap_distance_m[
        rows
      ] <- d[
        !keep
      ]
    }
  }

  df$retained_FINAL7 <-
    !is.na(
      df$assigned_cell_id
    )

  df
}

sara_point_audit <- assign_points_no_dedup(
  sara_input,
  final7_mask,
  valid_df,
  valid_m,
  max_snap_m = 1500
)

sara_point_audit$cell_already_in_GBIF_FINAL7 <-
  sara_point_audit$assigned_cell_id %in%
  gbif_final7$cell_id

gbif_original_sf <- sf::st_as_sf(
  gbif_final7,
  coords = c(
    "orig_decimal_longitude",
    "orig_decimal_latitude"
  ),
  crs = 4326,
  remove = FALSE
)

sara_sf <- sf::st_as_sf(
  sara_point_audit,
  coords = c(
    "sara_lon",
    "sara_lat"
  ),
  crs = 4326,
  remove = FALSE
)

gbif_original_m <- sf::st_transform(
  gbif_original_sf,
  distance_crs
)

sara_m <- sf::st_transform(
  sara_sf,
  distance_crs
)

nearest_retained_idx <- sf::st_nearest_feature(
  sara_m,
  gbif_original_m
)

nearest_retained_dist <- as.numeric(
  sf::st_distance(
    sara_m,
    gbif_original_m[
      nearest_retained_idx,
    ],
    by_element = TRUE
  )
)

sara_point_audit$nearest_RETAINED_GBIF_distance_m <-
  nearest_retained_dist

sara_point_audit$nearest_RETAINED_GBIF_id <-
  as.character(
    gbif_final7$gbif_id[
      nearest_retained_idx
    ]
  )

declared_gbif_sf <- sf::st_as_sf(
  eu_precoast,
  coords = c(
    "decimal_longitude",
    "decimal_latitude"
  ),
  crs = 4326,
  remove = FALSE
)

declared_gbif_m <- sf::st_transform(
  declared_gbif_sf,
  distance_crs
)

nearest_declared_idx <- sf::st_nearest_feature(
  sara_m,
  declared_gbif_m
)

nearest_declared_dist <- as.numeric(
  sf::st_distance(
    sara_m,
    declared_gbif_m[
      nearest_declared_idx,
    ],
    by_element = TRUE
  )
)

sara_point_audit$nearest_DECLARED_GBIF_distance_m <-
  nearest_declared_dist

sara_point_audit$nearest_DECLARED_GBIF_id <-
  as.character(
    eu_precoast$gbif_id[
      nearest_declared_idx
    ]
  )

sara_point_audit$candidate_independent_field_point <-
  sara_point_audit$retained_FINAL7 &
  !sara_point_audit$cell_already_in_GBIF_FINAL7 &
  is.finite(
    sara_point_audit$nearest_DECLARED_GBIF_distance_m
  ) &
  sara_point_audit$nearest_DECLARED_GBIF_distance_m >
  1000

candidate_independent_cells <- sara_point_audit %>%
  filter(
    candidate_independent_field_point
  ) %>%
  arrange(
    assigned_cell_id,
    nearest_DECLARED_GBIF_distance_m,
    sara_id
  ) %>%
  group_by(
    assigned_cell_id
  ) %>%
  slice(1) %>%
  ungroup()

italy_detail <- sara_point_audit %>%
  filter(
    country_code == "IT"
  ) %>%
  arrange(
    desc(
      candidate_independent_field_point
    ),
    nearest_DECLARED_GBIF_distance_m,
    sara_id
  )

sara_country_summary <- sara_point_audit %>%
  mutate(
    country_code =
      ifelse(
        is.na(
          country_code
        ) |
          country_code == "",
        "UNRESOLVED",
        country_code
      )
  ) %>%
  group_by(
    country_code
  ) %>%
  summarise(
    n_Sara_points = n(),
    n_retained_FINAL7 =
      sum(
        retained_FINAL7,
        na.rm = TRUE
      ),
    n_cells_already_GBIF =
      sum(
        retained_FINAL7 &
          cell_already_in_GBIF_FINAL7,
        na.rm = TRUE
      ),
    n_candidate_independent_points =
      sum(
        candidate_independent_field_point,
        na.rm = TRUE
      ),
    .groups = "drop"
  ) %>%
  arrange(
    desc(
      n_Sara_points
    )
  )

# ------------------------------------------------------------
# 10. Save canonical outputs
# ------------------------------------------------------------

readr::write_csv(
  gbif_final7,
  file.path(
    work_dir,
    "occurrences_REV2_FINAL7_GBIF.csv"
  )
)

readr::write_csv(
  cell_compare_summary,
  file.path(
    audit_dir,
    "05A_REV1_vs_REV2_cell_summary.csv"
  )
)

readr::write_csv(
  country_summary,
  file.path(
    audit_dir,
    "05B_REV2_FINAL7_country_counts.csv"
  )
)

readr::write_csv(
  target_country_summary,
  file.path(
    audit_dir,
    "05C_target_country_counts.csv"
  )
)

readr::write_csv(
  sara_point_audit,
  file.path(
    audit_dir,
    "05D_Sara_vs_REV2_FINAL7_point_audit.csv"
  )
)

readr::write_csv(
  sara_country_summary,
  file.path(
    audit_dir,
    "05E_Sara_vs_REV2_country_summary.csv"
  )
)

readr::write_csv(
  candidate_independent_cells,
  file.path(
    audit_dir,
    "05F_candidate_independent_Sara_cells.csv"
  )
)

readr::write_csv(
  italy_detail,
  file.path(
    audit_dir,
    "05G_Italy_Sara_point_detail.csv"
  )
)

summary_tbl <- tibble(
  metric = c(
    "Declared-filter records before Europe restriction",
    "Declared-filter records inside Europe bbox",
    "REV1 invaded cells on common CHELSA grid",
    "REV2 corrected FINAL7 GBIF cells",
    "New REV2 CHELSA-grid cells vs REV1",
    "REV1 CHELSA-grid cells absent REV2",
    "Sara points retained under corrected FINAL7",
    "Sara points in a corrected REV2 GBIF cell",
    "Sara candidate independent field points",
    "Sara candidate independent field cells",
    "Italian Sara points",
    "Italian Sara points retained under corrected FINAL7",
    "Italian candidate independent field points"
  ),
  value = c(
    nrow(
      precoast
    ),
    nrow(
      eu_precoast
    ),
    length(
      rev1_grid_cells
    ),
    nrow(
      gbif_final7
    ),
    length(
      setdiff(
        rev2_grid_cells,
        rev1_grid_cells
      )
    ),
    length(
      setdiff(
        rev1_grid_cells,
        rev2_grid_cells
      )
    ),
    sum(
      sara_point_audit$retained_FINAL7,
      na.rm = TRUE
    ),
    sum(
      sara_point_audit$retained_FINAL7 &
        sara_point_audit$cell_already_in_GBIF_FINAL7,
      na.rm = TRUE
    ),
    sum(
      sara_point_audit$candidate_independent_field_point,
      na.rm = TRUE
    ),
    nrow(
      candidate_independent_cells
    ),
    sum(
      sara_point_audit$country_code == "IT",
      na.rm = TRUE
    ),
    sum(
      sara_point_audit$country_code == "IT" &
        sara_point_audit$retained_FINAL7,
      na.rm = TRUE
    ),
    sum(
      sara_point_audit$country_code == "IT" &
        sara_point_audit$candidate_independent_field_point,
      na.rm = TRUE
    )
  )
)

readr::write_csv(
  summary_tbl,
  file.path(
    audit_dir,
    "05H_occurrence_rebuild_summary.csv"
  )
)

readr::write_csv(
  rev1_absent_trace,
  file.path(
    audit_dir,
    "05I_REV1_cells_absent_REV2_trace.csv"
  )
)

readr::write_csv(
  rev1_absent_reason_summary,
  file.path(
    audit_dir,
    "05J_REV1_cells_absent_REV2_reason_summary.csv"
  )
)

readr::write_csv(
  domain_summary,
  file.path(
    audit_dir,
    "05K_corrected_coastal_domain_summary.csv"
  )
)

saveRDS(
  list(
    final7_vars = final7_vars,
    n_final7_cells = n_final7_cells,
    gbif_final7 = gbif_final7,
    gbif_final7_all = gbif_final7_all,
    cell_compare_summary = cell_compare_summary,
    country_summary = country_summary,
    target_country_summary = target_country_summary,
    sara_point_audit = sara_point_audit,
    sara_country_summary = sara_country_summary,
    candidate_independent_cells = candidate_independent_cells,
    italy_detail = italy_detail,
    rev1_absent_trace = rev1_absent_trace,
    rev1_absent_reason_summary = rev1_absent_reason_summary,
    summary = summary_tbl,
    domain_summary = domain_summary,
    bbox_file = bbox_file,
    grid_file = mask_file,
    mask_file = mask_file,
    coastal_belt_file = belt_file,
    coastal_buffer_m = coastal_buffer_m,
    distance_crs = distance_crs,
    raster_inclusion_rule =
      "cell_center_inside_terrestrial_1km_belt"
  ),
  file.path(
    work_dir,
    "05_REV2_occurrences_FINAL7_VALID.rds"
  )
)

# ------------------------------------------------------------
# 11. Update canonical workflow manifest
# ------------------------------------------------------------

manifest_file <- file.path(
  audit_dir,
  "REV2_CANONICAL_WORKFLOW.txt"
)

canonical_workflow <- c(
  "ROSA RUGOSA REV2 - CANONICAL WORKFLOW",
  "",
  "VALID SCRIPTS",
  "00_setup_REV2.R",
  "01_occurrence_audit_REV2.R",
  "02_model_input_audit_REV2.R",
  "03_reconcile_Ribotta_with_REV1.R",
  "04_audit_REV1_coastal_mask_support.R",
  "05_rebuild_REV2_occurrences_FINAL7.R",
  "06_freeze_REV2_occurrence_dataset.R",
  "07_build_REV2_background_FINAL7.R",
  "",
  "IMPORTANT",
  "Steps 03-04 audit the REV1 workflow and diagnose the former ALL10/FCF coastal-mask problem.",
  "Step 05 rebuilds the REV2 European domain from the unmasked current-climate raster using a 1,000-m terrestrial coastal belt in EPSG:3035.",
  "Raster cells are included when their centres lie inside the corrected coastal belt; the former EPSG:3857 belt is not used.",
  "Occurrence snapping up to 1,500 m is also evaluated in EPSG:3035.",
  "Sara/Ribotta data are audited in step 05 but are not added to calibration there.",
  "No model is fitted in steps 00-07."
)

writeLines(
  canonical_workflow,
  manifest_file
)

# ------------------------------------------------------------
# 12. Console report
# ------------------------------------------------------------

cat(
  "\n============================================================\n"
)

cat(
  "REV2 - CORRECTED TRUE-1KM FINAL7 OCCURRENCE REBUILD\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "Corrected coastal-domain summary:\n"
)

print(
  as.data.frame(
    domain_summary
  ),
  row.names = FALSE
)

cat(
  "\nCore occurrence summary:\n"
)

print(
  as.data.frame(
    summary_tbl
  ),
  row.names = FALSE
)

cat(
  "\nREV1 vs corrected REV2 common-grid comparison:\n"
)

print(
  as.data.frame(
    cell_compare_summary
  ),
  row.names = FALSE
)

cat(
  "\nTarget-country corrected REV2 counts:\n"
)

print(
  as.data.frame(
    target_country_summary
  ),
  row.names = FALSE
)

cat(
  "\nREV1 cells absent from corrected REV2 - reason summary:\n"
)

print(
  as.data.frame(
    rev1_absent_reason_summary
  ),
  row.names = FALSE
)

cat(
  "\nItalian Sara/Ribotta detail:\n"
)

print(
  as.data.frame(
    italy_detail %>%
      select(
        sara_id,
        sara_lon,
        sara_lat,
        retained_FINAL7,
        assigned_cell_id,
        cell_already_in_GBIF_FINAL7,
        nearest_RETAINED_GBIF_distance_m,
        nearest_RETAINED_GBIF_id,
        nearest_DECLARED_GBIF_distance_m,
        nearest_DECLARED_GBIF_id,
        candidate_independent_field_point
      )
  ),
  row.names = FALSE
)

cat(
  "\nNo Sara/Ribotta record has been added to calibration in step 05.\n"
)

cat(
  "No model has been fitted. No REV1 file has been modified.\n"
)

cat(
  "============================================================\n"
)
