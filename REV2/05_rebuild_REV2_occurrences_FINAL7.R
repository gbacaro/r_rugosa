# ============================================================
# ROSA RUGOSA - REVISION 2
# 05_rebuild_REV2_occurrences_FINAL7.R
#
# Purpose
# -------
# Rebuild the European invaded occurrence set using:
#   - the occurrence-quality filters explicitly declared in REV1;
#   - the geographic coastal raster;
#   - ONLY the seven predictors retained in the final SDM;
#   - the original <=1500 m snapping rule;
#   - one occurrence per raster cell.
#
# Sara/Ribotta data are NOT added to calibration here. They are
# reconciled against the rebuilt GBIF set to identify genuinely
# independent candidate field observations for later review.
#
# No model is fitted.
# No REV1 file is modified.
# ============================================================

rm(list = ls())
gc()

required_pkgs <- c(
  "terra", "sf", "dplyr", "readr", "tibble", "tidyr"
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

old_05 <- c(
  list.files(
    audit_dir,
    pattern = "^05[A-Z]_",
    full.names = TRUE
  ),
  file.path(
    work_dir,
    "05_REV2_occurrences_FINAL7_VALID.rds"
  ),
  file.path(
    work_dir,
    "occurrences_REV2_FINAL7_GBIF.csv"
  )
)

old_05 <- old_05[
  file.exists(old_05)
]

if (length(old_05) > 0) {
  unlink(
    old_05,
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

# Search the official current coastal stack.
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
    out[file.exists(out)]
  )
}

coast_candidates <- find_exact_file(
  search_roots,
  "current_europe_coast1km_selectedvars.tif"
)

if (length(coast_candidates) == 0) {
  stop(
    "Could not find current_europe_coast1km_selectedvars.tif."
  )
}

preferred <- grep(
  "data_raw[/\\\\]climate_chelsa[/\\\\]current",
  coast_candidates,
  ignore.case = TRUE
)

if (length(preferred) > 0) {
  coast_file <- coast_candidates[
    preferred[1]
  ]
} else {
  coast_file <- coast_candidates[1]
}

# ------------------------------------------------------------
# 3. Build canonical FINAL7 mask
# ------------------------------------------------------------

clim <- terra::rast(
  coast_file
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
  names(clim)
)

if (length(missing_final7) > 0) {
  stop(
    "Current coastal raster is missing FINAL7 predictors: ",
    paste(
      missing_final7,
      collapse = ", "
    )
  )
}

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

final7_mask <- complete_mask(
  clim,
  final7_vars
)

n_final7_cells <- as.integer(
  terra::global(
    final7_mask,
    "sum",
    na.rm = TRUE
  )[1, 1]
)

expected_final7_cells <- 210748L

if (n_final7_cells != expected_final7_cells) {
  stop(
    "FINAL7 mask contains ",
    n_final7_cells,
    " cells instead of the audited 210,748. ",
    "Do not continue."
  )
}

# ------------------------------------------------------------
# 4. Read declared-filter GBIF occurrence archive
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

# Same European extent used by the REV1 pipeline.
eu_ext_vec <- c(
  -15,
  45,
  30,
  72
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
# 5. Canonical occurrence-to-mask assignment
# ------------------------------------------------------------

assign_occurrences_to_mask <- function(
  df,
  mask_rast,
  max_snap_m = 1500
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

  vals <- terra::values(
    mask_rast,
    mat = FALSE
  )

  exact_ok <-
    !is.na(cell_guess) &
    !is.na(
      vals[cell_guess]
    )

  out_list <- list()

  if (any(exact_ok)) {

    exact <- df[
      exact_ok,
      ,
      drop = FALSE
    ]

    exact$cell_id <-
      cell_guess[
        exact_ok
      ]

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

    exact$snap_type <-
      "exact_valid_cell"

    out_list[["exact"]] <-
      exact
  }

  remaining <- df[
    !exact_ok,
    ,
    drop = FALSE
  ]

  if (nrow(remaining) > 0) {

    valid_df <- terra::as.data.frame(
      mask_rast,
      xy = TRUE,
      cells = TRUE,
      na.rm = TRUE
    )[
      ,
      c(
        "cell",
        "x",
        "y"
      )
    ]

    valid_sf <- sf::st_as_sf(
      valid_df,
      coords = c(
        "x",
        "y"
      ),
      crs = 4326
    )

    rem_sf <- sf::st_as_sf(
      remaining,
      coords = c(
        "decimal_longitude",
        "decimal_latitude"
      ),
      crs = 4326,
      remove = FALSE
    )

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
        valid_df$cell[
          nearest
        ]

      snapped$decimal_longitude <-
        valid_df$x[
          nearest
        ]

      snapped$decimal_latitude <-
        valid_df$y[
          nearest
        ]

      snapped$snap_distance_m <-
        d[
          keep
        ]

      snapped$snap_type <-
        "nearest_valid_cell"

      out_list[["snapped"]] <-
        snapped
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

  # Deterministic one-record-per-cell aggregation.
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

  assigned$.source_row <- NULL

  assigned
}

gbif_final7 <- assign_occurrences_to_mask(
  eu_precoast,
  final7_mask,
  max_snap_m = 1500
)

gbif_final7$pa <- 1L

# ------------------------------------------------------------
# 6. Compare with REV1 final invaded occurrence set
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

if (!"cell_id" %in%
    names(rev1_occ)) {
  stop(
    "REV1 final invaded object lacks cell_id."
  )
}

rev1_cells <- unique(
  rev1_occ$cell_id
)

rev2_cells <- unique(
  gbif_final7$cell_id
)

cell_compare_summary <- tibble(
  metric = c(
    "REV1 final invaded cells",
    "REV2 FINAL7 GBIF cells",
    "Cells shared REV1 and REV2",
    "New REV2 cells absent REV1",
    "REV1 cells absent REV2"
  ),

  value = c(
    length(
      rev1_cells
    ),
    length(
      rev2_cells
    ),
    length(
      intersect(
        rev1_cells,
        rev2_cells
      )
    ),
    length(
      setdiff(
        rev2_cells,
        rev1_cells
      )
    ),
    length(
      setdiff(
        rev1_cells,
        rev2_cells
      )
    )
  )
)


# Trace the REV1 cells that disappear in the rebuilt REV2 set.
rev1_occ$gbif_id <- as.character(
  rev1_occ$gbif_id
)

absent_rev1_cells <- setdiff(
  rev1_cells,
  rev2_cells
)

rev1_absent_trace <- rev1_occ %>%
  filter(
    cell_id %in%
      absent_rev1_cells
  ) %>%
  mutate(
    in_declared_precoast =
      gbif_id %in%
      precoast$gbif_id,

    in_europe_precoast =
      gbif_id %in%
      eu_precoast$gbif_id
  )

# Reconstruct declared filter status for the missing REV1 records.
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

  if ("gbif_id" %in% names(flagged_rev1)) {

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
          !is.na(cc_pass) &
          as.logical(cc_pass),

        pass_uncertainty_5km =
          is.na(uncertainty_num) |
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
            is.na(pass_basis) ~
              "record_not_found_in_flagged_archive",

            !pass_basis ~
              "fails_declared_basisOfRecord_filter",

            !pass_cc ~
              "fails_declared_CoordinateCleaner_filter",

            !pass_uncertainty_5km ~
              "fails_declared_uncertainty_5km_filter",

            !in_declared_precoast ~
              "passes_declared_filters_but_removed_by_coordinate_dedup_or_other_precoast_step",

            !in_europe_precoast ~
              "outside_Europe_bbox",

            TRUE ~
              "record_enters_precoast_but_cell_not_retained_in_REV2"
          )
      )
  }
}

rev1_absent_reason_summary <- rev1_absent_trace %>%
  count(
    likely_reason,
    name = "n_REV1_cells",
    sort = TRUE
  )

# ------------------------------------------------------------
# 7. GBIF country summaries
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
# 8. Reconcile Sara/Ribotta against corrected FINAL7 GBIF set
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
      )
  )

sara_assigned <- assign_occurrences_to_mask(
  sara_input,
  final7_mask,
  max_snap_m = 1500
)

# assign_occurrences_to_mask aggregates one record per cell.
# For a point-level audit, assign all Sara points without deduplication.

assign_points_no_dedup <- function(
  df,
  mask_rast,
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

  vals <- terra::values(
    mask_rast,
    mat = FALSE
  )

  exact_ok <-
    !is.na(
      cell_guess
    ) &
    !is.na(
      vals[cell_guess]
    )

  df$assigned_cell_id <-
    NA_integer_

  df$snap_distance_m <-
    NA_real_

  df$snap_type <-
    "outside_1500m"

  if (any(exact_ok)) {

    df$assigned_cell_id[
      exact_ok
    ] <- cell_guess[
      exact_ok
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

    valid_df <- terra::as.data.frame(
      mask_rast,
      xy = TRUE,
      cells = TRUE,
      na.rm = TRUE
    )[
      ,
      c(
        "cell",
        "x",
        "y"
      )
    ]

    valid_sf <- sf::st_as_sf(
      valid_df,
      coords = c(
        "x",
        "y"
      ),
      crs = 4326
    )

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

      df$assigned_cell_id[
        rows
      ] <- valid_df$cell[
        idx[keep]
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
  max_snap_m = 1500
)

sara_point_audit$cell_already_in_GBIF_FINAL7 <-
  sara_point_audit$assigned_cell_id %in%
  rev2_cells

# Nearest retained FINAL7 GBIF occurrence using original coordinates.
# This is useful for model-cell context but is NOT used by itself to
# establish field-data independence.

gbif_sf <- sf::st_as_sf(
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

gbif_m <- sf::st_transform(
  gbif_sf,
  3035
)

sara_m <- sf::st_transform(
  sara_sf,
  3035
)

nearest_idx <- sf::st_nearest_feature(
  sara_m,
  gbif_m
)

nearest_dist <- as.numeric(
  sf::st_distance(
    sara_m,
    gbif_m[
      nearest_idx,
    ],
    by_element = TRUE
  )
)

sara_point_audit$nearest_RETAINED_GBIF_distance_m <-
  nearest_dist

sara_point_audit$nearest_RETAINED_GBIF_id <-
  as.character(
    gbif_final7$gbif_id[
      nearest_idx
    ]
  )

# Independence must be assessed against ALL GBIF records that passed
# the declared pre-coastal filters, before one-record-per-cell thinning.

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
  3035
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

# Candidate independent cells: one representative Sara point per cell.
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
    candidate_independent_field_point %>%
      desc(),
    nearest_DECLARED_GBIF_distance_m,
    sara_id
  )

# ------------------------------------------------------------
# 9. Sara summaries
# ------------------------------------------------------------

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

summary_tbl <- tibble(
  metric = c(
    "Declared-filter records before Europe restriction",
    "Declared-filter records inside Europe bbox",
    "REV1 final invaded cells",
    "REV2 FINAL7 GBIF cells",
    "REV2 cells new vs REV1",
    "REV1 cells absent REV2",
    "Sara points retained under FINAL7",
    "Sara points in a REV2 GBIF cell",
    "Sara candidate independent field points",
    "Sara candidate independent field cells",
    "Italian Sara points",
    "Italian Sara points retained under FINAL7",
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
      rev1_cells
    ),
    length(
      rev2_cells
    ),
    length(
      setdiff(
        rev2_cells,
        rev1_cells
      )
    ),
    length(
      setdiff(
        rev1_cells,
        rev2_cells
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

saveRDS(
  list(
    final7_vars = final7_vars,
    n_final7_cells = n_final7_cells,
    gbif_final7 = gbif_final7,
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
    coast_file = coast_file
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
  "",
  "SUPERSEDED / INVALIDATED",
  "03_reconcile_Ribotta_with_REV1_CORRECTED.R",
  "04_trace_missing_Ribotta_cells.R",
  "05_audit_Ribotta_geography_by_country.R",
  "",
  "IMPORTANT",
  "REV2 occurrence reconstruction uses the quality filters declared in the manuscript and the seven-predictor current coastal mask.",
  "Sara/Ribotta independence is assessed against all GBIF records passing the declared filters before spatial thinning; records are not yet added to model calibration.",
  "No model is refitted in steps 00-05."
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
  "REV2 - FINAL7 OCCURRENCE REBUILD COMPLETE\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "Core occurrence summary:\n"
)

print(
  as.data.frame(
    summary_tbl
  ),
  row.names = FALSE
)

cat(
  "\nREV1 vs REV2 cell comparison:\n"
)

print(
  as.data.frame(
    cell_compare_summary
  ),
  row.names = FALSE
)

cat(
  "\nTarget-country REV2 FINAL7 counts:\n"
)

print(
  as.data.frame(
    target_country_summary
  ),
  row.names = FALSE
)

cat(
  "\nREV1 cells absent from REV2 - reason summary:\n"
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
  "\nNo Sara/Ribotta record has been added to calibration.\n"
)

cat(
  "No model has been fitted. No REV1 file has been modified.\n"
)

cat(
  "============================================================\n"
)
