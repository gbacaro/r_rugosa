# ============================================================
# ROSA RUGOSA - REVISION 2
# 04_audit_REV1_coastal_mask_support.R
#
# Purpose
# -------
# Audit why the official REV1 raster-valid European coastal mask
# excludes large parts of France, Italy and the United Kingdom.
#
# The script reconstructs and compares:
#   A) geographic 1-km coastal belt
#   B) all-10-variable REV1 complete-case mask
#   C) final 7-predictor complete-case mask
#   D) per-variable current-climate spatial support
#
# It identifies which predictor(s) create current spatial gaps and
# quantifies their effects by country and for Sara/Ribotta points.
#
# No occurrence is added.
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

if (!exists("audit_dir")) audit_dir <- file.path(rev2_dir, "03_audit")
if (!exists("work_dir"))  work_dir  <- file.path(rev2_dir, "04_work")

dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

# Clean ONLY outputs belonging to canonical step 04.
old_04 <- c(
  list.files(
    audit_dir,
    pattern = "^04[A-Z]_",
    full.names = TRUE
  ),
  file.path(
    work_dir,
    "04_REV1_coastal_mask_support_VALID.rds"
  )
)

old_04 <- old_04[file.exists(old_04)]

if (length(old_04) > 0) {
  unlink(old_04, force = TRUE)
}

# ------------------------------------------------------------
# 2. Search official REV1 assets
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

find_exact_file <- function(roots, filename) {

  ans <- unlist(
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

  unique(
    ans[file.exists(ans)]
  )
}

choose_preferred <- function(
  files,
  preferred_pattern = NULL
) {

  if (length(files) == 0) {
    return(NA_character_)
  }

  if (!is.null(preferred_pattern)) {

    idx <- grep(
      preferred_pattern,
      files,
      ignore.case = TRUE
    )

    if (length(idx) > 0) {
      return(files[idx[1]])
    }
  }

  files[1]
}

belt_candidates <- find_exact_file(
  search_roots,
  "europe_coastal_belt_1km_geographic.rds"
)

mask10_candidates <- find_exact_file(
  search_roots,
  "europe_valid_coastal_cells_1km.rds"
)

coast_stack_candidates <- find_exact_file(
  search_roots,
  "current_europe_coast1km_selectedvars.tif"
)

bbox_stack_candidates <- find_exact_file(
  search_roots,
  "current_europe_bbox_selectedvars.tif"
)

belt_file <- choose_preferred(
  belt_candidates,
  "data_clean"
)

mask10_file <- choose_preferred(
  mask10_candidates,
  "data_clean"
)

coast_stack_file <- choose_preferred(
  coast_stack_candidates,
  "data_raw[/\\\\]climate_chelsa[/\\\\]current"
)

bbox_stack_file <- choose_preferred(
  bbox_stack_candidates,
  "data_raw[/\\\\]climate_chelsa[/\\\\]current"
)

asset_inventory <- tibble(
  asset = c(
    "geographic coastal belt",
    "REV1 all-10 complete-case mask",
    "current coastal climate stack",
    "current Europe bbox climate stack"
  ),
  path = c(
    belt_file,
    mask10_file,
    coast_stack_file,
    bbox_stack_file
  ),
  exists = c(
    !is.na(belt_file),
    !is.na(mask10_file),
    !is.na(coast_stack_file),
    !is.na(bbox_stack_file)
  )
)

write_csv(
  asset_inventory,
  file.path(
    audit_dir,
    "04A_asset_inventory.csv"
  )
)

print(
  asset_inventory,
  n = nrow(asset_inventory)
)

if (is.na(coast_stack_file)) {
  stop(
    "current_europe_coast1km_selectedvars.tif was not found. ",
    "This audit requires the current 10-layer coastal stack."
  )
}

if (is.na(bbox_stack_file)) {
  stop(
    "current_europe_bbox_selectedvars.tif was not found. ",
    "This audit requires the unmasked Europe bbox stack to distinguish ",
    "geographic coastal-mask loss from predictor-specific NA values."
  )
}

# ------------------------------------------------------------
# 3. Load current coastal and bbox climate stacks
# ------------------------------------------------------------

clim10 <- terra::rast(
  coast_stack_file
)

bbox10 <- terra::rast(
  bbox_stack_file
)

expected_vars <- c(
  "bio3", "bio8", "bio9", "bio10",
  "bio15", "bio18", "bio19",
  "fcf", "gsl", "gst"
)

missing_coast_vars <- setdiff(
  expected_vars,
  names(clim10)
)

missing_bbox_vars <- setdiff(
  expected_vars,
  names(bbox10)
)

if (length(missing_coast_vars) > 0) {
  stop(
    "Current coastal stack is missing: ",
    paste(missing_coast_vars, collapse = ", ")
  )
}

if (length(missing_bbox_vars) > 0) {
  stop(
    "Current bbox stack is missing: ",
    paste(missing_bbox_vars, collapse = ", ")
  )
}

clim10 <- clim10[[expected_vars]]
bbox10 <- bbox10[[expected_vars]]

if (!terra::compareGeom(
  clim10[[1]],
  bbox10[[1]],
  stopOnError = FALSE
)) {
  stop(
    "Current coastal and bbox stacks do not share raster geometry."
  )
}

final7_vars <- c(
  "bio3", "bio8", "bio9",
  "bio15", "bio18",
  "gsl", "gst"
)

# ------------------------------------------------------------
# 4. Optional saved REV1 mask
#
# The saved RDS is not required. When absent, the ALL10 mask is
# reconstructed directly from the official current coastal TIFF.
# ------------------------------------------------------------

mask10_saved <- NULL
n_mask10_saved <- NA_integer_

if (!is.na(mask10_file)) {

  mask10_saved <- readRDS(
    mask10_file
  )

  if (!inherits(mask10_saved, "SpatRaster")) {
    stop(
      "Saved REV1 mask exists but is not a SpatRaster."
    )
  }

  if (!terra::compareGeom(
    mask10_saved,
    clim10[[1]],
    stopOnError = FALSE
  )) {
    stop(
      "Saved REV1 mask and current coastal stack do not share raster geometry."
    )
  }

  n_mask10_saved <- as.integer(
    terra::global(
      mask10_saved,
      "sum",
      na.rm = TRUE
    )[1, 1]
  )
}

# ------------------------------------------------------------
# 5. Reconstruct complete-case masks
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

mask_all10_rebuilt <- complete_mask(
  clim10,
  expected_vars
)

mask_final7 <- complete_mask(
  clim10,
  final7_vars
)

n_all10_rebuilt <- as.integer(
  terra::global(
    mask_all10_rebuilt,
    "sum",
    na.rm = TRUE
  )[1, 1]
)

n_final7 <- as.integer(
  terra::global(
    mask_final7,
    "sum",
    na.rm = TRUE
  )[1, 1]
)

all10_rebuild_matches_saved <-
  if (is.na(n_mask10_saved)) {
    NA
  } else {
    n_all10_rebuilt == n_mask10_saved
  }

# The canonical REV1 analysis reported 79,579 ALL10 cells.
all10_matches_reported_REV1 <-
  n_all10_rebuilt == 79579L

if (!all10_matches_reported_REV1) {
  stop(
    "Rebuilt ALL10 mask contains ",
    n_all10_rebuilt,
    " cells instead of the REV1-reported 79,579. ",
    "Do not interpret downstream diagnostics."
  )
}

# ------------------------------------------------------------
# 6. Per-variable spatial support
# ------------------------------------------------------------

support_rows <- lapply(
  expected_vars,
  function(v) {

    valid <- !is.na(
      clim10[[v]]
    )

    n_valid <- as.integer(
      terra::global(
        valid,
        "sum",
        na.rm = TRUE
      )[1, 1]
    )

    tibble(
      variable = v,
      valid_cells = n_valid
    )
  }
)

variable_support <- bind_rows(
  support_rows
)

max_single_support <- max(
  variable_support$valid_cells,
  na.rm = TRUE
)

variable_support <- variable_support %>%
  mutate(
    pct_of_max_single_support =
      100 * valid_cells /
      max_single_support,

    cells_lost_vs_max_single =
      max_single_support -
      valid_cells
  ) %>%
  arrange(
    valid_cells
  )

write_csv(
  variable_support,
  file.path(
    audit_dir,
    "04B_current_variable_support.csv"
  )
)


# ------------------------------------------------------------
# 6B. BBOX vs COASTAL support by variable
# ------------------------------------------------------------

bbox_coast_support <- bind_rows(
  lapply(
    expected_vars,
    function(v) {

      bbox_valid <- !is.na(bbox10[[v]])
      coast_valid <- !is.na(clim10[[v]])

      n_bbox <- as.integer(
        terra::global(
          bbox_valid,
          "sum",
          na.rm = TRUE
        )[1, 1]
      )

      n_coast <- as.integer(
        terra::global(
          coast_valid,
          "sum",
          na.rm = TRUE
        )[1, 1]
      )

      tibble(
        variable = v,
        bbox_valid_cells = n_bbox,
        coastal_valid_cells = n_coast,
        coastal_pct_of_bbox =
          100 * n_coast / n_bbox
      )
    }
  )
)

write_csv(
  bbox_coast_support,
  file.path(
    audit_dir,
    "04B2_bbox_vs_coastal_variable_support.csv"
  )
)

# ------------------------------------------------------------
# 7. Incremental exclusion diagnostic
# ------------------------------------------------------------

# Start from final 7 predictors, then add each excluded/current-only
# predictor separately. This shows which variable causes spatial loss.

extra_vars <- setdiff(
  expected_vars,
  final7_vars
)

incremental_rows <- list()

incremental_rows[["FINAL7"]] <- tibble(
  predictor_set = "FINAL7",
  valid_cells = n_final7,
  cells_lost_vs_FINAL7 = 0L,
  pct_retained_vs_FINAL7 = 100
)

for (v in extra_vars) {

  m <- complete_mask(
    clim10,
    c(
      final7_vars,
      v
    )
  )

  n <- as.integer(
    terra::global(
      m,
      "sum",
      na.rm = TRUE
    )[1, 1]
  )

  incremental_rows[[v]] <- tibble(
    predictor_set = paste0(
      "FINAL7 + ",
      v
    ),
    valid_cells = n,
    cells_lost_vs_FINAL7 =
      n_final7 - n,
    pct_retained_vs_FINAL7 =
      100 * n / n_final7
  )
}

incremental_rows[["ALL10"]] <- tibble(
  predictor_set = "ALL10",
  valid_cells = n_all10_rebuilt,
  cells_lost_vs_FINAL7 =
    n_final7 - n_all10_rebuilt,
  pct_retained_vs_FINAL7 =
    100 * n_all10_rebuilt / n_final7
)

incremental_support <- bind_rows(
  incremental_rows
) %>%
  arrange(
    valid_cells
  )

write_csv(
  incremental_support,
  file.path(
    audit_dir,
    "04C_incremental_predictor_mask_effect.csv"
  )
)

# ------------------------------------------------------------
# 8. Load canonical step-03 reconciliation
# ------------------------------------------------------------

step03_file <- file.path(
  work_dir,
  "03_Ribotta_REV1_reconciliation_VALID.rds"
)

if (!file.exists(step03_file)) {
  stop(
    "Canonical step-03 RDS not found: ",
    step03_file
  )
}

step03 <- readRDS(
  step03_file
)

sara_points <- step03$sara_points

sara_sf <- sf::st_as_sf(
  sara_points,
  coords = c(
    "sara_lon",
    "sara_lat"
  ),
  crs = 4326,
  remove = FALSE
)

# ------------------------------------------------------------
# 9. Assign Sara points independently to ALL10 and FINAL7 masks
# ------------------------------------------------------------

assign_to_mask <- function(
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

  vals <- terra::values(
    mask_rast,
    mat = FALSE
  )

  exact_ok <-
    !is.na(cell_guess) &
    !is.na(vals[cell_guess])

  assigned <- rep(
    NA_integer_,
    nrow(points_sf)
  )

  dist_m <- rep(
    NA_real_,
    nrow(points_sf)
  )

  type <- rep(
    "outside_1500m",
    nrow(points_sf)
  )

  if (any(exact_ok)) {

    assigned[exact_ok] <-
      cell_guess[exact_ok]

    dist_m[exact_ok] <- 0

    type[exact_ok] <-
      "exact_valid_cell"
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
      coords = c("x", "y"),
      crs = 4326
    )

    valid_m <- sf::st_transform(
      valid_sf,
      3857
    )

    rem_m <- sf::st_transform(
      points_sf[
        rem,
      ],
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

    dist_m[rem] <- d

    keep <- is.finite(d) &
      d <= max_snap_m

    if (any(keep)) {

      rr <- rem[keep]

      assigned[rr] <-
        valid_df$cell[
          idx[keep]
        ]

      type[rr] <-
        "nearest_valid_cell"
    }
  }

  tibble(
    assigned_cell = assigned,
    distance_m = dist_m,
    assignment_type = type,
    retained =
      !is.na(assigned)
  )
}

sara_all10 <- assign_to_mask(
  sara_sf,
  mask_all10_rebuilt,
  1500
)

sara_final7 <- assign_to_mask(
  sara_sf,
  mask_final7,
  1500
)

sara_mask_compare <- bind_cols(
  sara_points %>%
    select(
      sara_id,
      sara_lon,
      sara_lat,
      country_code,
      country_name,
      nearest_flagged_distance_m,
      gbif_match_1km
    ),

  sara_all10 %>%
    rename(
      all10_cell = assigned_cell,
      all10_distance_m = distance_m,
      all10_assignment = assignment_type,
      all10_retained = retained
    ),

  sara_final7 %>%
    rename(
      final7_cell = assigned_cell,
      final7_distance_m = distance_m,
      final7_assignment = assignment_type,
      final7_retained = retained
    )
)

sara_mask_compare <- sara_mask_compare %>%
  mutate(
    recovered_by_FINAL7 =
      !all10_retained &
      final7_retained,

    retained_both =
      all10_retained &
      final7_retained,

    outside_both =
      !all10_retained &
      !final7_retained
  )

write_csv(
  sara_mask_compare,
  file.path(
    audit_dir,
    "04D_Sara_ALL10_vs_FINAL7_mask.csv"
  )
)

# ------------------------------------------------------------
# 9B. Direct predictor extraction at Sara/Ribotta coordinates
#
# This distinguishes:
# - climate available in the unmasked Europe bbox;
# - cell retained by the geographic coastal TIFF;
# - predictor-specific NA inside the coastal TIFF.
# ------------------------------------------------------------

sara_vect <- terra::vect(
  sara_points,
  geom = c("sara_lon", "sara_lat"),
  crs = "EPSG:4326"
)

bbox_at_sara <- terra::extract(
  bbox10,
  sara_vect,
  ID = FALSE
)

coast_at_sara <- terra::extract(
  clim10,
  sara_vect,
  ID = FALSE
)

names(bbox_at_sara) <- paste0(
  "bbox_",
  names(bbox_at_sara)
)

names(coast_at_sara) <- paste0(
  "coast_",
  names(coast_at_sara)
)

sara_predictor_support <- bind_cols(
  sara_points %>%
    select(
      sara_id,
      sara_lon,
      sara_lat,
      country_code,
      country_name,
      nearest_flagged_distance_m,
      gbif_match_1km
    ),
  bbox_at_sara,
  coast_at_sara
)

for (v in expected_vars) {

  sara_predictor_support[[
    paste0("bbox_valid_", v)
  ]] <- !is.na(
    sara_predictor_support[[
      paste0("bbox_", v)
    ]]
  )

  sara_predictor_support[[
    paste0("coast_valid_", v)
  ]] <- !is.na(
    sara_predictor_support[[
      paste0("coast_", v)
    ]]
  )
}

sara_predictor_support$n_bbox_valid_10 <- rowSums(
  as.data.frame(
    sara_predictor_support[
      paste0(
        "bbox_valid_",
        expected_vars
      )
    ]
  )
)

sara_predictor_support$n_coast_valid_10 <- rowSums(
  as.data.frame(
    sara_predictor_support[
      paste0(
        "coast_valid_",
        expected_vars
      )
    ]
  )
)

sara_predictor_support$n_bbox_valid_FINAL7 <- rowSums(
  as.data.frame(
    sara_predictor_support[
      paste0(
        "bbox_valid_",
        final7_vars
      )
    ]
  )
)

sara_predictor_support$n_coast_valid_FINAL7 <- rowSums(
  as.data.frame(
    sara_predictor_support[
      paste0(
        "coast_valid_",
        final7_vars
      )
    ]
  )
)

sara_predictor_support$bbox_FINAL7_complete <-
  sara_predictor_support$n_bbox_valid_FINAL7 ==
  length(final7_vars)

sara_predictor_support$coast_FINAL7_complete <-
  sara_predictor_support$n_coast_valid_FINAL7 ==
  length(final7_vars)

sara_predictor_support$bbox_ALL10_complete <-
  sara_predictor_support$n_bbox_valid_10 ==
  length(expected_vars)

sara_predictor_support$coast_ALL10_complete <-
  sara_predictor_support$n_coast_valid_10 ==
  length(expected_vars)

write_csv(
  sara_predictor_support,
  file.path(
    audit_dir,
    "04D2_Sara_predictor_support_bbox_vs_coast.csv"
  )
)

critical_predictor_support <- sara_predictor_support %>%
  filter(
    country_code %in%
      c("IT", "FR", "GB")
  ) %>%
  arrange(
    country_code,
    sara_id
  )

write_csv(
  critical_predictor_support,
  file.path(
    audit_dir,
    "04D3_Italy_France_UK_predictor_support.csv"
  )
)

critical_support_summary <- critical_predictor_support %>%
  group_by(
    country_code,
    country_name
  ) %>%
  summarise(
    n_points = n(),

    n_bbox_FINAL7_complete =
      sum(
        bbox_FINAL7_complete,
        na.rm = TRUE
      ),

    n_coast_FINAL7_complete =
      sum(
        coast_FINAL7_complete,
        na.rm = TRUE
      ),

    n_bbox_ALL10_complete =
      sum(
        bbox_ALL10_complete,
        na.rm = TRUE
      ),

    n_coast_ALL10_complete =
      sum(
        coast_ALL10_complete,
        na.rm = TRUE
      ),

    n_bbox_valid_fcf =
      sum(
        bbox_valid_fcf,
        na.rm = TRUE
      ),

    n_coast_valid_fcf =
      sum(
        coast_valid_fcf,
        na.rm = TRUE
      ),

    .groups = "drop"
  )

write_csv(
  critical_support_summary,
  file.path(
    audit_dir,
    "04D4_critical_country_predictor_support_summary.csv"
  )
)

# ------------------------------------------------------------
# 10. Reviewer-country comparison
# ------------------------------------------------------------

reviewer_codes <- c(
  "GB", "FR", "IT",
  "ES", "PT",
  "SI", "HR", "BA", "ME", "AL", "GR",
  "TR", "BG", "RO",
  "CY", "MT"
)

reviewer_mask_compare <- sara_mask_compare %>%
  filter(
    country_code %in%
      reviewer_codes
  ) %>%
  group_by(
    country_code,
    country_name
  ) %>%
  summarise(
    n_points = n(),

    n_retained_ALL10 =
      sum(
        all10_retained,
        na.rm = TRUE
      ),

    n_retained_FINAL7 =
      sum(
        final7_retained,
        na.rm = TRUE
      ),

    n_recovered_by_FINAL7 =
      sum(
        recovered_by_FINAL7,
        na.rm = TRUE
      ),

    n_GBIF_within_1km =
      sum(
        gbif_match_1km,
        na.rm = TRUE
      ),

    min_ALL10_distance_m =
      min(
        all10_distance_m,
        na.rm = TRUE
      ),

    min_FINAL7_distance_m =
      min(
        final7_distance_m,
        na.rm = TRUE
      ),

    median_ALL10_distance_m =
      median(
        all10_distance_m,
        na.rm = TRUE
      ),

    median_FINAL7_distance_m =
      median(
        final7_distance_m,
        na.rm = TRUE
      ),

    .groups = "drop"
  )

for (nm in c(
  "min_ALL10_distance_m",
  "min_FINAL7_distance_m",
  "median_ALL10_distance_m",
  "median_FINAL7_distance_m"
)) {

  reviewer_mask_compare[[nm]][
    !is.finite(
      reviewer_mask_compare[[nm]]
    )
  ] <- NA_real_
}

write_csv(
  reviewer_mask_compare,
  file.path(
    audit_dir,
    "04E_reviewer_country_ALL10_vs_FINAL7.csv"
  )
)

# ------------------------------------------------------------
# 11. Geographic belt diagnostic, if saved belt is available
# ------------------------------------------------------------

belt_summary <- tibble(
  belt_available = !is.na(belt_file),
  belt_area_km2 = NA_real_,
  belt_source = ifelse(
    is.na(belt_file),
    "saved belt RDS unavailable; using coastal TIFF support diagnostics",
    "saved geographic belt RDS"
  )
)

if (!is.na(belt_file)) {

  belt <- readRDS(
    belt_file
  )

  if (
    inherits(belt, "sf") ||
    inherits(belt, "sfc")
  ) {

    belt_area <- as.numeric(
      sum(
        sf::st_area(
          sf::st_transform(
            belt,
            6933
          )
        )
      )
    ) / 1e6

    belt_summary <- tibble(
      belt_available = TRUE,
      belt_area_km2 = belt_area,
      belt_source = "saved geographic belt RDS"
    )
  }
}

write_csv(
  belt_summary,
  file.path(
    audit_dir,
    "04F_geographic_belt_summary.csv"
  )
)

# ------------------------------------------------------------
# 12. Key synthesis
# ------------------------------------------------------------

summary_tbl <- tibble(
  metric = c(
    "REV1 saved ALL10 valid cells (if RDS available)",
    "Rebuilt ALL10 valid cells from coastal TIFF",
    "Saved vs rebuilt ALL10 count matches",
    "Rebuilt ALL10 matches REV1 reported 79,579",
    "FINAL7 current valid cells",
    "Cells recovered by FINAL7 vs ALL10",
    "FINAL7 percent of ALL10",
    "Sara points retained under ALL10",
    "Sara points retained under FINAL7",
    "Sara points recovered under FINAL7",
    "Sara unique cells under ALL10",
    "Sara unique cells under FINAL7"
  ),

  value = c(
    as.character(n_mask10_saved),
    as.character(n_all10_rebuilt),
    as.character(all10_rebuild_matches_saved),
    as.character(all10_matches_reported_REV1),
    as.character(n_final7),
    as.character(
      n_final7 -
        n_all10_rebuilt
    ),
    as.character(
      round(
        100 * n_final7 /
          n_all10_rebuilt,
        3
      )
    ),
    as.character(
      sum(
        sara_mask_compare$all10_retained,
        na.rm = TRUE
      )
    ),
    as.character(
      sum(
        sara_mask_compare$final7_retained,
        na.rm = TRUE
      )
    ),
    as.character(
      sum(
        sara_mask_compare$recovered_by_FINAL7,
        na.rm = TRUE
      )
    ),
    as.character(
      dplyr::n_distinct(
        sara_mask_compare$all10_cell[
          sara_mask_compare$all10_retained
        ]
      )
    ),
    as.character(
      dplyr::n_distinct(
        sara_mask_compare$final7_cell[
          sara_mask_compare$final7_retained
        ]
      )
    )
  )
)

write_csv(
  summary_tbl,
  file.path(
    audit_dir,
    "04G_mask_support_summary.csv"
  )
)

saveRDS(
  list(
    asset_inventory = asset_inventory,
    variable_support = variable_support,
    bbox_coast_support = bbox_coast_support,
    incremental_support = incremental_support,
    sara_predictor_support = sara_predictor_support,
    critical_support_summary = critical_support_summary,
    sara_mask_compare = sara_mask_compare,
    reviewer_mask_compare = reviewer_mask_compare,
    belt_summary = belt_summary,
    summary = summary_tbl,
    mask10_file = mask10_file,
    coast_stack_file = coast_stack_file
  ),
  file.path(
    work_dir,
    "04_REV1_coastal_mask_support_VALID.rds"
  )
)

# ------------------------------------------------------------
# 13. Update canonical workflow manifest
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
  "",
  "SUPERSEDED / INVALIDATED",
  "03_reconcile_Ribotta_with_REV1_CORRECTED.R",
  "04_trace_missing_Ribotta_cells.R",
  "05_audit_Ribotta_geography_by_country.R",
  "",
  "IMPORTANT",
  "The canonical step 03 reproduces the actual REV1 all-10-variable raster-valid mask.",
  "The canonical step 04 reconstructs the ALL10 mask directly from the official coastal TIFF when the saved mask RDS is unavailable, compares it with FINAL7 support, and checks BBOX-vs-COAST predictor values at Sara/Ribotta points.",
  "No model is refitted in steps 00-04."
)

writeLines(
  canonical_workflow,
  manifest_file
)

# ------------------------------------------------------------
# 14. Console report
# ------------------------------------------------------------

cat(
  "\n============================================================\n"
)

cat(
  "REV2 - REV1 COASTAL MASK SUPPORT AUDIT COMPLETE\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "Global mask summary:\n"
)

print(
  as.data.frame(
    summary_tbl
  ),
  row.names = FALSE
)

cat(
  "\nPer-variable current support:\n"
)

print(
  as.data.frame(
    variable_support
  ),
  row.names = FALSE
)

cat(
  "\nIncremental predictor effect relative to FINAL7:\n"
)

print(
  as.data.frame(
    incremental_support
  ),
  row.names = FALSE
)

cat(
  "\nReviewer-country ALL10 vs FINAL7:\n"
)

print(
  as.data.frame(
    reviewer_mask_compare
  ),
  row.names = FALSE
)


cat(
  "\nItaly / France / UK direct predictor support (BBOX vs COAST):\n"
)

print(
  as.data.frame(
    critical_support_summary
  ),
  row.names = FALSE
)

cat(
  "\nNo occurrence has been added. No model has been fitted.\n"
)

cat(
  "============================================================\n"
)
