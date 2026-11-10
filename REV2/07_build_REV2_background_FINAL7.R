# ============================================================
# ROSA RUGOSA - REVISION 2
# 07_build_REV2_background_FINAL7.R
#
# Purpose
# -------
# Build and freeze the European coastal background for REV2.
#
# The background follows the REV1 design:
#   - FINAL7 raster-valid European coastal domain;
#   - exclusion of cells within 1,000 m of any frozen presence;
#   - 10,000 background cells sampled at random.
#
# REV2 adds an explicit random seed and samples unique raster cells
# without replacement for full reproducibility.
#
# No occurrence is modified.
# No model is fitted.
# ============================================================

rm(list = ls())
gc()

required_pkgs <- c(
  "terra", "dplyr", "readr", "tibble"
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
library(dplyr)
library(readr)
library(tibble)

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

old_07 <- c(
  list.files(
    audit_dir,
    pattern = "^07[A-Z]_",
    full.names = TRUE
  ),
  file.path(
    work_dir,
    "background_REV2_FINAL7.csv"
  ),
  file.path(
    work_dir,
    "07_REV2_background_FINAL7_FROZEN.rds"
  )
)

old_07 <- old_07[
  file.exists(old_07)
]

if (length(old_07) > 0) {
  unlink(
    old_07,
    force = TRUE
  )
}

# ------------------------------------------------------------
# 2. Canonical inputs
# ------------------------------------------------------------

occ_file <- file.path(
  work_dir,
  "occurrences_REV2_FROZEN.csv"
)

step06_file <- file.path(
  work_dir,
  "06_REV2_occurrence_dataset_FROZEN.rds"
)

if (!file.exists(occ_file)) {
  stop(
    "Missing frozen occurrence file: ",
    occ_file,
    "\nRun 06_freeze_REV2_occurrence_dataset.R first."
  )
}

if (!file.exists(step06_file)) {
  stop(
    "Missing canonical step-06 RDS: ",
    step06_file
  )
}

occ <- readr::read_csv(
  occ_file,
  show_col_types = FALSE
)

step06 <- readRDS(
  step06_file
)

required_occ_cols <- c(
  "occurrence_id",
  "assigned_longitude",
  "assigned_latitude",
  "cell_id",
  "calibration_presence"
)

missing_occ_cols <- setdiff(
  required_occ_cols,
  names(occ)
)

if (length(missing_occ_cols) > 0) {
  stop(
    "Frozen occurrence file lacks: ",
    paste(
      missing_occ_cols,
      collapse = ", "
    )
  )
}

if (nrow(occ) != 8709L) {
  stop(
    "Expected 8,709 frozen presences, found ",
    nrow(occ),
    "."
  )
}

if (dplyr::n_distinct(
  occ$cell_id
) != 8709L) {
  stop(
    "Frozen occurrence cells are not unique."
  )
}

if (!all(
  occ$calibration_presence
)) {
  stop(
    "Frozen occurrence file contains non-calibration records."
  )
}

final7_vars <- step06$final7_vars

if (!identical(
  final7_vars,
  c(
    "bio3",
    "bio8",
    "bio9",
    "bio15",
    "bio18",
    "gsl",
    "gst"
  )
)) {
  stop(
    "Unexpected FINAL7 predictor definition in step 06."
  )
}

# ------------------------------------------------------------
# 3. Locate official current coastal raster
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

clim <- terra::rast(
  coast_file
)

missing_final7 <- setdiff(
  final7_vars,
  names(clim)
)

if (length(missing_final7) > 0) {
  stop(
    "Current coastal raster lacks FINAL7 variables: ",
    paste(
      missing_final7,
      collapse = ", "
    )
  )
}

clim7 <- clim[
  [final7_vars]
]

# ------------------------------------------------------------
# 4. Rebuild FINAL7 valid-cell mask
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

if (n_final7_cells != 210748L) {
  stop(
    "FINAL7 mask should contain 210,748 cells, found ",
    n_final7_cells,
    "."
  )
}

# ------------------------------------------------------------
# 5. Build 1-km presence-exclusion mask
# ------------------------------------------------------------

presence_cells <- as.integer(
  occ$cell_id
)

if (any(
  is.na(
    presence_cells
  )
)) {
  stop(
    "Frozen occurrence set contains NA cell_id."
  )
}

mask_vals <- terra::values(
  final7_mask,
  mat = FALSE
)

if (any(
  is.na(
    mask_vals[
      presence_cells
    ]
  )
)) {
  stop(
    "At least one frozen occurrence cell is outside the FINAL7 mask."
  )
}

presence_raster <- final7_mask
terra::values(
  presence_raster
) <- NA_real_

presence_raster[
  presence_cells
] <- 1

distance_to_presence <- terra::distance(
  presence_raster
)

background_exclusion_m <- 1000

eligible_mask <- terra::ifel(
  final7_mask == 1 &
    distance_to_presence >
      background_exclusion_m,
  1,
  NA
)

eligible_df <- terra::as.data.frame(
  eligible_mask,
  xy = TRUE,
  cells = TRUE,
  na.rm = TRUE
)

eligible_df <- eligible_df[
  ,
  c(
    "cell",
    "x",
    "y"
  )
]

names(eligible_df) <- c(
  "cell_id",
  "decimal_longitude",
  "decimal_latitude"
)

eligible_df$cell_id <- as.integer(
  eligible_df$cell_id
)

n_eligible <- nrow(
  eligible_df
)

if (n_eligible < 10000L) {
  stop(
    "Only ",
    n_eligible,
    " background-eligible cells remain after 1-km exclusion."
  )
}

if (any(
  eligible_df$cell_id %in%
    presence_cells
)) {
  stop(
    "Presence cells remain in the background-eligible pool."
  )
}

# ------------------------------------------------------------
# 6. Deterministic sample of 10,000 unique background cells
# ------------------------------------------------------------

background_n <- 10000L
background_seed <- 12345L

set.seed(
  background_seed
)

sample_idx <- sample(
  seq_len(
    n_eligible
  ),
  size = background_n,
  replace = FALSE
)

background <- eligible_df[
  sample_idx,
  ,
  drop = FALSE
] %>%
  arrange(
    cell_id
  )

if (nrow(background) != background_n) {
  stop(
    "Background sample does not contain exactly 10,000 rows."
  )
}

if (dplyr::n_distinct(
  background$cell_id
) != background_n) {
  stop(
    "Background sample contains duplicated raster cells."
  )
}

if (any(
  background$cell_id %in%
    presence_cells
)) {
  stop(
    "Background sample overlaps frozen presence cells."
  )
}

# ------------------------------------------------------------
# 7. Extract FINAL7 predictors and nearest-presence distance
# ------------------------------------------------------------

bg_xy <- as.matrix(
  background[
    ,
    c(
      "decimal_longitude",
      "decimal_latitude"
    )
  ]
)

bg_pred <- terra::extract(
  clim7,
  bg_xy,
  ID = FALSE
)

background <- bind_cols(
  background,
  as.data.frame(
    bg_pred
  )
)

bg_dist <- terra::extract(
  distance_to_presence,
  bg_xy,
  ID = FALSE
)

if (ncol(bg_dist) != 1) {
  stop(
    "Unexpected distance-raster extraction."
  )
}

background$distance_to_nearest_presence_m <-
  as.numeric(
    bg_dist[
      ,
      1
    ]
  )

background$pa <- 0L
background$source <- "background"
background$background_seed <- background_seed

if (any(
  !is.finite(
    background$distance_to_nearest_presence_m
  )
)) {
  stop(
    "Background sample contains non-finite distances to presences."
  )
}

if (any(
  background$distance_to_nearest_presence_m <=
    background_exclusion_m
)) {
  stop(
    "At least one background cell violates the 1-km exclusion rule."
  )
}

if (any(
  !stats::complete.cases(
    background[
      ,
      final7_vars,
      drop = FALSE
    ]
  )
)) {
  stop(
    "Background sample contains missing FINAL7 predictor values."
  )
}

# ------------------------------------------------------------
# 8. Extract FINAL7 predictors for frozen presences
# ------------------------------------------------------------

presence_xy <- as.matrix(
  occ[
    ,
    c(
      "assigned_longitude",
      "assigned_latitude"
    )
  ]
)

presence_pred <- terra::extract(
  clim7,
  presence_xy,
  ID = FALSE
)

presence_model <- bind_cols(
  occ,
  as.data.frame(
    presence_pred
  )
)

presence_model$pa <- 1L

if (any(
  !stats::complete.cases(
    presence_model[
      ,
      final7_vars,
      drop = FALSE
    ]
  )
)) {
  stop(
    "Frozen presence set contains missing FINAL7 predictor values."
  )
}

# ------------------------------------------------------------
# 9. Combined model-input table
# ------------------------------------------------------------

presence_model_min <- presence_model %>%
  transmute(
    pa = 1L,
    cell_id =
      as.integer(
        cell_id
      ),
    decimal_longitude =
      as.numeric(
        assigned_longitude
      ),
    decimal_latitude =
      as.numeric(
        assigned_latitude
      ),
    source =
      as.character(
        source
      ),
    occurrence_id =
      as.character(
        occurrence_id
      ),
    across(
      all_of(
        final7_vars
      )
    )
  )

background_model_min <- background %>%
  transmute(
    pa = 0L,
    cell_id =
      as.integer(
        cell_id
      ),
    decimal_longitude,
    decimal_latitude,
    source,
    occurrence_id =
      NA_character_,
    across(
      all_of(
        final7_vars
      )
    )
  )

model_input <- bind_rows(
  presence_model_min,
  background_model_min
)

if (nrow(model_input) != 18709L) {
  stop(
    "Expected 18,709 rows in combined model input, found ",
    nrow(model_input),
    "."
  )
}

if (sum(
  model_input$pa == 1
) != 8709L) {
  stop(
    "Combined model input does not contain 8,709 presences."
  )
}

if (sum(
  model_input$pa == 0
) != 10000L) {
  stop(
    "Combined model input does not contain 10,000 background cells."
  )
}

# ------------------------------------------------------------
# 10. Audit summaries
# ------------------------------------------------------------

summary_tbl <- tibble(
  metric = c(
    "FINAL7 valid coastal cells",
    "Frozen presence cells",
    "Background exclusion distance m",
    "Eligible background cells after exclusion",
    "Background seed",
    "Sampled background cells",
    "Unique sampled background cells",
    "Minimum sampled distance to presence m",
    "Median sampled distance to presence m",
    "Maximum sampled distance to presence m",
    "Combined model-input rows",
    "Combined presence rows",
    "Combined background rows"
  ),

  value = c(
    n_final7_cells,
    nrow(
      occ
    ),
    background_exclusion_m,
    n_eligible,
    background_seed,
    nrow(
      background
    ),
    dplyr::n_distinct(
      background$cell_id
    ),
    min(
      background$distance_to_nearest_presence_m
    ),
    stats::median(
      background$distance_to_nearest_presence_m
    ),
    max(
      background$distance_to_nearest_presence_m
    ),
    nrow(
      model_input
    ),
    sum(
      model_input$pa == 1
    ),
    sum(
      model_input$pa == 0
    )
  )
)

# ------------------------------------------------------------
# 11. Save canonical outputs
# ------------------------------------------------------------

readr::write_csv(
  background,
  file.path(
    work_dir,
    "background_REV2_FINAL7.csv"
  )
)

readr::write_csv(
  model_input,
  file.path(
    work_dir,
    "model_input_REV2_FINAL7.csv"
  )
)

readr::write_csv(
  summary_tbl,
  file.path(
    audit_dir,
    "07A_background_summary.csv"
  )
)

predictor_summary <- bind_rows(
  lapply(
    final7_vars,
    function(v) {

      tibble(
        predictor = v,
        group = c(
          "presence",
          "background"
        ),
        n = c(
          nrow(
            presence_model
          ),
          nrow(
            background
          )
        ),
        mean = c(
          mean(
            presence_model[[v]],
            na.rm = TRUE
          ),
          mean(
            background[[v]],
            na.rm = TRUE
          )
        ),
        sd = c(
          stats::sd(
            presence_model[[v]],
            na.rm = TRUE
          ),
          stats::sd(
            background[[v]],
            na.rm = TRUE
          )
        ),
        min = c(
          min(
            presence_model[[v]],
            na.rm = TRUE
          ),
          min(
            background[[v]],
            na.rm = TRUE
          )
        ),
        median = c(
          stats::median(
            presence_model[[v]],
            na.rm = TRUE
          ),
          stats::median(
            background[[v]],
            na.rm = TRUE
          )
        ),
        max = c(
          max(
            presence_model[[v]],
            na.rm = TRUE
          ),
          max(
            background[[v]],
            na.rm = TRUE
          )
        )
      )
    }
  )
)

readr::write_csv(
  predictor_summary,
  file.path(
    audit_dir,
    "07B_presence_background_predictor_summary.csv"
  )
)

saveRDS(
  list(
    background = background,
    presence_model = presence_model,
    model_input = model_input,
    summary = summary_tbl,
    predictor_summary = predictor_summary,
    final7_vars = final7_vars,
    background_n = background_n,
    background_seed = background_seed,
    background_exclusion_m = background_exclusion_m,
    n_eligible_background_cells = n_eligible,
    coast_file = coast_file
  ),
  file.path(
    work_dir,
    "07_REV2_background_FINAL7_FROZEN.rds"
  )
)

# ------------------------------------------------------------
# 12. Update canonical workflow manifest
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
  "SUPERSEDED / INVALIDATED",
  "03_reconcile_Ribotta_with_REV1_CORRECTED.R",
  "04_trace_missing_Ribotta_cells.R",
  "05_audit_Ribotta_geography_by_country.R",
  "",
  "IMPORTANT",
  "Step 06 freezes the invaded calibration set at 8,709 unique cells: 8,705 GBIF + 4 independent Brussa field cells.",
  "Step 07 samples 10,000 unique FINAL7 background cells after excluding locations within 1,000 m of any frozen presence.",
  "Background sampling is reproducible with random seed 12345.",
  "Background records are reference/background points; their role will be handled separately for presence-background and binary classification algorithms.",
  "No model is fitted in steps 00-07."
)

writeLines(
  canonical_workflow,
  manifest_file
)

# ------------------------------------------------------------
# 13. Console report
# ------------------------------------------------------------

cat(
  "\n============================================================\n"
)

cat(
  "REV2 - FINAL7 BACKGROUND FROZEN\n"
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
  "\nBackground file:\n",
  file.path(
    work_dir,
    "background_REV2_FINAL7.csv"
  ),
  "\n"
)

cat(
  "\nCombined model-input file:\n",
  file.path(
    work_dir,
    "model_input_REV2_FINAL7.csv"
  ),
  "\n"
)

cat(
  "\nNo occurrence has been modified. No model has been fitted.\n"
)

cat(
  "============================================================\n"
)
