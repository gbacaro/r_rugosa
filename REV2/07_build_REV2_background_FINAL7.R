# ============================================================
# ROSA RUGOSA - REVISION 2
# 07_build_REV2_background_FINAL7.R
#
# Purpose
# -------
# Build and freeze the European background from the corrected
# true-1km FINAL7 domain produced by steps 05-06.
#
# Background design:
#   - corrected FINAL7 raster-valid coastal domain;
#   - exclusion of cells within 1,000 m of any frozen presence;
#   - 10,000 unique background cells sampled without replacement;
#   - explicit random seed for reproducibility.
#
# No occurrence is modified and no model is fitted.
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
    "model_input_REV2_FINAL7.csv"
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
# 2. Canonical corrected inputs
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
    "\nRun corrected 06_freeze_REV2_occurrence_dataset.R first."
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

required_step06 <- c(
  "final7_vars",
  "n_final7_cells",
  "frozen_n",
  "bbox_file",
  "mask_file",
  "distance_crs",
  "coastal_buffer_m",
  "raster_inclusion_rule"
)

missing_step06 <- required_step06[
  !vapply(
    required_step06,
    function(x) {
      !is.null(
        step06[[x]]
      )
    },
    logical(1)
  )
]

if (length(missing_step06) > 0) {
  stop(
    "Step-06 RDS lacks: ",
    paste(
      missing_step06,
      collapse = ", "
    )
  )
}

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

frozen_n <- as.integer(
  step06$frozen_n
)

if (nrow(occ) != frozen_n) {
  stop(
    "Frozen occurrence CSV contains ",
    nrow(occ),
    " rows but step 06 records ",
    frozen_n,
    "."
  )
}

if (dplyr::n_distinct(
  occ$cell_id
) != frozen_n) {
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

expected_final7_vars <- c(
  "bio3",
  "bio8",
  "bio9",
  "bio15",
  "bio18",
  "gsl",
  "gst"
)

if (!identical(
  final7_vars,
  expected_final7_vars
)) {
  stop(
    "Unexpected FINAL7 predictor definition in step 06."
  )
}

bbox_file <- step06$bbox_file
mask_file <- step06$mask_file

if (!file.exists(bbox_file)) {
  stop(
    "Corrected European bbox raster is missing: ",
    bbox_file
  )
}

if (!file.exists(mask_file)) {
  stop(
    "Corrected FINAL7 mask is missing: ",
    mask_file
  )
}

# ------------------------------------------------------------
# 3. Climate stack and corrected FINAL7 mask
# ------------------------------------------------------------

clim_bbox <- terra::rast(
  bbox_file
)

missing_final7 <- setdiff(
  final7_vars,
  names(clim_bbox)
)

if (length(missing_final7) > 0) {
  stop(
    "European bbox raster lacks FINAL7 variables: ",
    paste(
      missing_final7,
      collapse = ", "
    )
  )
}

clim7 <- terra::subset(
  clim_bbox,
  final7_vars
)

if (terra::nlyr(clim7) !=
    length(final7_vars)) {
  stop(
    "FINAL7 layer selection returned ",
    terra::nlyr(clim7),
    " layers instead of ",
    length(final7_vars),
    "."
  )
}

names(clim7) <- final7_vars

final7_mask <- terra::rast(
  mask_file
)

n_final7_cells <- as.integer(
  terra::global(
    final7_mask,
    "sum",
    na.rm = TRUE
  )[1, 1]
)

if (n_final7_cells !=
    as.integer(
      step06$n_final7_cells
    )) {
  stop(
    "Corrected FINAL7 mask size differs from step 06."
  )
}

# ------------------------------------------------------------
# 4. Validate presence cells on corrected mask
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

if (any(
  presence_cells < 1L |
    presence_cells >
      terra::ncell(
        final7_mask
      )
)) {
  stop(
    "Frozen occurrence set contains invalid corrected-mask cell IDs."
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
    "At least one frozen occurrence cell is outside the corrected FINAL7 mask."
  )
}

# ------------------------------------------------------------
# 5. Exclude cells within 1 km of frozen presences
# ------------------------------------------------------------

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

if (nrow(background) !=
    background_n) {
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
# 7. Extract FINAL7 predictors and distance diagnostics
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

bg_pred <- as.data.frame(
  terra::extract(
    clim7,
    bg_xy
  )
)

if (ncol(bg_pred) !=
    length(final7_vars)) {
  stop(
    "Background predictor extraction returned ",
    ncol(bg_pred),
    " columns instead of ",
    length(final7_vars),
    "."
  )
}

names(bg_pred) <- final7_vars

background <- bind_cols(
  background,
  bg_pred
)

bg_dist <- terra::extract(
  distance_to_presence,
  bg_xy
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
background$background_seed <-
  background_seed

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

presence_pred <- as.data.frame(
  terra::extract(
    clim7,
    presence_xy
  )
)

if (ncol(presence_pred) !=
    length(final7_vars)) {
  stop(
    "Presence predictor extraction returned ",
    ncol(presence_pred),
    " columns instead of ",
    length(final7_vars),
    "."
  )
}

names(presence_pred) <- final7_vars

presence_model <- bind_cols(
  occ,
  presence_pred
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

expected_model_rows <- frozen_n +
  background_n

if (nrow(model_input) !=
    expected_model_rows) {
  stop(
    "Expected ",
    expected_model_rows,
    " rows in combined model input, found ",
    nrow(model_input),
    "."
  )
}

if (sum(
  model_input$pa == 1
) != frozen_n) {
  stop(
    "Combined model input presence count is inconsistent with step 06."
  )
}

if (sum(
  model_input$pa == 0
) != background_n) {
  stop(
    "Combined model input does not contain 10,000 background cells."
  )
}

# ------------------------------------------------------------
# 10. Audit summaries
# ------------------------------------------------------------

summary_tbl <- tibble(
  metric = c(
    "Corrected FINAL7 valid coastal cells",
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
    frozen_n,
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
    n_final7_cells = n_final7_cells,
    frozen_n = frozen_n,
    bbox_file = bbox_file,
    mask_file = mask_file,
    distance_crs = step06$distance_crs,
    coastal_buffer_m = step06$coastal_buffer_m,
    raster_inclusion_rule =
      step06$raster_inclusion_rule
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
  "IMPORTANT",
  "Step 05 rebuilds the corrected 1,000-m terrestrial European coastal domain from the unmasked current-climate raster in EPSG:3035.",
  "Step 06 freezes all corrected GBIF cells plus four independently verified Brussa field cells.",
  "Step 07 samples 10,000 unique corrected-domain background cells after a 1,000-m exclusion from frozen presences.",
  "Background sampling is reproducible with random seed 12345.",
  "Background records are reference/background points; their algorithm-specific role is handled at model fitting.",
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
  "REV2 - CORRECTED FINAL7 BACKGROUND FROZEN\n"
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
