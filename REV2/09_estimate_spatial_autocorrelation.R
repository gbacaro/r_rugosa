# ============================================================
# ROSA RUGOSA - REVISION 2
# 09_estimate_spatial_autocorrelation.R
#
# Purpose
# -------
# Estimate the spatial autocorrelation range of the seven
# environmental predictors over the corrected REV2 coastal
# domain and use it as an exploratory guide for spatial
# cross-validation block size.
#
# The analysis is performed on a reproducible sample of raster
# cells transformed to EPSG:3035 so that variogram distances are
# measured in metres.
#
# This step does not fit an SDM, does not select a weighting
# strategy, and does not modify occurrences or background data.
# ============================================================

rm(list = ls())
gc()

required_pkgs <- c(
  "terra", "sf", "dplyr", "readr", "tibble", "blockCV"
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
    paste(
      missing_pkgs,
      collapse = ", "
    ),
    "\nInstall missing packages before running this script. ",
    "For blockCV use: install.packages(\"blockCV\")"
  )
}

library(terra)
library(sf)
library(dplyr)
library(readr)
library(tibble)
library(blockCV)

sf::sf_use_s2(TRUE)

if (utils::packageVersion(
  "blockCV"
) < "3.0.0") {
  stop(
    "blockCV >= 3.0.0 is required. Installed version: ",
    as.character(
      utils::packageVersion(
        "blockCV"
      )
    )
  )
}

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

old_09 <- c(
  list.files(
    audit_dir,
    pattern = "^09[A-Z]_",
    full.names = TRUE
  ),
  file.path(
    work_dir,
    "09_REV2_spatial_autocorrelation.rds"
  )
)

old_09 <- old_09[
  file.exists(
    old_09
  )
]

if (length(old_09) > 0) {
  unlink(
    old_09,
    force = TRUE
  )
}

# ------------------------------------------------------------
# 2. Canonical inputs
# ------------------------------------------------------------

step06_file <- file.path(
  work_dir,
  "06_REV2_occurrence_dataset_FROZEN.rds"
)

step08_file <- file.path(
  work_dir,
  "08_REV2_sampling_bias_audit.rds"
)

if (!file.exists(step06_file)) {
  stop(
    "Missing corrected step-06 RDS: ",
    step06_file
  )
}

if (!file.exists(step08_file)) {
  stop(
    "Missing step-08 sampling-bias audit RDS: ",
    step08_file
  )
}

step06 <- readRDS(
  step06_file
)

step08 <- readRDS(
  step08_file
)

final7_vars <- step06$final7_vars
distance_crs <- as.integer(
  step06$distance_crs
)

expected_final7 <- c(
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
  expected_final7
)) {
  stop(
    "Unexpected FINAL7 predictor definition."
  )
}

if (!is.finite(
  distance_crs
)) {
  stop(
    "Invalid projected CRS inherited from step 06."
  )
}

bbox_file <- step06$bbox_file
mask_file <- step06$mask_file

if (!file.exists(bbox_file)) {
  stop(
    "European bbox climate raster is missing: ",
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
# 3. Rebuild corrected FINAL7 environmental raster
# ------------------------------------------------------------

clim_bbox <- terra::rast(
  bbox_file
)

missing_final7 <- setdiff(
  final7_vars,
  names(
    clim_bbox
  )
)

if (length(missing_final7) > 0) {
  stop(
    "European bbox raster lacks FINAL7 predictors: ",
    paste(
      missing_final7,
      collapse = ", "
    )
  )
}

clim7_bbox <- terra::subset(
  clim_bbox,
  final7_vars
)

names(
  clim7_bbox
) <- final7_vars

final7_mask <- terra::rast(
  mask_file
)

clim7_domain <- terra::crop(
  clim7_bbox,
  final7_mask
)

if (!terra::compareGeom(
  clim7_domain,
  final7_mask,
  stopOnError = FALSE
)) {
  stop(
    "Corrected FINAL7 mask and climate raster are not geometrically aligned."
  )
}

clim7_domain <- terra::mask(
  clim7_domain,
  final7_mask
)

n_complete_cells <- as.integer(
  terra::global(
    !is.na(
      clim7_domain[[1]]
    ),
    "sum",
    na.rm = TRUE
  )[1, 1]
)

if (n_complete_cells !=
    as.integer(
      step06$n_final7_cells
    )) {
  stop(
    "Corrected climate domain contains ",
    n_complete_cells,
    " valid cells instead of ",
    step06$n_final7_cells,
    "."
  )
}

# ------------------------------------------------------------
# 4. Reproducible environmental sample
# ------------------------------------------------------------

autocor_seed <- 12345L
autocor_sample_n <- 5000L

set.seed(
  autocor_seed
)

sample_df <- terra::spatSample(
  clim7_domain,
  size = autocor_sample_n,
  method = "random",
  na.rm = TRUE,
  xy = TRUE,
  values = TRUE,
  as.df = TRUE
)

if (nrow(sample_df) !=
    autocor_sample_n) {
  stop(
    "Expected ",
    autocor_sample_n,
    " sampled environmental cells, obtained ",
    nrow(sample_df),
    "."
  )
}

required_sample_cols <- c(
  "x",
  "y",
  final7_vars
)

missing_sample_cols <- setdiff(
  required_sample_cols,
  names(
    sample_df
  )
)

if (length(missing_sample_cols) > 0) {
  stop(
    "Environmental sample lacks: ",
    paste(
      missing_sample_cols,
      collapse = ", "
    )
  )
}

if (any(
  !stats::complete.cases(
    sample_df[
      ,
      final7_vars,
      drop = FALSE
    ]
  )
)) {
  stop(
    "Environmental autocorrelation sample contains missing predictor values."
  )
}

sample_sf <- sf::st_as_sf(
  sample_df,
  coords = c(
    "x",
    "y"
  ),
  crs = sf::st_crs(
    terra::crs(
      clim7_domain,
      proj = TRUE
    )
  ),
  remove = FALSE
)

sample_m <- sf::st_transform(
  sample_sf,
  distance_crs
)

# ------------------------------------------------------------
# 5. Predictor spatial autocorrelation
# ------------------------------------------------------------

sac_predictors <- blockCV::cv_spatial_autocor(
  x = sample_m,
  column = final7_vars,
  plot = FALSE,
  progress = FALSE
)

predictor_range_m <- as.numeric(
  sac_predictors$range
)

if (length(predictor_range_m) != 1 ||
    !is.finite(
      predictor_range_m
    ) ||
    predictor_range_m <= 0) {
  stop(
    "blockCV did not return a valid median predictor autocorrelation range."
  )
}

range_table <- as.data.frame(
  sac_predictors$range_table
)

if (nrow(range_table) < 1) {
  stop(
    "blockCV returned an empty predictor range table."
  )
}

# ------------------------------------------------------------
# 6. Compare autocorrelation range with step-08 candidate scales
# ------------------------------------------------------------

candidate_scales_km <- as.numeric(
  step08$block_scales_km
)

if (length(candidate_scales_km) < 1 ||
    any(
      !is.finite(
        candidate_scales_km
      )
    )) {
  stop(
    "Step 08 does not contain valid candidate block scales."
  )
}

candidate_scales_m <-
  candidate_scales_km *
  1000

eligible_candidate <- candidate_scales_m[
  candidate_scales_m >=
    predictor_range_m
]

if (length(eligible_candidate) > 0) {
  exploratory_candidate_m <-
    min(
      eligible_candidate
    )
  candidate_rule <-
    "smallest_step08_scale_not_below_median_predictor_range"
} else {
  exploratory_candidate_m <-
    max(
      candidate_scales_m
    )
  candidate_rule <-
    "median_predictor_range_exceeds_all_step08_candidate_scales"
}

exploratory_candidate_km <-
  exploratory_candidate_m /
  1000

scale_comparison <- tibble(
  candidate_block_km =
    candidate_scales_km,
  candidate_block_m =
    candidate_scales_m,
  median_predictor_autocorrelation_range_m =
    predictor_range_m,
  block_minus_range_m =
    candidate_scales_m -
    predictor_range_m,
  block_not_below_range =
    candidate_scales_m >=
    predictor_range_m,
  exploratory_candidate =
    candidate_scales_m ==
    exploratory_candidate_m
)

selected_step08_scale <- step08$spatial_scale_summary %>%
  filter(
    block_km ==
      exploratory_candidate_km
  )

selected_weight_summary <- step08$candidate_weight_summary %>%
  filter(
    block_km ==
      exploratory_candidate_km
  )

selected_cv_quality <- step08$candidate_cv_scale_quality %>%
  filter(
    block_km ==
      exploratory_candidate_km
  )

decision_support <- tibble(
  metric = c(
    "blockCV version",
    "Environmental sample cells",
    "Autocorrelation sample seed",
    "Distance CRS EPSG",
    "Median predictor autocorrelation range m",
    "Median predictor autocorrelation range km",
    "Exploratory step-08 candidate block size km",
    "Candidate-selection rule",
    "Interpretation"
  ),
  value = c(
    as.character(
      utils::packageVersion(
        "blockCV"
      )
    ),
    as.character(
      autocor_sample_n
    ),
    as.character(
      autocor_seed
    ),
    as.character(
      distance_crs
    ),
    format(
      predictor_range_m,
      digits = 12,
      scientific = FALSE
    ),
    format(
      predictor_range_m /
        1000,
      digits = 12,
      scientific = FALSE
    ),
    as.character(
      exploratory_candidate_km
    ),
    candidate_rule,
    "Exploratory predictor-based guide only; residual autocorrelation should be checked after model fitting."
  )
)

# ------------------------------------------------------------
# 7. Save outputs
# ------------------------------------------------------------

readr::write_csv(
  range_table,
  file.path(
    audit_dir,
    "09A_predictor_spatial_autocorrelation_ranges.csv"
  )
)

readr::write_csv(
  scale_comparison,
  file.path(
    audit_dir,
    "09B_block_scale_autocorrelation_comparison.csv"
  )
)

readr::write_csv(
  decision_support,
  file.path(
    audit_dir,
    "09C_autocorrelation_decision_support.csv"
  )
)

readr::write_csv(
  selected_step08_scale,
  file.path(
    audit_dir,
    "09D_exploratory_candidate_scale_summary.csv"
  )
)

readr::write_csv(
  selected_weight_summary,
  file.path(
    audit_dir,
    "09E_exploratory_candidate_weight_summary.csv"
  )
)

readr::write_csv(
  selected_cv_quality,
  file.path(
    audit_dir,
    "09F_exploratory_candidate_cv_quality.csv"
  )
)

saveRDS(
  list(
    predictor_autocorrelation =
      sac_predictors,
    predictor_range_table =
      range_table,
    median_predictor_range_m =
      predictor_range_m,
    scale_comparison =
      scale_comparison,
    exploratory_candidate_km =
      exploratory_candidate_km,
    candidate_rule =
      candidate_rule,
    selected_step08_scale =
      selected_step08_scale,
    selected_weight_summary =
      selected_weight_summary,
    selected_cv_quality =
      selected_cv_quality,
    environmental_sample_n =
      autocor_sample_n,
    environmental_sample_seed =
      autocor_seed,
    distance_crs =
      distance_crs,
    note =
      "Exploratory predictor-based guide only; no final CV scale or weighting strategy selected and no model fitted."
  ),
  file.path(
    work_dir,
    "09_REV2_spatial_autocorrelation.rds"
  )
)

# ------------------------------------------------------------
# 8. Update canonical workflow manifest
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
  "08_audit_spatial_sampling_bias.R",
  "09_estimate_spatial_autocorrelation.R",
  "",
  "IMPORTANT",
  "Step 08 quantifies geographic sampling imbalance and creates candidate block scales and weights without selecting them.",
  "Step 09 estimates predictor spatial autocorrelation ranges on a reproducible equal-area sample using blockCV.",
  "The step-09 candidate scale is exploratory: predictor autocorrelation is a pre-fitting guide, not residual autocorrelation.",
  "No weighting strategy or final cross-validation scale is selected in step 09 and no SDM is fitted."
)

writeLines(
  canonical_workflow,
  manifest_file
)

# ------------------------------------------------------------
# 9. Console report
# ------------------------------------------------------------

cat(
  "\n============================================================\n"
)

cat(
  "REV2 - PREDICTOR SPATIAL AUTOCORRELATION AUDIT\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "Predictor autocorrelation ranges returned by blockCV:\n"
)

print(
  range_table,
  row.names = FALSE
)

cat(
  "\nMedian predictor autocorrelation range (m): ",
  format(
    predictor_range_m,
    digits = 12,
    scientific = FALSE
  ),
  "\n",
  sep = ""
)

cat(
  "Median predictor autocorrelation range (km): ",
  format(
    predictor_range_m /
      1000,
    digits = 12,
    scientific = FALSE
  ),
  "\n",
  sep = ""
)

cat(
  "\nComparison with step-08 candidate scales:\n"
)

print(
  as.data.frame(
    scale_comparison
  ),
  row.names = FALSE
)

cat(
  "\nExploratory candidate block size from predictor autocorrelation: ",
  exploratory_candidate_km,
  " km\n",
  sep = ""
)

cat(
  "Rule: ",
  candidate_rule,
  "\n",
  sep = ""
)

cat(
  "\nStep-08 diagnostics for the exploratory candidate scale:\n"
)

print(
  as.data.frame(
    selected_step08_scale
  ),
  row.names = FALSE
)

cat(
  "\nCandidate weight diagnostics at that scale:\n"
)

print(
  as.data.frame(
    selected_weight_summary
  ),
  row.names = FALSE
)

cat(
  "\nCandidate CV-balance diagnostics at that scale:\n"
)

print(
  as.data.frame(
    selected_cv_quality
  ),
  row.names = FALSE
)

cat(
  "\nThis is an exploratory pre-fitting guide. ",
  "Residual spatial autocorrelation must be checked after fitting.\n",
  sep = ""
)

cat(
  "No SDM has been fitted.\n"
)

cat(
  "============================================================\n"
)
