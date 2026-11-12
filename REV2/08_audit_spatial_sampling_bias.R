# ============================================================
# ROSA RUGOSA - REVISION 2
# 08_audit_spatial_sampling_bias.R
#
# Purpose
# -------
# Quantify geographic sampling imbalance in the frozen REV2
# calibration data and prepare candidate spatial corrections for
# later evaluation.
#
# This script:
#   - summarises occurrence concentration by country;
#   - quantifies clustering on 50, 100 and 200 km equal-area grids;
#   - calculates candidate inverse-density presence weights;
#   - creates deterministic 5-fold block assignments at each scale;
#   - reports fold balance for presences and background records.
#
# Candidate weights and folds are diagnostic outputs only.
# No SDM is fitted and no weighting strategy is selected here.
# ============================================================

rm(list = ls())
gc()

required_pkgs <- c(
  "sf", "dplyr", "readr", "tibble", "tidyr"
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
    )
  )
}

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

old_08 <- c(
  list.files(
    audit_dir,
    pattern = "^08[A-Z]_",
    full.names = TRUE
  ),
  file.path(
    work_dir,
    "08_REV2_sampling_bias_audit.rds"
  ),
  file.path(
    work_dir,
    "08_model_input_spatial_folds.csv"
  ),
  file.path(
    work_dir,
    "08_presence_spatial_weights.csv"
  )
)

old_08 <- old_08[
  file.exists(
    old_08
  )
]

if (length(old_08) > 0) {
  unlink(
    old_08,
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

step07_file <- file.path(
  work_dir,
  "07_REV2_background_FINAL7_FROZEN.rds"
)

if (!file.exists(step06_file)) {
  stop(
    "Missing corrected step-06 RDS: ",
    step06_file
  )
}

if (!file.exists(step07_file)) {
  stop(
    "Missing corrected step-07 RDS: ",
    step07_file
  )
}

step06 <- readRDS(
  step06_file
)

step07 <- readRDS(
  step07_file
)

occ <- as.data.frame(
  step06$occurrences_frozen
)

model_input <- as.data.frame(
  step07$model_input
)

final7_vars <- step07$final7_vars

if (nrow(occ) !=
    as.integer(
      step06$frozen_n
    )) {
  stop(
    "Occurrence count differs from the frozen step-06 count."
  )
}

if (sum(
  model_input$pa == 1
) != nrow(occ)) {
  stop(
    "Presence count in model input differs from frozen occurrences."
  )
}

if (sum(
  model_input$pa == 0
) != 10000L) {
  stop(
    "Step-07 model input does not contain 10,000 background records."
  )
}

required_occ_cols <- c(
  "occurrence_id",
  "source",
  "country_code",
  "assigned_longitude",
  "assigned_latitude",
  "cell_id"
)

missing_occ_cols <- setdiff(
  required_occ_cols,
  names(occ)
)

if (length(missing_occ_cols) > 0) {
  stop(
    "Frozen occurrence table lacks: ",
    paste(
      missing_occ_cols,
      collapse = ", "
    )
  )
}

required_model_cols <- c(
  "pa",
  "cell_id",
  "decimal_longitude",
  "decimal_latitude",
  "source",
  "occurrence_id",
  final7_vars
)

missing_model_cols <- setdiff(
  required_model_cols,
  names(model_input)
)

if (length(missing_model_cols) > 0) {
  stop(
    "Model-input table lacks: ",
    paste(
      missing_model_cols,
      collapse = ", "
    )
  )
}

# ------------------------------------------------------------
# 3. Equal-area coordinates
# ------------------------------------------------------------

distance_crs <- as.integer(
  step06$distance_crs
)

if (!is.finite(distance_crs)) {
  stop(
    "Missing valid projected CRS in step 06."
  )
}

model_sf <- sf::st_as_sf(
  model_input,
  coords = c(
    "decimal_longitude",
    "decimal_latitude"
  ),
  crs = 4326,
  remove = FALSE
)

model_m <- sf::st_transform(
  model_sf,
  distance_crs
)

xy_m <- sf::st_coordinates(
  model_m
)

model_input$x_m <- as.numeric(
  xy_m[
    ,
    1
  ]
)

model_input$y_m <- as.numeric(
  xy_m[
    ,
    2
  ]
)

if (any(
  !is.finite(
    model_input$x_m
  ) |
    !is.finite(
      model_input$y_m
    )
)) {
  stop(
    "Non-finite projected coordinates detected."
  )
}

# ------------------------------------------------------------
# 4. Country-level sampling concentration
# ------------------------------------------------------------

occ$country_code <- toupper(
  trimws(
    as.character(
      occ$country_code
    )
  )
)

occ$country_code[
  is.na(
    occ$country_code
  ) |
    occ$country_code == ""
] <- "UNRESOLVED"

country_counts <- occ %>%
  count(
    country_code,
    name = "n_presences",
    sort = TRUE
  ) %>%
  mutate(
    proportion =
      n_presences /
      sum(
        n_presences
      ),
    cumulative_proportion =
      cumsum(
        proportion
      ),
    rank =
      row_number()
  )

country_p <- country_counts$proportion

country_hhi <- sum(
  country_p^2
)

country_effective_n <- ifelse(
  country_hhi > 0,
  1 / country_hhi,
  NA_real_
)

country_shannon <- -sum(
  country_p *
    log(
      country_p
    )
)

top1_share <- sum(
  head(
    country_p,
    1
  )
)

top5_share <- sum(
  head(
    country_p,
    5
  )
)

top10_share <- sum(
  head(
    country_p,
    10
  )
)

country_metrics <- tibble(
  metric = c(
    "Total frozen presences",
    "Countries or territories represented",
    "Largest-country share",
    "Top-5 country share",
    "Top-10 country share",
    "Country HHI",
    "Effective number of equally represented countries",
    "Country Shannon entropy"
  ),
  value = c(
    nrow(
      occ
    ),
    nrow(
      country_counts
    ),
    top1_share,
    top5_share,
    top10_share,
    country_hhi,
    country_effective_n,
    country_shannon
  )
)

# ------------------------------------------------------------
# 5. Spatial block definitions
# ------------------------------------------------------------

block_scales_km <- c(
  50,
  100,
  200
)

k_folds <- 5L

make_block_id <- function(
    x,
    y,
    block_km
) {

  block_m <- block_km *
    1000

  bx <- floor(
    x /
      block_m
  )

  by <- floor(
    y /
      block_m
  )

  paste0(
    bx,
    "_",
    by
  )
}

# Greedy deterministic allocation of complete spatial blocks to folds.
# Blocks remain indivisible. The objective balances both presences and
# background records across the five folds.
assign_blocks_to_folds <- function(
    block_table,
    k = 5L
) {

  tab <- block_table %>%
    arrange(
      desc(
        n_presence
      ),
      desc(
        n_background
      ),
      block_id
    )

  fold_presence <- rep(
    0,
    k
  )

  fold_background <- rep(
    0,
    k
  )

  fold_blocks <- rep(
    0,
    k
  )

  target_presence <- sum(
    tab$n_presence
  ) /
    k

  target_background <- sum(
    tab$n_background
  ) /
    k

  assignment <- integer(
    nrow(
      tab
    )
  )

  for (i in seq_len(
    nrow(
      tab
    )
  )) {

    scores <- rep(
      NA_real_,
      k
    )

    for (f in seq_len(
      k
    )) {

      p_new <- fold_presence
      b_new <- fold_background
      k_new <- fold_blocks

      p_new[f] <- p_new[f] +
        tab$n_presence[i]

      b_new[f] <- b_new[f] +
        tab$n_background[i]

      k_new[f] <- k_new[f] +
        1

      p_term <- if (
        target_presence > 0
      ) {
        sum(
          (
            p_new -
              target_presence
          )^2
        ) /
          target_presence^2
      } else {
        0
      }

      b_term <- if (
        target_background > 0
      ) {
        sum(
          (
            b_new -
              target_background
          )^2
        ) /
          target_background^2
      } else {
        0
      }

      block_term <- 0.01 *
        sum(
          (
            k_new -
              mean(
                k_new
              )
          )^2
        )

      scores[f] <- p_term +
        b_term +
        block_term
    }

    best <- which(
      scores ==
        min(
          scores
        )
    )

    if (length(best) > 1) {

      total_load <- (
        fold_presence[
          best
        ] /
          max(
            1,
            target_presence
          )
      ) +
        (
          fold_background[
            best
          ] /
            max(
              1,
              target_background
            )
        )

      best <- best[
        which.min(
          total_load
        )
      ]
    }

    f <- best[1]

    assignment[i] <- f

    fold_presence[f] <-
      fold_presence[f] +
      tab$n_presence[i]

    fold_background[f] <-
      fold_background[f] +
      tab$n_background[i]

    fold_blocks[f] <-
      fold_blocks[f] +
      1
  }

  tab$fold <- assignment

  tab
}

scale_summaries <- list()
block_tables <- list()
fold_tables <- list()
weight_tables <- list()
model_fold_tables <- list()
weight_metric_tables <- list()

for (block_km in block_scales_km) {

  block_col <- paste0(
    "block_",
    block_km,
    "km"
  )

  model_input[[block_col]] <-
    make_block_id(
      model_input$x_m,
      model_input$y_m,
      block_km
    )

  block_counts <- model_input %>%
    group_by(
      block_id =
        .data[[block_col]]
    ) %>%
    summarise(
      n_presence =
        sum(
          pa == 1
        ),
      n_background =
        sum(
          pa == 0
        ),
      n_total =
        n(),
      .groups = "drop"
    )

  presence_blocks <- block_counts %>%
    filter(
      n_presence > 0
    ) %>%
    arrange(
      desc(
        n_presence
      )
    )

  presence_counts <- presence_blocks$n_presence

  n_presence_blocks <- nrow(
    presence_blocks
  )

  max_block_share <- max(
    presence_counts
  ) /
    nrow(
      occ
    )

  top10_n <- max(
    1L,
    ceiling(
      0.10 *
        n_presence_blocks
    )
  )

  top10_block_share <- sum(
    head(
      sort(
        presence_counts,
        decreasing = TRUE
      ),
      top10_n
    )
  ) /
    nrow(
      occ
    )

  singleton_share <- mean(
    presence_counts ==
      1
  )

  block_hhi <- sum(
    (
      presence_counts /
        sum(
          presence_counts
        )
    )^2
  )

  effective_blocks <- ifelse(
    block_hhi > 0,
    1 /
      block_hhi,
    NA_real_
  )

  scale_summaries[[
    as.character(
      block_km
    )
  ]] <- tibble(
    block_km = block_km,
    n_blocks_total =
      nrow(
        block_counts
      ),
    n_presence_blocks =
      n_presence_blocks,
    n_background_blocks =
      sum(
        block_counts$n_background >
          0
      ),
    n_blocks_with_both =
      sum(
        block_counts$n_presence >
          0 &
          block_counts$n_background >
          0
      ),
    mean_presences_per_occupied_block =
      mean(
        presence_counts
      ),
    median_presences_per_occupied_block =
      stats::median(
        presence_counts
      ),
    max_presences_in_one_block =
      max(
        presence_counts
      ),
    max_block_presence_share =
      max_block_share,
    top_10pct_blocks_presence_share =
      top10_block_share,
    singleton_presence_block_share =
      singleton_share,
    presence_block_HHI =
      block_hhi,
    effective_number_presence_blocks =
      effective_blocks
  )

  block_counts$block_km <-
    block_km

  block_tables[[
    as.character(
      block_km
    )
  ]] <- block_counts

  # Candidate presence weights.
  presence_weights <- model_input %>%
    filter(
      pa == 1
    ) %>%
    select(
      occurrence_id,
      cell_id,
      block_id =
        all_of(
          block_col
        )
    ) %>%
    left_join(
      presence_blocks %>%
        select(
          block_id,
          n_presence
        ),
      by = "block_id"
    ) %>%
    mutate(
      block_km =
        block_km,
      weight_inverse_raw =
        1 /
        n_presence,
      weight_sqrt_inverse_raw =
        1 /
        sqrt(
          n_presence
        ),
      weight_inverse =
        weight_inverse_raw /
        mean(
          weight_inverse_raw
        ),
      weight_sqrt_inverse =
        weight_sqrt_inverse_raw /
        mean(
          weight_sqrt_inverse_raw
        )
    )

  weight_tables[[
    as.character(
      block_km
    )
  ]] <- presence_weights

  weight_metrics <- bind_rows(
    lapply(
      c(
        "weight_inverse",
        "weight_sqrt_inverse"
      ),
      function(w) {

        x <- presence_weights[[
          w
        ]]

        tibble(
          block_km =
            block_km,
          weight_scheme =
            w,
          min_weight =
            min(
              x
            ),
          median_weight =
            stats::median(
              x
            ),
          max_weight =
            max(
              x
            ),
          effective_sample_size =
            (
              sum(
                x
              )^2
            ) /
            sum(
              x^2
            )
        )
      }
    )
  )

  weight_metric_tables[[
    as.character(
      block_km
    )
  ]] <- weight_metrics

  # Candidate five-fold block allocation.
  assigned_blocks <- assign_blocks_to_folds(
    block_counts,
    k =
      k_folds
  )

  assigned_blocks$block_km <-
    block_km

  fold_tables[[
    as.character(
      block_km
    )
  ]] <- assigned_blocks

  fold_lookup <- assigned_blocks %>%
    select(
      block_id,
      fold
    )

  fold_col <- paste0(
    "fold_",
    block_km,
    "km"
  )

  fold_vec <- tibble(
    row_id =
      seq_len(
        nrow(
          model_input
        )
      ),
    block_id =
      model_input[[
        block_col
      ]]
  ) %>%
    left_join(
      fold_lookup,
      by = "block_id"
    ) %>%
    arrange(
      row_id
    )

  if (any(
    is.na(
      fold_vec$fold
    )
  )) {
    stop(
      "Missing fold assignment at ",
      block_km,
      " km."
    )
  }

  model_input[[fold_col]] <-
    as.integer(
      fold_vec$fold
    )

  fold_balance <- model_input %>%
    group_by(
      fold =
        .data[[fold_col]]
    ) %>%
    summarise(
      n_presence =
        sum(
          pa == 1
        ),
      n_background =
        sum(
          pa == 0
        ),
      n_total =
        n(),
      .groups = "drop"
    ) %>%
    mutate(
      block_km =
        block_km
    )

  model_fold_tables[[
    as.character(
      block_km
    )
  ]] <- fold_balance
}

spatial_scale_summary <- bind_rows(
  scale_summaries
) %>%
  arrange(
    block_km
  )

all_block_counts <- bind_rows(
  block_tables
) %>%
  arrange(
    block_km,
    desc(
      n_presence
    ),
    block_id
  )

all_weights <- bind_rows(
  weight_tables
) %>%
  arrange(
    block_km,
    cell_id
  )

weight_summary <- bind_rows(
  weight_metric_tables
) %>%
  arrange(
    block_km,
    weight_scheme
  )

fold_balance <- bind_rows(
  model_fold_tables
) %>%
  arrange(
    block_km,
    fold
  )

all_fold_blocks <- bind_rows(
  fold_tables
) %>%
  arrange(
    block_km,
    fold,
    block_id
  )

# ------------------------------------------------------------
# 6. Fold-balance diagnostics
# ------------------------------------------------------------

fold_scale_quality <- fold_balance %>%
  group_by(
    block_km
  ) %>%
  summarise(
    min_presence_per_fold =
      min(
        n_presence
      ),
    max_presence_per_fold =
      max(
        n_presence
      ),
    presence_fold_CV =
      stats::sd(
        n_presence
      ) /
      mean(
        n_presence
      ),
    min_background_per_fold =
      min(
        n_background
      ),
    max_background_per_fold =
      max(
        n_background
      ),
    background_fold_CV =
      stats::sd(
        n_background
      ) /
      mean(
        n_background
      ),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 7. Save outputs
# ------------------------------------------------------------

readr::write_csv(
  country_counts,
  file.path(
    audit_dir,
    "08A_country_sampling_concentration.csv"
  )
)

readr::write_csv(
  country_metrics,
  file.path(
    audit_dir,
    "08B_country_bias_metrics.csv"
  )
)

readr::write_csv(
  spatial_scale_summary,
  file.path(
    audit_dir,
    "08C_spatial_block_scale_summary.csv"
  )
)

readr::write_csv(
  weight_summary,
  file.path(
    audit_dir,
    "08D_candidate_weight_summary.csv"
  )
)

readr::write_csv(
  fold_balance,
  file.path(
    audit_dir,
    "08E_candidate_cv_fold_balance.csv"
  )
)

readr::write_csv(
  fold_scale_quality,
  file.path(
    audit_dir,
    "08F_candidate_cv_scale_quality.csv"
  )
)

readr::write_csv(
  all_block_counts,
  file.path(
    audit_dir,
    "08G_spatial_block_counts.csv"
  )
)

readr::write_csv(
  all_fold_blocks,
  file.path(
    audit_dir,
    "08H_spatial_block_fold_assignment.csv"
  )
)

readr::write_csv(
  all_weights,
  file.path(
    work_dir,
    "08_presence_spatial_weights.csv"
  )
)

model_input_out <- model_input %>%
  select(
    -x_m,
    -y_m
  )

readr::write_csv(
  model_input_out,
  file.path(
    work_dir,
    "08_model_input_spatial_folds.csv"
  )
)

saveRDS(
  list(
    country_counts =
      country_counts,
    country_metrics =
      country_metrics,
    spatial_scale_summary =
      spatial_scale_summary,
    candidate_weight_summary =
      weight_summary,
    candidate_cv_fold_balance =
      fold_balance,
    candidate_cv_scale_quality =
      fold_scale_quality,
    block_counts =
      all_block_counts,
    block_fold_assignment =
      all_fold_blocks,
    presence_spatial_weights =
      all_weights,
    model_input_spatial_folds =
      model_input_out,
    block_scales_km =
      block_scales_km,
    k_folds =
      k_folds,
    distance_crs =
      distance_crs,
    note =
      "Diagnostic only: no weighting scheme or CV scale selected and no model fitted."
  ),
  file.path(
    work_dir,
    "08_REV2_sampling_bias_audit.rds"
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
  "",
  "IMPORTANT",
  "Step 05 defines the corrected 1,000-m terrestrial European coastal domain and corrected occurrence assignment in EPSG:3035.",
  "Step 06 freezes corrected GBIF cells plus four independently verified Brussa field cells.",
  "Step 07 samples 10,000 corrected-domain background cells outside 1,000 m from frozen presences.",
  "Step 08 audits geographic sampling imbalance at 50, 100 and 200 km and prepares candidate spatial weights and five-fold block assignments.",
  "Step 08 does not select a weighting strategy or CV scale and does not fit any model."
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
  "REV2 - SPATIAL SAMPLING-BIAS AUDIT\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "Country concentration metrics:\n"
)

print(
  as.data.frame(
    country_metrics
  ),
  row.names = FALSE
)

cat(
  "\nTop 15 countries or territories by frozen presence count:\n"
)

print(
  as.data.frame(
    head(
      country_counts,
      15
    )
  ),
  row.names = FALSE
)

cat(
  "\nSpatial block-scale summary:\n"
)

print(
  as.data.frame(
    spatial_scale_summary
  ),
  row.names = FALSE
)

cat(
  "\nCandidate weighting summary:\n"
)

print(
  as.data.frame(
    weight_summary
  ),
  row.names = FALSE
)

cat(
  "\nCandidate 5-fold CV balance:\n"
)

print(
  as.data.frame(
    fold_balance
  ),
  row.names = FALSE
)

cat(
  "\nCandidate CV scale-quality summary:\n"
)

print(
  as.data.frame(
    fold_scale_quality
  ),
  row.names = FALSE
)

cat(
  "\nDiagnostic only: no weighting scheme or CV scale has been selected.\n"
)

cat(
  "No SDM has been fitted.\n"
)

cat(
  "============================================================\n"
)
