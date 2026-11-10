# ============================================================
# ROSA RUGOSA - REVISION 2
# 06_freeze_REV2_occurrence_dataset.R
#
# Purpose
# -------
# Freeze the definitive invaded-range occurrence dataset for REV2.
#
# Calibration presences:
#   1) 8,705 GBIF cells rebuilt under the FINAL7 coastal mask;
#   2) 4 independent Brussa field cells from the Sara/Ribotta
#      dataset, confirmed as >1 km from all GBIF records passing
#      the declared filters and absent from the FINAL7 GBIF cells.
#
# Other Sara/Ribotta candidate cells are retained only as audit
# information and are NOT added to calibration.
#
# No background is sampled.
# No model is fitted.
# No REV1 file is modified.
# ============================================================

rm(list = ls())
gc()

required_pkgs <- c(
  "dplyr", "readr", "tibble", "tidyr"
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

library(dplyr)
library(readr)
library(tibble)
library(tidyr)

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

old_06 <- c(
  list.files(
    audit_dir,
    pattern = "^06[A-Z]_",
    full.names = TRUE
  ),
  file.path(
    work_dir,
    "06_REV2_occurrence_dataset_FROZEN.rds"
  ),
  file.path(
    work_dir,
    "occurrences_REV2_FROZEN.csv"
  )
)

old_06 <- old_06[
  file.exists(old_06)
]

if (length(old_06) > 0) {
  unlink(
    old_06,
    force = TRUE
  )
}

# ------------------------------------------------------------
# 2. Load canonical step-05 output
# ------------------------------------------------------------

step05_file <- file.path(
  work_dir,
  "05_REV2_occurrences_FINAL7_VALID.rds"
)

if (!file.exists(step05_file)) {
  stop(
    "Missing canonical step-05 RDS: ",
    step05_file,
    "\nRun 05_rebuild_REV2_occurrences_FINAL7.R first."
  )
}

step05 <- readRDS(
  step05_file
)

required_step05_objects <- c(
  "gbif_final7",
  "candidate_independent_cells",
  "sara_point_audit",
  "final7_vars",
  "n_final7_cells"
)

missing_objects <- required_step05_objects[
  !vapply(
    required_step05_objects,
    function(x) {
      !is.null(step05[[x]])
    },
    logical(1)
  )
]

if (length(missing_objects) > 0) {
  stop(
    "Step-05 RDS lacks: ",
    paste(
      missing_objects,
      collapse = ", "
    )
  )
}

gbif_final7 <- as.data.frame(
  step05$gbif_final7
)

candidate_cells <- as.data.frame(
  step05$candidate_independent_cells
)

sara_points <- as.data.frame(
  step05$sara_point_audit
)

final7_vars <- step05$final7_vars
n_final7_cells <- step05$n_final7_cells

# ------------------------------------------------------------
# 3. Validate canonical GBIF FINAL7 set
# ------------------------------------------------------------

required_gbif_cols <- c(
  "gbif_id",
  "cell_id",
  "country_code",
  "orig_decimal_longitude",
  "orig_decimal_latitude",
  "decimal_longitude",
  "decimal_latitude",
  "snap_distance_m",
  "snap_type"
)

missing_gbif_cols <- setdiff(
  required_gbif_cols,
  names(gbif_final7)
)

if (length(missing_gbif_cols) > 0) {
  stop(
    "Canonical GBIF FINAL7 set lacks: ",
    paste(
      missing_gbif_cols,
      collapse = ", "
    )
  )
}

gbif_final7$gbif_id <- as.character(
  gbif_final7$gbif_id
)

gbif_final7$cell_id <- as.integer(
  gbif_final7$cell_id
)

if (nrow(gbif_final7) != 8705L) {
  stop(
    "Expected 8,705 GBIF FINAL7 cells, found ",
    nrow(gbif_final7),
    ". Do not freeze dataset."
  )
}

if (dplyr::n_distinct(
  gbif_final7$cell_id
) != 8705L) {
  stop(
    "GBIF FINAL7 cell_id values are not unique."
  )
}

if (n_final7_cells != 210748L) {
  stop(
    "Expected FINAL7 mask size 210,748 cells, found ",
    n_final7_cells,
    ". Do not freeze dataset."
  )
}

# ------------------------------------------------------------
# 4. Select only independent Brussa field cells
# ------------------------------------------------------------

required_candidate_cols <- c(
  "sara_id",
  "sara_lon",
  "sara_lat",
  "country_code",
  "country_name",
  "assigned_cell_id",
  "snap_distance_m",
  "snap_type",
  "retained_FINAL7",
  "cell_already_in_GBIF_FINAL7",
  "nearest_DECLARED_GBIF_distance_m",
  "nearest_DECLARED_GBIF_id",
  "candidate_independent_field_point"
)

missing_candidate_cols <- setdiff(
  required_candidate_cols,
  names(candidate_cells)
)

if (length(missing_candidate_cols) > 0) {
  stop(
    "Candidate-cell table lacks: ",
    paste(
      missing_candidate_cols,
      collapse = ", "
    )
  )
}

brussa_cells <- candidate_cells %>%
  filter(
    country_code == "IT",
    retained_FINAL7,
    !cell_already_in_GBIF_FINAL7,
    candidate_independent_field_point,
    is.finite(
      nearest_DECLARED_GBIF_distance_m
    ),
    nearest_DECLARED_GBIF_distance_m > 1000
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

expected_brussa_cells <- c(
  22794523L,
  22801719L,
  22801720L,
  22801723L
)

if (nrow(brussa_cells) != 4L) {
  stop(
    "Expected 4 independent Brussa cells, found ",
    nrow(brussa_cells),
    ". Do not freeze dataset."
  )
}

if (!setequal(
  as.integer(
    brussa_cells$assigned_cell_id
  ),
  expected_brussa_cells
)) {
  stop(
    "Independent Italian candidate cells do not match the audited Brussa cell IDs."
  )
}

if (any(
  brussa_cells$assigned_cell_id %in%
    gbif_final7$cell_id
)) {
  stop(
    "At least one Brussa field cell overlaps the canonical GBIF FINAL7 set."
  )
}

# Keep all seven independent Italian observations as provenance evidence.
brussa_point_audit <- sara_points %>%
  filter(
    country_code == "IT",
    retained_FINAL7,
    !cell_already_in_GBIF_FINAL7,
    candidate_independent_field_point,
    is.finite(
      nearest_DECLARED_GBIF_distance_m
    ),
    nearest_DECLARED_GBIF_distance_m > 1000
  ) %>%
  arrange(
    assigned_cell_id,
    sara_id
  )

if (nrow(brussa_point_audit) != 7L) {
  stop(
    "Expected 7 independent Italian field observations, found ",
    nrow(brussa_point_audit),
    ". Do not freeze dataset."
  )
}

# ------------------------------------------------------------
# 5. Standardise calibration occurrence schema
# ------------------------------------------------------------

gbif_frozen <- gbif_final7 %>%
  transmute(
    occurrence_id =
      paste0(
        "GBIF_",
        gbif_id
      ),

    source =
      "GBIF",

    source_record_id =
      gbif_id,

    source_dataset =
      "GBIF download used in REV1/REV2 audit",

    country_code =
      as.character(
        country_code
      ),

    country_name =
      if ("country" %in% names(gbif_final7)) {
        as.character(country)
      } else {
        NA_character_
      },

    original_longitude =
      as.numeric(
        orig_decimal_longitude
      ),

    original_latitude =
      as.numeric(
        orig_decimal_latitude
      ),

    assigned_longitude =
      as.numeric(
        decimal_longitude
      ),

    assigned_latitude =
      as.numeric(
        decimal_latitude
      ),

    cell_id =
      as.integer(
        cell_id
      ),

    snap_distance_m =
      as.numeric(
        snap_distance_m
      ),

    snap_type =
      as.character(
        snap_type
      ),

    nearest_declared_GBIF_distance_m =
      0,

    nearest_declared_GBIF_id =
      gbif_id,

    field_observation =
      FALSE,

    calibration_presence =
      TRUE,

    provenance_note =
      "GBIF occurrence passing declared REV2 filters; one representative per FINAL7 coastal raster cell"
  )

brussa_frozen <- brussa_cells %>%
  transmute(
    occurrence_id =
      paste0(
        "BRUSSA_CELL_",
        assigned_cell_id
      ),

    source =
      "Brussa_field",

    source_record_id =
      paste0(
        "SARA_",
        sara_id
      ),

    source_dataset =
      "Sara/Ribotta field dataset - Brussa",

    country_code =
      "IT",

    country_name =
      "Italy",

    original_longitude =
      as.numeric(
        sara_lon
      ),

    original_latitude =
      as.numeric(
        sara_lat
      ),

    assigned_longitude =
      as.numeric(
        sara_lon
      ),

    assigned_latitude =
      as.numeric(
        sara_lat
      ),

    cell_id =
      as.integer(
        assigned_cell_id
      ),

    snap_distance_m =
      as.numeric(
        snap_distance_m
      ),

    snap_type =
      as.character(
        snap_type
      ),

    nearest_declared_GBIF_distance_m =
      as.numeric(
        nearest_DECLARED_GBIF_distance_m
      ),

    nearest_declared_GBIF_id =
      as.character(
        nearest_DECLARED_GBIF_id
      ),

    field_observation =
      TRUE,

    calibration_presence =
      TRUE,

    provenance_note =
      "Independent Brussa field observation; cell absent from FINAL7 GBIF set and >1 km from all GBIF records passing declared filters"
  )

occ_frozen <- bind_rows(
  gbif_frozen,
  brussa_frozen
) %>%
  arrange(
    cell_id,
    source,
    occurrence_id
  )

# ------------------------------------------------------------
# 6. Hard validation of frozen occurrence set
# ------------------------------------------------------------

if (nrow(occ_frozen) != 8709L) {
  stop(
    "Frozen occurrence set should contain 8,709 rows, found ",
    nrow(occ_frozen),
    "."
  )
}

if (dplyr::n_distinct(
  occ_frozen$cell_id
) != 8709L) {
  stop(
    "Frozen occurrence set does not contain exactly one presence per raster cell."
  )
}

if (anyDuplicated(
  occ_frozen$occurrence_id
) > 0) {
  stop(
    "Frozen occurrence_id values are not unique."
  )
}

if (sum(
  occ_frozen$source == "GBIF"
) != 8705L) {
  stop(
    "Frozen set does not contain exactly 8,705 GBIF presences."
  )
}

if (sum(
  occ_frozen$source == "Brussa_field"
) != 4L) {
  stop(
    "Frozen set does not contain exactly 4 Brussa field cells."
  )
}

if (!all(
  occ_frozen$calibration_presence
)) {
  stop(
    "Frozen occurrence set contains non-calibration records."
  )
}

# ------------------------------------------------------------
# 7. Audit tables
# ------------------------------------------------------------

source_summary <- occ_frozen %>%
  count(
    source,
    name = "n_cells"
  ) %>%
  arrange(
    source
  )

country_summary <- occ_frozen %>%
  mutate(
    country_code =
      ifelse(
        is.na(country_code) |
          country_code == "",
        "UNRESOLVED",
        country_code
      )
  ) %>%
  count(
    country_code,
    source,
    name = "n_cells"
  ) %>%
  arrange(
    country_code,
    source
  )

freeze_summary <- tibble(
  metric = c(
    "FINAL7 valid coastal cells",
    "Frozen GBIF calibration cells",
    "Independent Brussa field observations",
    "Independent Brussa field cells",
    "Frozen total calibration presences",
    "Unique frozen raster cells",
    "Other Sara/Ribotta candidate cells excluded from calibration"
  ),

  value = c(
    n_final7_cells,
    sum(
      occ_frozen$source == "GBIF"
    ),
    nrow(
      brussa_point_audit
    ),
    sum(
      occ_frozen$source == "Brussa_field"
    ),
    nrow(
      occ_frozen
    ),
    dplyr::n_distinct(
      occ_frozen$cell_id
    ),
    nrow(
      candidate_cells
    ) -
      nrow(
        brussa_cells
      )
  )
)

# ------------------------------------------------------------
# 8. Save canonical frozen outputs
# ------------------------------------------------------------

readr::write_csv(
  occ_frozen,
  file.path(
    work_dir,
    "occurrences_REV2_FROZEN.csv"
  )
)

readr::write_csv(
  freeze_summary,
  file.path(
    audit_dir,
    "06A_frozen_occurrence_summary.csv"
  )
)

readr::write_csv(
  source_summary,
  file.path(
    audit_dir,
    "06B_frozen_occurrence_source_counts.csv"
  )
)

readr::write_csv(
  country_summary,
  file.path(
    audit_dir,
    "06C_frozen_occurrence_country_counts.csv"
  )
)

readr::write_csv(
  brussa_point_audit,
  file.path(
    audit_dir,
    "06D_Brussa_independent_point_audit.csv"
  )
)

readr::write_csv(
  brussa_frozen,
  file.path(
    audit_dir,
    "06E_Brussa_cells_added_to_calibration.csv"
  )
)

other_sara_candidates <- candidate_cells %>%
  filter(
    !(
      country_code == "IT" &
        assigned_cell_id %in%
          brussa_cells$assigned_cell_id
    )
  ) %>%
  mutate(
    included_in_calibration =
      FALSE,

    exclusion_from_calibration_reason =
      "Independent provenance not confirmed as Brussa field survey; retained as audit only"
  )

readr::write_csv(
  other_sara_candidates,
  file.path(
    audit_dir,
    "06F_other_Sara_candidate_cells_audit_only.csv"
  )
)

saveRDS(
  list(
    occurrences_frozen = occ_frozen,
    gbif_frozen = gbif_frozen,
    brussa_frozen = brussa_frozen,
    brussa_point_audit = brussa_point_audit,
    other_sara_candidates = other_sara_candidates,
    source_summary = source_summary,
    country_summary = country_summary,
    freeze_summary = freeze_summary,
    final7_vars = final7_vars,
    n_final7_cells = n_final7_cells,
    frozen_n = nrow(occ_frozen)
  ),
  file.path(
    work_dir,
    "06_REV2_occurrence_dataset_FROZEN.rds"
  )
)

# ------------------------------------------------------------
# 9. Update canonical workflow manifest
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
  "",
  "SUPERSEDED / INVALIDATED",
  "03_reconcile_Ribotta_with_REV1_CORRECTED.R",
  "04_trace_missing_Ribotta_cells.R",
  "05_audit_Ribotta_geography_by_country.R",
  "",
  "IMPORTANT",
  "REV2 occurrence reconstruction uses the quality filters declared in the manuscript and the seven-predictor current coastal mask.",
  "Step 05 established the canonical 8,705-cell GBIF FINAL7 set and audited Sara/Ribotta independence before spatial thinning.",
  "Step 06 freezes the calibration occurrence set at 8,709 unique cells: 8,705 GBIF + 4 independent Brussa field cells.",
  "Other Sara/Ribotta candidate cells remain audit-only and are not used for calibration.",
  "No background is sampled and no model is fitted in steps 00-06."
)

writeLines(
  canonical_workflow,
  manifest_file
)

# ------------------------------------------------------------
# 10. Console report
# ------------------------------------------------------------

cat(
  "\n============================================================\n"
)

cat(
  "REV2 - OCCURRENCE DATASET FROZEN\n"
)

cat(
  "============================================================\n\n"
)

print(
  as.data.frame(
    freeze_summary
  ),
  row.names = FALSE
)

cat(
  "\nSources in frozen calibration set:\n"
)

print(
  as.data.frame(
    source_summary
  ),
  row.names = FALSE
)

cat(
  "\nBrussa cells added to calibration:\n"
)

print(
  as.data.frame(
    brussa_frozen %>%
      select(
        occurrence_id,
        source_record_id,
        original_longitude,
        original_latitude,
        cell_id,
        nearest_declared_GBIF_distance_m,
        nearest_declared_GBIF_id
      )
  ),
  row.names = FALSE
)

cat(
  "\nFrozen occurrence file:\n",
  file.path(
    work_dir,
    "occurrences_REV2_FROZEN.csv"
  ),
  "\n"
)

cat(
  "\nNo background has been sampled. No model has been fitted.\n"
)

cat(
  "============================================================\n"
)
