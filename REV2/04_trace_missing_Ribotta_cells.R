# ============================================================
# ROSA RUGOSA - REVISION 2
# 04_trace_missing_Ribotta_cells.R
#
# Purpose
# -------
# Explain why Sara/Ribotta coastal cells that are absent from the
# final REV1 invaded dataset did not reach the final occurrence set.
# The script:
#   1) fixes country assignment using the nearest GBIF source record;
#   2) traces the nearest GBIF record through flagged -> clean ->
#      coastal -> main_occ -> occ;
#   3) identifies the likely exclusion stage/reason;
#   4) isolates the genuinely independent Sara/Ribotta cells for
#      manual inspection (e.g. field observations not in GBIF).
#
# No occurrence is added.
# No model is fitted.
# No REV1 file is modified.
# ============================================================

rm(list = ls())
gc()

required_pkgs <- c("dplyr", "readr", "tibble", "sf")
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
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
library(sf)

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

reconciliation_file <- file.path(
  work_dir,
  "03_Ribotta_REV1_reconciliation.rds"
)

if (!file.exists(reconciliation_file)) {
  stop(
    "Missing reconciliation object: ",
    reconciliation_file,
    "\nRun 03_reconcile_Ribotta_with_REV1.R first."
  )
}

rdata_revision <- file.path(
  rev1_dir,
  "dati_per_revisione.RData"
)

if (!file.exists(rdata_revision)) {
  stop("Cannot find REV1 RData: ", rdata_revision)
}

# ------------------------------------------------------------
# 2. Load objects
# ------------------------------------------------------------

rec <- readRDS(reconciliation_file)

point_audit <- rec$point_audit
cell_audit  <- rec$cell_audit

rev1 <- new.env(parent = baseenv())
load(rdata_revision, envir = rev1)

get_rev1 <- function(name, required = TRUE) {
  if (!exists(name, envir = rev1, inherits = FALSE)) {
    if (required) stop("Required REV1 object not found: ", name)
    return(NULL)
  }
  get(name, envir = rev1, inherits = FALSE)
}

flagged  <- get_rev1("flagged")
clean    <- get_rev1("clean")
coastal  <- get_rev1("coastal")
main_occ <- get_rev1("main_occ")
occ      <- get_rev1("occ")

# ------------------------------------------------------------
# 3. Helpers
# ------------------------------------------------------------

pick_col <- function(df, candidates, required = FALSE) {

  nms <- names(df)

  for (candidate in candidates) {
    hit <- which(tolower(nms) == tolower(candidate))
    if (length(hit) > 0) return(nms[hit[1]])
  }

  for (candidate in candidates) {
    hit <- grep(candidate, nms, ignore.case = TRUE, value = TRUE)
    if (length(hit) > 0) return(hit[1])
  }

  if (required) {
    stop(
      "Could not find required field. Candidates: ",
      paste(candidates, collapse = ", "),
      "\nAvailable fields: ",
      paste(nms, collapse = ", ")
    )
  }

  NA_character_
}

first_existing_value <- function(df, row, candidates) {

  cl <- pick_col(df, candidates, required = FALSE)

  if (is.na(cl)) return(NA_character_)

  as.character(df[[cl]][row])
}

normalize_country <- function(x) {

  y <- toupper(trimws(as.character(x)))

  # ISO-2 already present
  iso2 <- grepl("^[A-Z]{2}$", y)
  out <- rep(NA_character_, length(y))
  out[iso2] <- y[iso2]

  map <- c(
    "UNITED KINGDOM" = "GB",
    "GREAT BRITAIN" = "GB",
    "GBR" = "GB",
    "UK" = "GB",

    "FRANCE" = "FR",
    "FRA" = "FR",

    "ITALY" = "IT",
    "ITALIA" = "IT",
    "ITA" = "IT",

    "SPAIN" = "ES",
    "ESPANA" = "ES",
    "ESPAÑA" = "ES",
    "ESP" = "ES",

    "PORTUGAL" = "PT",
    "PRT" = "PT",

    "GREECE" = "GR",
    "GRC" = "GR",

    "CROATIA" = "HR",
    "HRV" = "HR",

    "SLOVENIA" = "SI",
    "SVN" = "SI",

    "MONTENEGRO" = "ME",
    "MNE" = "ME",

    "ALBANIA" = "AL",
    "ALB" = "AL",

    "BOSNIA AND HERZEGOVINA" = "BA",
    "BIH" = "BA",

    "TURKEY" = "TR",
    "TÜRKIYE" = "TR",
    "TURKIYE" = "TR",
    "TUR" = "TR",

    "BULGARIA" = "BG",
    "BGR" = "BG",

    "ROMANIA" = "RO",
    "ROU" = "RO",

    "CYPRUS" = "CY",
    "CYP" = "CY",

    "MALTA" = "MT",
    "MLT" = "MT",

    "GERMANY" = "DE",
    "DEU" = "DE",

    "DENMARK" = "DK",
    "DNK" = "DK",

    "NORWAY" = "NO",
    "NOR" = "NO",

    "SWEDEN" = "SE",
    "SWE" = "SE",

    "FINLAND" = "FI",
    "FIN" = "FI",

    "POLAND" = "PL",
    "POL" = "PL",

    "NETHERLANDS" = "NL",
    "NLD" = "NL",

    "BELGIUM" = "BE",
    "BEL" = "BE",

    "IRELAND" = "IE",
    "IRL" = "IE",

    "ESTONIA" = "EE",
    "EST" = "EE",

    "LATVIA" = "LV",
    "LVA" = "LV",

    "LITHUANIA" = "LT",
    "LTU" = "LT",

    "ICELAND" = "IS",
    "ISL" = "IS",

    "UKRAINE" = "UA",
    "UKR" = "UA",

    "GEORGIA" = "GE",
    "GEO" = "GE",

    "RUSSIA" = "RU",
    "RUSSIAN FEDERATION" = "RU",
    "RUS" = "RU"
  )

  unresolved <- is.na(out)

  if (any(unresolved)) {

    z <- y[unresolved]

    matched <- unname(
      map[
        match(
          z,
          names(map)
        )
      ]
    )

    out[unresolved] <- matched
  }

  out
}

region_group <- function(code) {

  dplyr::case_when(
    code == "GB" ~ "United Kingdom",
    code == "FR" ~ "France",
    code == "IT" ~ "Italy",
    code %in% c(
      "ES", "PT", "SI", "HR", "BA",
      "ME", "AL", "GR", "CY", "MT"
    ) ~ "Mediterranean",
    code %in% c("RO", "BG", "TR", "UA", "GE", "RU") ~ "Black Sea",
    TRUE ~ "Other Europe"
  )
}

record_id_col <- function(df) {

  pick_col(
    df,
    c(
      "gbifID",
      "gbif_id",
      "key",
      "occurrenceID",
      "occurrence_id"
    ),
    required = FALSE
  )
}

id_flagged <- record_id_col(flagged)
id_clean   <- record_id_col(clean)
id_coastal <- record_id_col(coastal)
id_main    <- record_id_col(main_occ)
id_occ     <- record_id_col(occ)

has_id_trace <- !any(
  is.na(
    c(
      id_flagged,
      id_clean,
      id_coastal,
      id_main,
      id_occ
    )
  )
)

# ------------------------------------------------------------
# 4. Restrict to Sara/Ribotta cells absent from final REV1
# ------------------------------------------------------------

missing_cells <- cell_audit %>%
  filter(
    sara_only_cell
  )

missing_points <- point_audit %>%
  filter(
    within_valid_coastal_domain,
    domain_cell_id %in% missing_cells$domain_cell_id
  )

# Representative point per missing cell:
# the Sara point closest to a flagged GBIF record.

representatives <- missing_points %>%
  group_by(domain_cell_id) %>%
  arrange(
    flagged_distance_m,
    .by_group = TRUE
  ) %>%
  slice(1) %>%
  ungroup()

# ------------------------------------------------------------
# 5. Attach nearest flagged GBIF record metadata
# ------------------------------------------------------------

trace_rows <- vector(
  "list",
  nrow(representatives)
)

for (i in seq_len(nrow(representatives))) {

  r <- representatives[i, ]

  idx <- r$flagged_nearest_index

  gbif_near <- is.finite(r$flagged_distance_m) &&
    r$flagged_distance_m <= 1000 &&
    !is.na(idx) &&
    idx >= 1 &&
    idx <= nrow(flagged)

  if (gbif_near) {

    rec_id <- if (!is.na(id_flagged)) {
      as.character(flagged[[id_flagged]][idx])
    } else {
      NA_character_
    }

    country_raw <- first_existing_value(
      flagged,
      idx,
      c(
        "countryCode",
        "country_code",
        "country",
        "countrycode"
      )
    )

    basis <- first_existing_value(
      flagged,
      idx,
      c(
        "basisOfRecord",
        "basis_of_record",
        "basis"
      )
    )

    uncertainty <- suppressWarnings(
      as.numeric(
        first_existing_value(
          flagged,
          idx,
          c(
            "coordinateUncertaintyInMeters",
            "coordinate_uncertainty",
            "uncertainty"
          )
        )
      )
    )

    cultivated <- first_existing_value(
      flagged,
      idx,
      c(
        "likely_cultivated",
        "cultivated",
        "is_cultivated"
      )
    )

    # Exact record tracing by GBIF/occurrence identifier when available.
    if (has_id_trace && !is.na(rec_id) && nzchar(rec_id)) {

      in_clean <- rec_id %in%
        as.character(clean[[id_clean]])

      in_coastal <- rec_id %in%
        as.character(coastal[[id_coastal]])

      in_main <- rec_id %in%
        as.character(main_occ[[id_main]])

      in_final <- rec_id %in%
        as.character(occ[[id_occ]])

    } else {

      in_clean <- NA
      in_coastal <- NA
      in_main <- NA
      in_final <- NA
    }

    country_code <- normalize_country(
      country_raw
    )

    # Identify likely cleaning reason.
    allowed_basis <- c(
      "HUMAN_OBSERVATION",
      "OBSERVATION",
      "PRESERVED_SPECIMEN",
      "HUMANOBSERVATION",
      "PRESERVEDSPECIMEN"
    )

    basis_norm <- toupper(
      gsub(
        "[^A-Z_]",
        "",
        basis
      )
    )

    basis_allowed <- if (
      is.na(basis) ||
      !nzchar(basis)
    ) {
      NA
    } else {
      basis_norm %in%
        allowed_basis
    }

    likely_cultivated <- toupper(
      as.character(cultivated)
    ) %in%
      c("TRUE", "T", "1", "YES")

    if (isTRUE(in_final)) {

      fate <- "nearest_GBIF_record_is_final_but_in_different_cell"

    } else if (identical(in_clean, FALSE)) {

      if (
        !is.na(basis_allowed) &&
        !basis_allowed
      ) {

        fate <- "excluded_before_old_clean_basisOfRecord"

      } else if (
        is.finite(uncertainty) &&
        uncertainty > 10000
      ) {

        fate <- "excluded_before_old_clean_uncertainty_gt_10km"

      } else if (
        likely_cultivated
      ) {

        fate <- "excluded_before_old_clean_likely_cultivated"

      } else {

        fate <- "excluded_before_old_clean_other_coordinate_or_QC_filter"
      }

    } else if (isTRUE(in_clean) && identical(in_coastal, FALSE)) {

      fate <- "passed_old_clean_but_not_retained_in_old_coastal_set"

    } else if (isTRUE(in_coastal) && identical(in_main, FALSE)) {

      fate <- "old_coastal_record_removed_in_REV1_correction"

    } else if (isTRUE(in_main) && identical(in_final, FALSE)) {

      fate <- "corrected_coastal_record_not_in_final_invaded_set"

    } else {

      fate <- "GBIF_within_1km_but_exact_record_trace_unresolved"
    }

  } else {

    rec_id <- NA_character_
    country_raw <- NA_character_
    country_code <- normalize_country(
      r$country_code
    )
    basis <- NA_character_
    uncertainty <- NA_real_
    cultivated <- NA_character_
    in_clean <- FALSE
    in_coastal <- FALSE
    in_main <- FALSE
    in_final <- FALSE
    fate <- "no_GBIF_flagged_record_within_1km"
  }

  trace_rows[[i]] <- tibble(
    domain_cell_id = r$domain_cell_id,

    sara_id = r$sara_id,

    sara_lon = r$sara_lon,

    sara_lat = r$sara_lat,

    n_sara_points_in_cell =
      missing_cells$n_sara_points[
        match(
          r$domain_cell_id,
          missing_cells$domain_cell_id
        )
      ],

    nearest_flagged_distance_m =
      r$flagged_distance_m,

    gbif_record_id =
      rec_id,

    country_raw =
      country_raw,

    country_code =
      country_code,

    region_group =
      region_group(
        country_code
      ),

    basisOfRecord =
      basis,

    coordinate_uncertainty_m =
      uncertainty,

    likely_cultivated =
      cultivated,

    exact_record_in_old_clean =
      in_clean,

    exact_record_in_old_coastal =
      in_coastal,

    exact_record_in_corrected_allrange =
      in_main,

    exact_record_in_final_invaded =
      in_final,

    fate =
      fate
  )
}

trace_tbl <- bind_rows(
  trace_rows
)

# If country is still unresolved, retain the country code already assigned
# during the reconciliation script as a fallback.

fallback_code <- normalize_country(
  representatives$country_code
)

use_fallback <- (
  is.na(trace_tbl$country_code) |
    trace_tbl$country_code == ""
) &
  !is.na(fallback_code)

trace_tbl$country_code[
  use_fallback
] <- fallback_code[
  use_fallback
]

trace_tbl$region_group <-
  region_group(
    trace_tbl$country_code
  )

# ------------------------------------------------------------
# 6. Summaries
# ------------------------------------------------------------

fate_summary <- trace_tbl %>%
  count(
    fate,
    sort = TRUE,
    name = "n_cells"
  )

country_summary <- trace_tbl %>%
  mutate(
    country_code = ifelse(
      is.na(country_code) |
        country_code == "",
      "UNRESOLVED",
      country_code
    )
  ) %>%
  group_by(country_code) %>%
  summarise(
    missing_sara_cells = n(),

    cells_with_GBIF_within_1km =
      sum(
        nearest_flagged_distance_m <= 1000,
        na.rm = TRUE
      ),

    cells_without_GBIF_within_1km =
      sum(
        !is.finite(nearest_flagged_distance_m) |
          nearest_flagged_distance_m > 1000,
        na.rm = TRUE
      ),

    old_clean_absences =
      sum(
        exact_record_in_old_clean %in% FALSE,
        na.rm = TRUE
      ),

    old_coastal_absences =
      sum(
        exact_record_in_old_clean %in% TRUE &
          exact_record_in_old_coastal %in% FALSE,
        na.rm = TRUE
      ),

    .groups = "drop"
  ) %>%
  arrange(
    desc(
      missing_sara_cells
    )
  )

reviewer_codes <- c(
  "GB", "FR", "IT",
  "ES", "PT",
  "GR", "HR", "SI", "ME", "AL",
  "TR", "BG", "RO",
  "CY", "MT"
)

reviewer_summary <- country_summary %>%
  filter(
    country_code %in%
      reviewer_codes
  )

independent_cells <- trace_tbl %>%
  filter(
    fate ==
      "no_GBIF_flagged_record_within_1km"
  ) %>%
  arrange(
    country_code,
    sara_lat,
    sara_lon
  )

# ------------------------------------------------------------
# 7. Additional diagnostic of old-clean exclusion fields
# ------------------------------------------------------------

qc_candidates <- grep(
  paste(
    c(
      "^cc_",
      "centroid",
      "institution",
      "capital",
      "zero",
      "equal",
      "sea",
      "clean",
      "flag"
    ),
    collapse = "|"
  ),
  names(flagged),
  ignore.case = TRUE,
  value = TRUE
)

qc_detail <- NULL

if (
  length(qc_candidates) > 0 &&
  any(
    trace_tbl$fate ==
      "excluded_before_old_clean_other_coordinate_or_QC_filter"
  )
) {

  unresolved_trace <- trace_tbl %>%
    filter(
      fate ==
        "excluded_before_old_clean_other_coordinate_or_QC_filter",
      !is.na(gbif_record_id)
    )

  if (
    nrow(unresolved_trace) > 0 &&
    !is.na(id_flagged)
  ) {

    idx <- match(
      unresolved_trace$gbif_record_id,
      as.character(
        flagged[[id_flagged]]
      )
    )

    good <- !is.na(idx)

    if (any(good)) {

      qc_detail <- bind_cols(
        unresolved_trace[good, ],
        flagged[
          idx[good],
          qc_candidates,
          drop = FALSE
        ]
      )
    }
  }
}

# ------------------------------------------------------------
# 8. Save outputs
# ------------------------------------------------------------

write_csv(
  trace_tbl,
  file.path(
    audit_dir,
    "04A_missing_Ribotta_cells_trace.csv"
  )
)

write_csv(
  fate_summary,
  file.path(
    audit_dir,
    "04B_missing_cell_fate_summary.csv"
  )
)

write_csv(
  country_summary,
  file.path(
    audit_dir,
    "04C_missing_cells_by_country.csv"
  )
)

write_csv(
  reviewer_summary,
  file.path(
    audit_dir,
    "04D_reviewer_countries_missing_cells.csv"
  )
)

write_csv(
  independent_cells,
  file.path(
    audit_dir,
    "04E_independent_Sara_Ribotta_cells.csv"
  )
)

if (!is.null(qc_detail)) {

  write_csv(
    qc_detail,
    file.path(
      audit_dir,
      "04F_unresolved_old_clean_QC_details.csv"
    )
  )
}

saveRDS(
  list(
    trace = trace_tbl,
    fate_summary = fate_summary,
    country_summary = country_summary,
    reviewer_summary = reviewer_summary,
    independent_cells = independent_cells,
    qc_fields = qc_candidates
  ),
  file.path(
    work_dir,
    "04_missing_Ribotta_cells_trace.rds"
  )
)

# ------------------------------------------------------------
# 9. Console report
# ------------------------------------------------------------

cat(
  "\n============================================================\n"
)

cat(
  "REV2 - TRACE OF MISSING RIBOTTA CELLS COMPLETE\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "Sara/Ribotta coastal cells absent from final REV1:",
  nrow(trace_tbl),
  "\n\n"
)

cat(
  "Fate summary:\n"
)

print(
  fate_summary,
  n = nrow(fate_summary)
)

cat(
  "\nReviewer-focused country summary:\n"
)

print(
  reviewer_summary,
  n = nrow(reviewer_summary)
)

cat(
  "\nPotentially independent Sara/Ribotta cells (no flagged GBIF record within 1 km):\n"
)

print(
  independent_cells,
  n = nrow(independent_cells)
)

cat(
  "\nIdentifier tracing available:",
  has_id_trace,
  "\n"
)

cat(
  "\nQC fields detected in flagged object:\n"
)

print(
  qc_candidates
)

cat(
  "\nKey outputs:\n"
)

cat(
  file.path(
    audit_dir,
    "04B_missing_cell_fate_summary.csv"
  ),
  "\n"
)

cat(
  file.path(
    audit_dir,
    "04D_reviewer_countries_missing_cells.csv"
  ),
  "\n"
)

cat(
  file.path(
    audit_dir,
    "04E_independent_Sara_Ribotta_cells.csv"
  ),
  "\n"
)

cat(
  "\nNo occurrence has been added. No model has been fitted.\n"
)

cat(
  "============================================================\n"
)
