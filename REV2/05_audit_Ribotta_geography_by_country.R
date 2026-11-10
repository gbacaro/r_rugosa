# ============================================================
# ROSA RUGOSA - REVISION 2
# 05_audit_Ribotta_geography_by_country.R
#
# Purpose
# -------
# Audit ALL Sara/Ribotta points by country, including points that
# fall outside the current raster-valid coastal domain.
#
# This script answers the reviewer-critical questions:
# - Are Italian / French / UK / Mediterranean Ribotta points present?
# - Are they represented in the REV1 final occurrence set?
# - If not, are they excluded because they fall outside the current
#   raster-valid coastal domain?
# - Are any likely independent field records absent from GBIF?
#
# No occurrence is added.
# No model is fitted.
# No REV1 file is modified.
# ============================================================

rm(list = ls())
gc()

required_pkgs <- c("dplyr", "readr", "tibble")
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

rec_file <- file.path(
  work_dir,
  "03_Ribotta_REV1_reconciliation.rds"
)

if (!file.exists(rec_file)) {
  stop(
    "Missing reconciliation object. Run 03_reconcile_Ribotta_with_REV1.R first."
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
# 2. Load reconciliation + REV1
# ------------------------------------------------------------

rec <- readRDS(rec_file)
point_audit <- rec$point_audit

rev1 <- new.env(parent = baseenv())
load(rdata_revision, envir = rev1)

if (!exists("flagged", envir = rev1, inherits = FALSE)) {
  stop("REV1 object 'flagged' not found.")
}

flagged <- get("flagged", envir = rev1)

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

normalize_country <- function(x) {

  y <- toupper(trimws(as.character(x)))
  out <- rep(NA_character_, length(y))

  iso2 <- grepl("^[A-Z]{2}$", y)
  out[iso2] <- y[iso2]

  map <- c(
    "UNITED KINGDOM" = "GB", "GREAT BRITAIN" = "GB", "GBR" = "GB", "UK" = "GB",
    "FRANCE" = "FR", "FRA" = "FR",
    "ITALY" = "IT", "ITALIA" = "IT", "ITA" = "IT",
    "SPAIN" = "ES", "ESPANA" = "ES", "ESPAÑA" = "ES", "ESP" = "ES",
    "PORTUGAL" = "PT", "PRT" = "PT",
    "GREECE" = "GR", "GRC" = "GR",
    "CROATIA" = "HR", "HRV" = "HR",
    "SLOVENIA" = "SI", "SVN" = "SI",
    "MONTENEGRO" = "ME", "MNE" = "ME",
    "ALBANIA" = "AL", "ALB" = "AL",
    "BOSNIA AND HERZEGOVINA" = "BA", "BIH" = "BA",
    "TURKEY" = "TR", "TÜRKIYE" = "TR", "TURKIYE" = "TR", "TUR" = "TR",
    "BULGARIA" = "BG", "BGR" = "BG",
    "ROMANIA" = "RO", "ROU" = "RO",
    "CYPRUS" = "CY", "CYP" = "CY",
    "MALTA" = "MT", "MLT" = "MT",
    "GERMANY" = "DE", "DEU" = "DE",
    "DENMARK" = "DK", "DNK" = "DK",
    "NORWAY" = "NO", "NOR" = "NO",
    "SWEDEN" = "SE", "SWE" = "SE",
    "FINLAND" = "FI", "FIN" = "FI",
    "POLAND" = "PL", "POL" = "PL",
    "NETHERLANDS" = "NL", "NLD" = "NL",
    "BELGIUM" = "BE", "BEL" = "BE",
    "IRELAND" = "IE", "IRL" = "IE",
    "ESTONIA" = "EE", "EST" = "EE",
    "LATVIA" = "LV", "LVA" = "LV",
    "LITHUANIA" = "LT", "LTU" = "LT",
    "ICELAND" = "IS", "ISL" = "IS",
    "UKRAINE" = "UA", "UKR" = "UA",
    "GEORGIA" = "GE", "GEO" = "GE",
    "RUSSIA" = "RU", "RUSSIAN FEDERATION" = "RU", "RUS" = "RU"
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

country_name <- function(code) {

  dplyr::recode(
    code,
    GB = "United Kingdom",
    FR = "France",
    IT = "Italy",
    ES = "Spain",
    PT = "Portugal",
    GR = "Greece",
    HR = "Croatia",
    SI = "Slovenia",
    ME = "Montenegro",
    AL = "Albania",
    BA = "Bosnia and Herzegovina",
    TR = "Türkiye",
    BG = "Bulgaria",
    RO = "Romania",
    CY = "Cyprus",
    MT = "Malta",
    DE = "Germany",
    DK = "Denmark",
    NO = "Norway",
    SE = "Sweden",
    FI = "Finland",
    PL = "Poland",
    NL = "Netherlands",
    BE = "Belgium",
    IE = "Ireland",
    EE = "Estonia",
    LV = "Latvia",
    LT = "Lithuania",
    IS = "Iceland",
    UA = "Ukraine",
    GE = "Georgia",
    RU = "Russia",
    .default = code
  )
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

# ------------------------------------------------------------
# 4. Identify GBIF metadata columns
# ------------------------------------------------------------

flag_country_col <- pick_col(
  flagged,
  c(
    "countryCode",
    "country_code",
    "country",
    "countrycode"
  ),
  required = FALSE
)

flag_id_col <- pick_col(
  flagged,
  c(
    "gbifID",
    "gbif_id",
    "key",
    "occurrenceID",
    "occurrence_id"
  ),
  required = FALSE
)

flag_basis_col <- pick_col(
  flagged,
  c(
    "basisOfRecord",
    "basis_of_record",
    "basis"
  ),
  required = FALSE
)

flag_unc_col <- pick_col(
  flagged,
  c(
    "coordinateUncertaintyInMeters",
    "coordinate_uncertainty",
    "uncertainty"
  ),
  required = FALSE
)

# ------------------------------------------------------------
# 5. Assign country to every Sara point using nearest GBIF record
# ------------------------------------------------------------

pa <- point_audit

pa$country_from_nearest_gbif <- NA_character_
pa$nearest_gbif_id <- NA_character_
pa$nearest_gbif_basis <- NA_character_
pa$nearest_gbif_uncertainty_m <- NA_real_

valid_idx <- !is.na(pa$flagged_nearest_index) &
  pa$flagged_nearest_index >= 1 &
  pa$flagged_nearest_index <= nrow(flagged)

if (!is.na(flag_country_col)) {

  raw_country <- rep(NA_character_, nrow(pa))

  raw_country[valid_idx] <-
    as.character(
      flagged[[flag_country_col]][
        pa$flagged_nearest_index[valid_idx]
      ]
    )

  pa$country_from_nearest_gbif <-
    normalize_country(
      raw_country
    )
}

if (!is.na(flag_id_col)) {

  pa$nearest_gbif_id[valid_idx] <-
    as.character(
      flagged[[flag_id_col]][
        pa$flagged_nearest_index[valid_idx]
      ]
    )
}

if (!is.na(flag_basis_col)) {

  pa$nearest_gbif_basis[valid_idx] <-
    as.character(
      flagged[[flag_basis_col]][
        pa$flagged_nearest_index[valid_idx]
      ]
    )
}

if (!is.na(flag_unc_col)) {

  pa$nearest_gbif_uncertainty_m[valid_idx] <-
    suppressWarnings(
      as.numeric(
        flagged[[flag_unc_col]][
          pa$flagged_nearest_index[valid_idx]
        ]
      )
    )
}

# Prefer the nearest GBIF country only when it is reasonably close.
# First use <=1 km, then <=5 km. Otherwise retain the country assignment
# inherited from the previous reconciliation if available.

existing_code <- normalize_country(
  pa$country_code
)

pa$country_code_final <- existing_code

use_1km <- is.finite(pa$flagged_distance_m) &
  pa$flagged_distance_m <= 1000 &
  !is.na(pa$country_from_nearest_gbif)

pa$country_code_final[use_1km] <-
  pa$country_from_nearest_gbif[use_1km]

use_5km <- (
  is.na(pa$country_code_final) |
    pa$country_code_final == ""
) &
  is.finite(pa$flagged_distance_m) &
  pa$flagged_distance_m <= 5000 &
  !is.na(pa$country_from_nearest_gbif)

pa$country_code_final[use_5km] <-
  pa$country_from_nearest_gbif[use_5km]

pa$country_name_final <-
  country_name(
    pa$country_code_final
  )

pa$region_group_final <-
  region_group(
    pa$country_code_final
  )

# ------------------------------------------------------------
# 6. Status variables
# ------------------------------------------------------------

pa$inside_domain_1p5km <-
  is.finite(pa$domain_distance_m) &
  pa$domain_distance_m <= 1500

pa$outside_domain_1p5km <-
  !pa$inside_domain_1p5km

pa$represented_final_within_1km <-
  is.finite(pa$final_distance_m) &
  pa$final_distance_m <= 1000

pa$represented_final_within_5km <-
  is.finite(pa$final_distance_m) &
  pa$final_distance_m <= 5000

pa$likely_independent_from_GBIF <-
  !is.finite(pa$flagged_distance_m) |
  pa$flagged_distance_m > 1000

# ------------------------------------------------------------
# 7. Country-level summary of ALL Sara points
# ------------------------------------------------------------

country_summary_all <- pa %>%
  mutate(
    country_code_final = ifelse(
      is.na(country_code_final) |
        country_code_final == "",
      "UNRESOLVED",
      country_code_final
    )
  ) %>%
  group_by(
    country_code_final,
    country_name_final
  ) %>%
  summarise(
    n_sara_points = n(),

    n_points_GBIF_within_100m =
      sum(
        flagged_distance_m <= 100,
        na.rm = TRUE
      ),

    n_points_GBIF_within_1km =
      sum(
        flagged_distance_m <= 1000,
        na.rm = TRUE
      ),

    n_points_GBIF_within_5km =
      sum(
        flagged_distance_m <= 5000,
        na.rm = TRUE
      ),

    n_points_inside_current_domain =
      sum(
        inside_domain_1p5km,
        na.rm = TRUE
      ),

    n_points_outside_current_domain =
      sum(
        outside_domain_1p5km,
        na.rm = TRUE
      ),

    n_points_final_within_1km =
      sum(
        represented_final_within_1km,
        na.rm = TRUE
      ),

    median_distance_to_domain_m =
      median(
        domain_distance_m,
        na.rm = TRUE
      ),

    max_distance_to_domain_m =
      max(
        domain_distance_m,
        na.rm = TRUE
      ),

    .groups = "drop"
  ) %>%
  arrange(
    desc(
      n_sara_points
    )
  )

country_summary_all$median_distance_to_domain_m[
  !is.finite(
    country_summary_all$median_distance_to_domain_m
  )
] <- NA_real_

country_summary_all$max_distance_to_domain_m[
  !is.finite(
    country_summary_all$max_distance_to_domain_m
  )
] <- NA_real_

# ------------------------------------------------------------
# 8. Cell-level summary by country
# ------------------------------------------------------------

valid_cells_by_country <- pa %>%
  filter(
    inside_domain_1p5km,
    !is.na(domain_cell_id)
  ) %>%
  mutate(
    country_code_final = ifelse(
      is.na(country_code_final) |
        country_code_final == "",
      "UNRESOLVED",
      country_code_final
    )
  ) %>%
  group_by(
    country_code_final,
    country_name_final,
    domain_cell_id
  ) %>%
  summarise(
    n_sara_points = n(),

    represented_in_final_REV1 =
      any(
        cell_represented_in_final_REV1,
        na.rm = TRUE
      ),

    min_GBIF_distance_m =
      min(
        flagged_distance_m,
        na.rm = TRUE
      ),

    .groups = "drop"
  )

valid_cells_by_country$min_GBIF_distance_m[
  !is.finite(
    valid_cells_by_country$min_GBIF_distance_m
  )
] <- NA_real_

cell_country_summary <- valid_cells_by_country %>%
  group_by(
    country_code_final,
    country_name_final
  ) %>%
  summarise(
    n_unique_valid_coastal_cells =
      n(),

    n_cells_already_final_REV1 =
      sum(
        represented_in_final_REV1,
        na.rm = TRUE
      ),

    n_cells_absent_final_REV1 =
      sum(
        !represented_in_final_REV1,
        na.rm = TRUE
      ),

    n_cells_without_GBIF_within_1km =
      sum(
        is.na(min_GBIF_distance_m) |
          min_GBIF_distance_m > 1000,
        na.rm = TRUE
      ),

    .groups = "drop"
  ) %>%
  arrange(
    desc(
      n_unique_valid_coastal_cells
    )
  )

# ------------------------------------------------------------
# 9. Reviewer-focused countries
# ------------------------------------------------------------

reviewer_codes <- c(
  "GB", "FR", "IT",
  "ES", "PT",
  "SI", "HR", "BA", "ME", "AL", "GR",
  "TR", "BG", "RO",
  "CY", "MT"
)

reviewer_point_summary <- country_summary_all %>%
  filter(
    country_code_final %in%
      reviewer_codes
  )

reviewer_cell_summary <- cell_country_summary %>%
  filter(
    country_code_final %in%
      reviewer_codes
  )

reviewer_points <- pa %>%
  filter(
    country_code_final %in%
      reviewer_codes
  ) %>%
  arrange(
    country_code_final,
    sara_lat,
    sara_lon
  )

# ------------------------------------------------------------
# 10. Italy / France / UK detail
# ------------------------------------------------------------

critical_points <- pa %>%
  filter(
    country_code_final %in%
      c("IT", "FR", "GB")
  ) %>%
  select(
    sara_id,
    sara_lon,
    sara_lat,
    country_code_final,
    country_name_final,
    flagged_distance_m,
    nearest_gbif_id,
    nearest_gbif_basis,
    nearest_gbif_uncertainty_m,
    domain_distance_m,
    inside_domain_1p5km,
    domain_cell_id,
    final_distance_m,
    represented_final_within_1km,
    represented_final_within_5km,
    cell_represented_in_final_REV1,
    likely_independent_from_GBIF,
    everything()
  ) %>%
  arrange(
    country_code_final,
    inside_domain_1p5km,
    domain_distance_m
  )

# ------------------------------------------------------------
# 11. Candidate field observations absent from GBIF
# ------------------------------------------------------------

independent_candidates_all <- pa %>%
  filter(
    likely_independent_from_GBIF
  ) %>%
  arrange(
    country_code_final,
    flagged_distance_m
  )

independent_candidates_review <-
  independent_candidates_all %>%
  filter(
    country_code_final %in%
      reviewer_codes
  )

# ------------------------------------------------------------
# 12. Search Sara attributes for possible site/source labels
# ------------------------------------------------------------

ignore_fields <- c(
  "sara_id",
  "sara_lon",
  "sara_lat",
  "country_code",
  "country_name",
  "country_from_nearest_gbif",
  "country_code_final",
  "country_name_final",
  "region_group",
  "region_group_final",
  "nearest_gbif_id",
  "nearest_gbif_basis"
)

char_fields <- names(pa)[
  vapply(
    pa,
    function(x) is.character(x) || is.factor(x),
    logical(1)
  )
]

char_fields <- setdiff(
  char_fields,
  ignore_fields
)

site_keyword_hits <- list()

keywords <- c(
  "brussa",
  "caorle",
  "venez",
  "ital",
  "ribotta",
  "adriatic",
  "mediterr"
)

for (cl in char_fields) {

  vals <- as.character(
    pa[[cl]]
  )

  hit <- Reduce(
    `|`,
    lapply(
      keywords,
      function(k) {
        grepl(
          k,
          vals,
          ignore.case = TRUE
        )
      }
    )
  )

  hit[is.na(hit)] <- FALSE

  if (any(hit)) {

    site_keyword_hits[[cl]] <-
      tibble(
        field = cl,
        sara_id = pa$sara_id[hit],
        value = vals[hit],
        sara_lon = pa$sara_lon[hit],
        sara_lat = pa$sara_lat[hit],
        country_code_final =
          pa$country_code_final[hit],
        domain_distance_m =
          pa$domain_distance_m[hit],
        flagged_distance_m =
          pa$flagged_distance_m[hit]
      )
  }
}

if (length(site_keyword_hits) > 0) {

  site_keyword_hits <- bind_rows(
    site_keyword_hits
  )

} else {

  site_keyword_hits <- tibble(
    field = character(),
    sara_id = integer(),
    value = character(),
    sara_lon = double(),
    sara_lat = double(),
    country_code_final = character(),
    domain_distance_m = double(),
    flagged_distance_m = double()
  )
}

# ------------------------------------------------------------
# 13. Save outputs
# ------------------------------------------------------------

write_csv(
  country_summary_all,
  file.path(
    audit_dir,
    "05A_all_Sara_points_by_country.csv"
  )
)

write_csv(
  cell_country_summary,
  file.path(
    audit_dir,
    "05B_valid_coastal_cells_by_country.csv"
  )
)

write_csv(
  reviewer_point_summary,
  file.path(
    audit_dir,
    "05C_reviewer_country_point_summary.csv"
  )
)

write_csv(
  reviewer_cell_summary,
  file.path(
    audit_dir,
    "05D_reviewer_country_cell_summary.csv"
  )
)

write_csv(
  reviewer_points,
  file.path(
    audit_dir,
    "05E_all_reviewer_country_points.csv"
  )
)

write_csv(
  critical_points,
  file.path(
    audit_dir,
    "05F_Italy_France_UK_points.csv"
  )
)

write_csv(
  independent_candidates_all,
  file.path(
    audit_dir,
    "05G_all_candidate_nonGBIF_points.csv"
  )
)

write_csv(
  independent_candidates_review,
  file.path(
    audit_dir,
    "05H_reviewer_country_candidate_nonGBIF_points.csv"
  )
)

write_csv(
  site_keyword_hits,
  file.path(
    audit_dir,
    "05I_Sara_attribute_site_keyword_hits.csv"
  )
)

saveRDS(
  list(
    point_audit = pa,
    country_summary_all = country_summary_all,
    cell_country_summary = cell_country_summary,
    reviewer_point_summary = reviewer_point_summary,
    reviewer_cell_summary = reviewer_cell_summary,
    critical_points = critical_points,
    independent_candidates = independent_candidates_all,
    site_keyword_hits = site_keyword_hits
  ),
  file.path(
    work_dir,
    "05_Ribotta_geography_audit.rds"
  )
)

# ------------------------------------------------------------
# 14. Console report
# ------------------------------------------------------------

cat(
  "\n============================================================\n"
)

cat(
  "REV2 - RIBOTTA GEOGRAPHIC AUDIT COMPLETE\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "Reviewer-focused POINT summary:\n"
)

print(
  reviewer_point_summary,
  n = nrow(reviewer_point_summary)
)

cat(
  "\nReviewer-focused CELL summary:\n"
)

print(
  reviewer_cell_summary,
  n = nrow(reviewer_cell_summary)
)

cat(
  "\nCandidate non-GBIF points in reviewer-focused countries:\n"
)

print(
  independent_candidates_review,
  n = nrow(independent_candidates_review)
)

cat(
  "\nSara attribute keyword hits (Brussa/Caorle/Italy/etc.):\n"
)

print(
  site_keyword_hits,
  n = min(
    nrow(site_keyword_hits),
    50
  )
)

cat(
  "\nKey files:\n"
)

cat(
  file.path(
    audit_dir,
    "05C_reviewer_country_point_summary.csv"
  ),
  "\n"
)

cat(
  file.path(
    audit_dir,
    "05D_reviewer_country_cell_summary.csv"
  ),
  "\n"
)

cat(
  file.path(
    audit_dir,
    "05F_Italy_France_UK_points.csv"
  ),
  "\n"
)

cat(
  file.path(
    audit_dir,
    "05H_reviewer_country_candidate_nonGBIF_points.csv"
  ),
  "\n"
)

cat(
  "\nNo occurrence has been added. No model has been fitted.\n"
)

cat(
  "============================================================\n"
)
