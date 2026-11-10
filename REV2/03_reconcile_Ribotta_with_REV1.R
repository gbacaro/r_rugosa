# ============================================================
# ROSA RUGOSA - REVISION 2
# 03_reconcile_Ribotta_with_REV1.R
#
# Reconcile Sara/Ribotta (2021) occurrence points with:
# - REV1 raw/clean/coastal occurrence stages
# - final invaded coastal occurrences
# - fixed raster-valid European coastal domain
#
# No occurrence is added to the model dataset.
# No model is refitted.
# No REV1 file is modified.
# ============================================================

rm(list = ls())
gc()

required_pkgs <- c("sf", "dplyr", "readr", "tibble")
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing packages: ",
    paste(missing_pkgs, collapse = ", "),
    ". Install them before running this script."
  )
}

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

if (!exists("input_dir")) input_dir <- file.path(rev2_dir, "01_input")
if (!exists("audit_dir")) audit_dir <- file.path(rev2_dir, "03_audit")
if (!exists("work_dir"))  work_dir  <- file.path(rev2_dir, "04_work")

dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

rdata_revision <- file.path(
  rev1_dir,
  "dati_per_revisione.RData"
)

if (!file.exists(rdata_revision)) {
  stop("Cannot find REV1 RData: ", rdata_revision)
}

# ------------------------------------------------------------
# 2. Load REV1 objects in a private environment
# ------------------------------------------------------------

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
eu_df    <- get_rev1("eu_df")

world_eu <- get_rev1("world_eu", required = FALSE)
if (is.null(world_eu)) {
  world_eu <- get_rev1("world", required = FALSE)
}

# ------------------------------------------------------------
# 3. Helper functions
# ------------------------------------------------------------

pick_col <- function(df, candidates, required = TRUE) {

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
      "Could not identify required column. Candidates: ",
      paste(candidates, collapse = ", "),
      "\nAvailable fields: ",
      paste(nms, collapse = ", ")
    )
  }

  NA_character_
}

make_sf_points <- function(df, label) {

  lon_col <- pick_col(
    df,
    c(
      "decimal_longitude",
      "decimallongitude",
      "longitude",
      "lon",
      "x"
    )
  )

  lat_col <- pick_col(
    df,
    c(
      "decimal_latitude",
      "decimallatitude",
      "latitude",
      "lat",
      "y"
    )
  )

  lon <- suppressWarnings(as.numeric(df[[lon_col]]))
  lat <- suppressWarnings(as.numeric(df[[lat_col]]))

  ok <- is.finite(lon) &
    is.finite(lat) &
    abs(lon) <= 180 &
    abs(lat) <= 90

  if (!any(ok)) {
    stop(label, ": no valid longitude/latitude values found.")
  }

  out <- df[ok, , drop = FALSE]
  out$.source_row <- which(ok)
  out$.lon <- lon[ok]
  out$.lat <- lat[ok]

  sf::st_as_sf(
    out,
    coords = c(".lon", ".lat"),
    crs = 4326,
    remove = FALSE
  )
}

nearest_distance <- function(query_sf, target_sf, prefix) {

  q <- sf::st_transform(query_sf, 3035)
  t <- sf::st_transform(target_sf, 3035)

  idx <- sf::st_nearest_feature(q, t)

  dist <- as.numeric(
    sf::st_distance(
      q,
      t[idx, ],
      by_element = TRUE
    )
  )

  tibble(
    !!paste0(prefix, "_nearest_index") := idx,
    !!paste0(prefix, "_distance_m") := dist
  )
}

normalize_country_code <- function(x) {

  y <- toupper(trimws(as.character(x)))

  y[y %in% c("UNITED KINGDOM", "GREAT BRITAIN", "UK")] <- "GB"
  y[y %in% c("FRANCE", "FRENCH REPUBLIC")] <- "FR"
  y[y %in% c("ITALY", "ITALIA")] <- "IT"
  y[y %in% c("SPAIN", "ESPAÑA", "ESPANA")] <- "ES"
  y[y == "PORTUGAL"] <- "PT"
  y[y == "GREECE"] <- "GR"
  y[y == "CROATIA"] <- "HR"
  y[y == "SLOVENIA"] <- "SI"
  y[y == "MONTENEGRO"] <- "ME"
  y[y == "ALBANIA"] <- "AL"
  y[y %in% c("TURKEY", "TÜRKIYE", "TURKIYE")] <- "TR"
  y[y == "BULGARIA"] <- "BG"
  y[y == "ROMANIA"] <- "RO"
  y[y == "CYPRUS"] <- "CY"
  y[y == "MALTA"] <- "MT"

  y
}

region_group <- function(code) {

  dplyr::case_when(
    code == "GB" ~ "United Kingdom",
    code == "FR" ~ "France",
    code == "IT" ~ "Italy",
    code %in% c(
      "ES", "PT", "FR", "IT", "SI", "HR", "BA",
      "ME", "AL", "GR", "TR", "CY", "MT", "MC"
    ) ~ "Mediterranean",
    code %in% c("RO", "BG", "TR", "UA", "GE", "RU") ~ "Black Sea",
    TRUE ~ "Other Europe"
  )
}

# ------------------------------------------------------------
# 4. Locate the Sara/Ribotta shapefile
# ------------------------------------------------------------

shp_files <- list.files(
  rev2_dir,
  pattern = "\\.shp$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(shp_files) == 0) {

  zip_candidates <- list.files(
    rev2_dir,
    pattern = "\\.zip$",
    recursive = FALSE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (length(zip_candidates) == 0) {
    stop("No Sara/Ribotta shapefile or ZIP archive found.")
  }

  sara_extract_dir <- file.path(
    input_dir,
    "Sara_Ribotta2021"
  )

  dir.create(
    sara_extract_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  unzip(
    zip_candidates[1],
    exdir = sara_extract_dir,
    overwrite = TRUE
  )

  shp_files <- list.files(
    sara_extract_dir,
    pattern = "\\.shp$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
}

if (length(shp_files) == 0) {
  stop("No .shp file found after extraction.")
}

shp_inventory <- lapply(
  shp_files,
  function(f) {

    x <- try(
      sf::st_read(f, quiet = TRUE),
      silent = TRUE
    )

    if (inherits(x, "try-error")) return(NULL)

    tibble(
      file = f,
      n = nrow(x),
      geometry = paste(
        unique(as.character(sf::st_geometry_type(x))),
        collapse = ";"
      )
    )
  }
)

shp_inventory <- dplyr::bind_rows(shp_inventory)

if (nrow(shp_inventory) == 0) {
  stop("All shapefiles failed to read.")
}

sara_file <- shp_inventory$file[
  which.max(shp_inventory$n)
]

sara_raw <- sf::st_read(
  sara_file,
  quiet = TRUE
)

if (is.na(sf::st_crs(sara_raw))) {
  stop("Sara/Ribotta shapefile has no CRS.")
}

sara_raw <- sf::st_transform(
  sara_raw,
  4326
)

sara <- suppressWarnings(
  sf::st_cast(
    sara_raw,
    "POINT"
  )
)

sara$sara_id <- seq_len(nrow(sara))

sara_xy <- sf::st_coordinates(sara)
sara$sara_lon <- sara_xy[, "X"]
sara$sara_lat <- sara_xy[, "Y"]

# ------------------------------------------------------------
# 5. Spatial objects for REV1 stages
# ------------------------------------------------------------

flagged_sf <- make_sf_points(flagged, "flagged")
clean_sf   <- make_sf_points(clean, "clean")
coastal_sf <- make_sf_points(coastal, "coastal")
main_sf    <- make_sf_points(main_occ, "main_occ")
occ_sf     <- make_sf_points(occ, "occ")
eu_sf      <- make_sf_points(eu_df, "eu_df")

eu_sf$domain_cell_id <- seq_len(nrow(eu_sf))

# ------------------------------------------------------------
# 6. Country assignment
# ------------------------------------------------------------

sara$country_code <- NA_character_
sara$country_name <- NA_character_

if (
  !is.null(world_eu) &&
  inherits(world_eu, "sf")
) {

  world_poly <- sf::st_transform(
    world_eu,
    4326
  )

  world_att <- sf::st_drop_geometry(
    world_poly
  )

  code_col <- pick_col(
    world_att,
    c(
      "iso_a2",
      "countrycode",
      "country_code",
      "adm0_a3"
    ),
    required = FALSE
  )

  name_col <- pick_col(
    world_att,
    c(
      "admin",
      "name",
      "country",
      "sovereignt"
    ),
    required = FALSE
  )

  keep_cols <- c(code_col, name_col)
  keep_cols <- keep_cols[!is.na(keep_cols)]

  if (length(keep_cols) > 0) {

    joined <- suppressWarnings(
      sf::st_join(
        sara,
        world_poly[, keep_cols, drop = FALSE],
        left = TRUE,
        join = sf::st_intersects
      )
    )

    if (
      !is.na(code_col) &&
      code_col %in% names(joined)
    ) {
      sara$country_code <- normalize_country_code(
        joined[[code_col]]
      )
    }

    if (
      !is.na(name_col) &&
      name_col %in% names(joined)
    ) {
      sara$country_name <- as.character(
        joined[[name_col]]
      )
    }
  }
}

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

# ------------------------------------------------------------
# 7. Point-level reconciliation
# ------------------------------------------------------------

near_flagged <- nearest_distance(
  sara,
  flagged_sf,
  "flagged"
)

near_clean <- nearest_distance(
  sara,
  clean_sf,
  "clean"
)

near_coastal <- nearest_distance(
  sara,
  coastal_sf,
  "coastal"
)

near_main <- nearest_distance(
  sara,
  main_sf,
  "main"
)

near_final <- nearest_distance(
  sara,
  occ_sf,
  "final"
)

near_domain <- nearest_distance(
  sara,
  eu_sf,
  "domain"
)

point_audit <- dplyr::bind_cols(
  sf::st_drop_geometry(sara),
  near_flagged,
  near_clean,
  near_coastal,
  near_main,
  near_final,
  near_domain
)

if (!is.na(flag_country_col)) {

  borrow <- (
    is.na(point_audit$country_code) |
      point_audit$country_code == ""
  ) &
    is.finite(point_audit$flagged_distance_m) &
    point_audit$flagged_distance_m <= 1000

  if (any(borrow)) {

    nearest_country <- flagged[[flag_country_col]][
      point_audit$flagged_nearest_index[borrow]
    ]

    point_audit$country_code[borrow] <-
      normalize_country_code(
        nearest_country
      )
  }
}

point_audit$region_group <-
  region_group(
    point_audit$country_code
  )

# The REV1 rule allowed assignment to the nearest raster-valid coastal
# cell when the occurrence was within 1,500 m.

point_audit$within_valid_coastal_domain <-
  is.finite(point_audit$domain_distance_m) &
  point_audit$domain_distance_m <= 1500

point_audit$domain_cell_id <- ifelse(
  point_audit$within_valid_coastal_domain,
  eu_sf$domain_cell_id[
    point_audit$domain_nearest_index
  ],
  NA_integer_
)

point_audit$gbif_within_100m <-
  point_audit$flagged_distance_m <= 100

point_audit$gbif_within_500m <-
  point_audit$flagged_distance_m <= 500

point_audit$gbif_within_1km <-
  point_audit$flagged_distance_m <= 1000

point_audit$gbif_within_5km <-
  point_audit$flagged_distance_m <= 5000

point_audit$old_clean_within_1km <-
  point_audit$clean_distance_m <= 1000

point_audit$old_coastal_within_1km <-
  point_audit$coastal_distance_m <= 1000

point_audit$corrected_allrange_within_1km <-
  point_audit$main_distance_m <= 1000

point_audit$final_invaded_within_1km <-
  point_audit$final_distance_m <= 1000

point_audit$final_invaded_within_5km <-
  point_audit$final_distance_m <= 5000

point_audit$candidate_nonGBIF_point <-
  point_audit$within_valid_coastal_domain &
  (
    is.na(point_audit$gbif_within_1km) |
      !point_audit$gbif_within_1km
  )

# ------------------------------------------------------------
# 8. Reconstruct final REV1 invaded cells
# ------------------------------------------------------------

occ_to_domain <- nearest_distance(
  occ_sf,
  eu_sf,
  "domain"
)

occ_cells <- tibble(
  occ_row = seq_len(nrow(occ_sf)),
  domain_distance_m =
    occ_to_domain$domain_distance_m,
  domain_cell_id =
    eu_sf$domain_cell_id[
      occ_to_domain$domain_nearest_index
    ]
) %>%
  filter(
    is.finite(domain_distance_m),
    domain_distance_m <= 1500
  )

final_cells <- sort(
  unique(
    occ_cells$domain_cell_id
  )
)

point_audit$cell_represented_in_final_REV1 <-
  point_audit$domain_cell_id %in%
  final_cells

point_audit$sara_adds_new_cell_vs_REV1 <-
  point_audit$within_valid_coastal_domain &
  !point_audit$cell_represented_in_final_REV1

# ------------------------------------------------------------
# 9. Cell-level reconciliation
# ------------------------------------------------------------

sara_valid <- point_audit %>%
  filter(
    within_valid_coastal_domain,
    !is.na(domain_cell_id)
  )

cell_audit <- sara_valid %>%
  group_by(domain_cell_id) %>%
  summarise(
    n_sara_points = n(),

    min_gbif_distance_m =
      suppressWarnings(
        min(
          flagged_distance_m,
          na.rm = TRUE
        )
      ),

    any_gbif_within_100m =
      any(
        gbif_within_100m,
        na.rm = TRUE
      ),

    any_gbif_within_500m =
      any(
        gbif_within_500m,
        na.rm = TRUE
      ),

    any_gbif_within_1km =
      any(
        gbif_within_1km,
        na.rm = TRUE
      ),

    any_old_clean_within_1km =
      any(
        old_clean_within_1km,
        na.rm = TRUE
      ),

    any_old_coastal_within_1km =
      any(
        old_coastal_within_1km,
        na.rm = TRUE
      ),

    any_final_invaded_within_1km =
      any(
        final_invaded_within_1km,
        na.rm = TRUE
      ),

    country_code =
      dplyr::first(
        country_code[
          !is.na(country_code) &
            country_code != ""
        ],
        default = NA_character_
      ),

    region_group =
      dplyr::first(
        region_group[
          !is.na(region_group) &
            region_group != ""
        ],
        default = NA_character_
      ),

    .groups = "drop"
  ) %>%
  mutate(
    represented_in_final_REV1 =
      domain_cell_id %in% final_cells,

    sara_only_cell =
      !represented_in_final_REV1,

    candidate_independent_cell =
      sara_only_cell &
      !any_gbif_within_1km
  )

cell_audit$min_gbif_distance_m[
  !is.finite(
    cell_audit$min_gbif_distance_m
  )
] <- NA_real_

# ------------------------------------------------------------
# 10. Country and regional summaries
# ------------------------------------------------------------

country_summary <- cell_audit %>%
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
    sara_points =
      sum(n_sara_points),

    sara_unique_valid_coastal_cells =
      n(),

    cells_already_in_final_REV1 =
      sum(
        represented_in_final_REV1
      ),

    sara_only_cells =
      sum(
        sara_only_cell
      ),

    candidate_independent_cells =
      sum(
        candidate_independent_cell
      ),

    .groups = "drop"
  ) %>%
  arrange(
    desc(
      sara_unique_valid_coastal_cells
    )
  )

regional_summary <- cell_audit %>%
  mutate(
    region_group = ifelse(
      is.na(region_group) |
        region_group == "",
      "Unresolved",
      region_group
    )
  ) %>%
  group_by(region_group) %>%
  summarise(
    sara_points =
      sum(n_sara_points),

    sara_unique_valid_coastal_cells =
      n(),

    cells_already_in_final_REV1 =
      sum(
        represented_in_final_REV1
      ),

    sara_only_cells =
      sum(
        sara_only_cell
      ),

    candidate_independent_cells =
      sum(
        candidate_independent_cell
      ),

    .groups = "drop"
  ) %>%
  arrange(
    desc(
      sara_unique_valid_coastal_cells
    )
  )

reviewer_codes <- c(
  "GB", "FR", "IT",
  "ES", "PT",
  "GR", "HR", "SI", "ME", "AL",
  "TR", "BG", "RO",
  "CY", "MT"
)

reviewer_country_summary <-
  country_summary %>%
  filter(
    country_code %in%
      reviewer_codes
  )

# ------------------------------------------------------------
# 11. Candidate independent Sara/Ribotta information
# ------------------------------------------------------------

candidate_points <- point_audit %>%
  filter(
    within_valid_coastal_domain,
    sara_adds_new_cell_vs_REV1
  ) %>%
  arrange(
    flagged_distance_m
  )

candidate_independent_points <-
  point_audit %>%
  filter(
    within_valid_coastal_domain,
    sara_adds_new_cell_vs_REV1,
    candidate_nonGBIF_point
  ) %>%
  arrange(
    desc(
      flagged_distance_m
    )
  )

candidate_cells <- cell_audit %>%
  filter(
    sara_only_cell
  ) %>%
  arrange(
    desc(
      candidate_independent_cell
    ),
    desc(
      n_sara_points
    )
  )

# ------------------------------------------------------------
# 12. Global summary
# ------------------------------------------------------------

summary_tbl <- tibble(
  metric = c(
    "Sara source point geometries",

    "Sara points within 1.5 km of raster-valid coastal domain",

    "Sara unique raster-valid coastal cells",

    "Sara coastal cells represented in final REV1 invaded dataset",

    "Sara coastal cells absent from final REV1 invaded dataset",

    "Sara-only cells with a GBIF flagged record within 1 km",

    "Sara-only cells without a GBIF flagged record within 1 km",

    "Sara points with GBIF flagged record within 100 m",

    "Sara points with GBIF flagged record within 500 m",

    "Sara points with GBIF flagged record within 1 km",

    "Sara points with GBIF flagged record within 5 km",

    "Final REV1 invaded cells reconstructed",

    "Final REV1 invaded occurrence records"
  ),

  value = c(
    nrow(point_audit),

    sum(
      point_audit$within_valid_coastal_domain,
      na.rm = TRUE
    ),

    n_distinct(
      sara_valid$domain_cell_id
    ),

    sum(
      cell_audit$represented_in_final_REV1,
      na.rm = TRUE
    ),

    sum(
      cell_audit$sara_only_cell,
      na.rm = TRUE
    ),

    sum(
      cell_audit$sara_only_cell &
        cell_audit$any_gbif_within_1km,
      na.rm = TRUE
    ),

    sum(
      cell_audit$candidate_independent_cell,
      na.rm = TRUE
    ),

    sum(
      point_audit$gbif_within_100m,
      na.rm = TRUE
    ),

    sum(
      point_audit$gbif_within_500m,
      na.rm = TRUE
    ),

    sum(
      point_audit$gbif_within_1km,
      na.rm = TRUE
    ),

    sum(
      point_audit$gbif_within_5km,
      na.rm = TRUE
    ),

    length(
      final_cells
    ),

    nrow(
      occ
    )
  )
)

# ------------------------------------------------------------
# 13. Save outputs
# ------------------------------------------------------------

readr::write_csv(
  summary_tbl,
  file.path(
    audit_dir,
    "03A_reconciliation_summary.csv"
  )
)

readr::write_csv(
  point_audit,
  file.path(
    audit_dir,
    "03B_sara_point_reconciliation.csv"
  )
)

readr::write_csv(
  cell_audit,
  file.path(
    audit_dir,
    "03C_sara_cell_reconciliation.csv"
  )
)

readr::write_csv(
  country_summary,
  file.path(
    audit_dir,
    "03D_country_cell_summary.csv"
  )
)

readr::write_csv(
  regional_summary,
  file.path(
    audit_dir,
    "03E_regional_cell_summary.csv"
  )
)

readr::write_csv(
  reviewer_country_summary,
  file.path(
    audit_dir,
    "03F_reviewer_country_summary.csv"
  )
)

readr::write_csv(
  candidate_points,
  file.path(
    audit_dir,
    "03G_sara_only_cell_points.csv"
  )
)

readr::write_csv(
  candidate_independent_points,
  file.path(
    audit_dir,
    "03H_candidate_nonGBIF_sara_points.csv"
  )
)

readr::write_csv(
  candidate_cells,
  file.path(
    audit_dir,
    "03I_sara_only_cells.csv"
  )
)

saveRDS(
  list(
    summary = summary_tbl,
    point_audit = point_audit,
    cell_audit = cell_audit,
    country_summary = country_summary,
    regional_summary = regional_summary,
    final_cells = final_cells,
    sara_shapefile = sara_file
  ),
  file.path(
    work_dir,
    "03_Ribotta_REV1_reconciliation.rds"
  )
)

# ------------------------------------------------------------
# 14. Console report
# ------------------------------------------------------------

cat(
  "\n============================================================\n"
)

cat(
  "REV2 - RIBOTTA / REV1 RECONCILIATION COMPLETE\n"
)

cat(
  "============================================================\n\n"
)

print(
  summary_tbl,
  n = nrow(summary_tbl)
)

cat(
  "\nReviewer-focused country summary:\n"
)

print(
  reviewer_country_summary,
  n = nrow(reviewer_country_summary)
)

cat(
  "\nRegional summary:\n"
)

print(
  regional_summary,
  n = nrow(regional_summary)
)

cat(
  "\nSara shapefile used:\n",
  sara_file,
  "\n",
  sep = ""
)

cat(
  "\nKey outputs:\n"
)

cat(
  file.path(
    audit_dir,
    "03A_reconciliation_summary.csv"
  ),
  "\n"
)

cat(
  file.path(
    audit_dir,
    "03F_reviewer_country_summary.csv"
  ),
  "\n"
)

cat(
  file.path(
    audit_dir,
    "03H_candidate_nonGBIF_sara_points.csv"
  ),
  "\n"
)

cat(
  file.path(
    audit_dir,
    "03I_sara_only_cells.csv"
  ),
  "\n"
)

cat(
  "\nNo occurrence has been added to the model dataset.\n"
)

cat(
  "No model has been refitted.\n"
)

cat(
  "No REV1 file has been modified.\n"
)

cat(
  "============================================================\n"
)
