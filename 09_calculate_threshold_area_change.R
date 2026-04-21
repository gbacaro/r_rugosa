# 09_calculate_threshold_area_change.R
# Quantify current and future coastal area above the mean cross-validated suitability threshold.

suppressPackageStartupMessages({
  library(terra)
  library(readr)
  library(dplyr)
  library(tibble)
})

# ------------------------------------------------------------
# 1. Local project settings
# ------------------------------------------------------------
project_dir <- getwd()

dir_models <- file.path(project_dir, "models")
dir_tab    <- file.path(project_dir, "outputs", "tables")

current_sum_file  <- file.path(dir_tab, "table_current_prediction_summary.csv")
future_sum_file   <- file.path(dir_tab, "table_future_projection_summary.csv")
current_mean_file <- file.path(dir_models, "current_ensemble_mean.tif")
future_proj_dir   <- file.path(dir_models, "future_projections")

# ------------------------------------------------------------
# 2. Helper functions
# ------------------------------------------------------------
message_section <- function(txt) {
  cat("\n", paste(rep("=", nchar(txt) + 4), collapse = ""), "\n", sep = "")
  cat("  ", txt, "\n", sep = "")
  cat(paste(rep("=", nchar(txt) + 4), collapse = ""), "\n\n", sep = "")
}

fmt_num <- function(x, digits = 3) {
  if (length(x) == 0 || all(is.na(x))) return("NA")
  format(round(as.numeric(x), digits), nsmall = digits, trim = TRUE, scientific = FALSE)
}

fmt_num4 <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return("NA")
  format(round(as.numeric(x), 4), nsmall = 4, trim = TRUE, scientific = FALSE)
}

fmt_int <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return("NA")
  format(round(as.numeric(x), 0), big.mark = ",", scientific = FALSE, trim = TRUE)
}

print_kv <- function(name, value) {
  cat(sprintf("%-60s %s\n", paste0(name, ":"), value))
}

safe_global_sum <- function(r) {
  out <- tryCatch(terra::global(r, "sum", na.rm = TRUE), error = function(e) NULL)
  if (is.null(out)) return(NA_real_)
  as.numeric(out[1, 1])
}

# ------------------------------------------------------------
# 3. Read threshold and rasters
# ------------------------------------------------------------
if (!file.exists(current_sum_file))  stop("Missing file: ", current_sum_file)
if (!file.exists(future_sum_file))   stop("Missing file: ", future_sum_file)
if (!file.exists(current_mean_file)) stop("Missing file: ", current_mean_file)
if (!dir.exists(future_proj_dir))    stop("Missing directory: ", future_proj_dir)

current_sum <- readr::read_csv(current_sum_file, show_col_types = FALSE)
future_sum  <- readr::read_csv(future_sum_file, show_col_types = FALSE)

thr <- current_sum$threshold_mean[1]
if (!is.finite(thr)) stop("Threshold is missing or not finite in table_current_prediction_summary.csv")

current_r <- terra::rast(current_mean_file)

# ------------------------------------------------------------
# 4. Build valid coastal domain mask
# ------------------------------------------------------------
# Valid domain = cells where current ensemble has non-NA values
valid_domain_mask <- !is.na(current_r)

# Cell area in km2
area_r <- terra::cellSize(current_r, unit = "km")

# Restrict area raster to valid coastal cells only
area_valid <- terra::mask(area_r, current_r)

valid_domain_area_km2 <- safe_global_sum(area_valid)

# Sanity check from previous summary
expected_cell_prop <- current_sum$proportion_cells_above_threshold[1] * 100

# ------------------------------------------------------------
# 5. Current suitable area
# ------------------------------------------------------------
message_section("Current suitable coastal area above threshold")

current_bin <- current_r >= thr
current_bin <- terra::mask(current_bin, current_r)

current_suitable_area_km2 <- safe_global_sum(current_bin * area_valid)
current_suitable_area_pct <- 100 * current_suitable_area_km2 / valid_domain_area_km2

print_kv("Fixed threshold used", fmt_num4(thr))
print_kv("Valid coastal domain area (km2)", fmt_int(valid_domain_area_km2))
print_kv("Current suitable area (km2)", fmt_int(current_suitable_area_km2))
print_kv("Current suitable area (% of valid coastal domain)", fmt_num(current_suitable_area_pct, 2))
print_kv("Expected proportion above threshold from cell counts (%)", fmt_num(expected_cell_prop, 2))

# ------------------------------------------------------------
# 6. Loop over future scenarios
# ------------------------------------------------------------
message_section("Future suitable coastal area by scenario")

res_list <- list()

for (i in seq_len(nrow(future_sum))) {
  gcm <- future_sum$gcm[i]
  ssp <- future_sum$ssp[i]
  tt  <- future_sum$time_period[i]
  
  scen_tag <- paste(gcm, ssp, tt, sep = "_")
  fut_file <- file.path(
    future_proj_dir,
    scen_tag,
    paste0("future_ensemble_mean_", scen_tag, ".tif")
  )
  
  if (!file.exists(fut_file)) {
    warning("Missing file: ", fut_file)
    next
  }
  
  fut_r <- terra::rast(fut_file)
  
  if (!terra::compareGeom(current_r, fut_r, stopOnError = FALSE)) {
    stop("Geometry mismatch between current raster and future raster: ", fut_file)
  }
  
  # Force all future calculations onto the same valid coastal domain
  fut_r_masked <- terra::mask(fut_r, current_r)
  fut_bin <- fut_r_masked >= thr
  fut_bin <- terra::mask(fut_bin, current_r)
  
  fut_area_km2 <- safe_global_sum(fut_bin * area_valid)
  fut_area_pct <- 100 * fut_area_km2 / valid_domain_area_km2
  
  delta_area_km2 <- fut_area_km2 - current_suitable_area_km2
  delta_area_pct <- 100 * delta_area_km2 / current_suitable_area_km2
  
  res_list[[scen_tag]] <- tibble::tibble(
    gcm = gcm,
    ssp = ssp,
    time_period = tt,
    threshold = thr,
    valid_domain_area_km2 = valid_domain_area_km2,
    current_suitable_area_km2 = current_suitable_area_km2,
    current_suitable_area_pct = current_suitable_area_pct,
    future_suitable_area_km2  = fut_area_km2,
    future_suitable_area_pct  = fut_area_pct,
    delta_area_km2            = delta_area_km2,
    delta_area_pct            = delta_area_pct
  )
}

area_tbl <- dplyr::bind_rows(res_list)

# ------------------------------------------------------------
# 7. Save table
# ------------------------------------------------------------
out_file <- file.path(dir_tab, "table_suitable_area_above_threshold.csv")
readr::write_csv(area_tbl, out_file)

# ------------------------------------------------------------
# 8. Print scenario table
# ------------------------------------------------------------
print(area_tbl, n = nrow(area_tbl))

# ------------------------------------------------------------
# 9. Aggregated by SSP and time period
# ------------------------------------------------------------
message_section("Suitable area aggregated by SSP and time period")

agg_tbl <- area_tbl %>%
  group_by(ssp, time_period) %>%
  summarise(
    mean_future_area_km2 = mean(future_suitable_area_km2, na.rm = TRUE),
    mean_future_area_pct = mean(future_suitable_area_pct, na.rm = TRUE),
    mean_delta_area_km2  = mean(delta_area_km2, na.rm = TRUE),
    mean_delta_area_pct  = mean(delta_area_pct, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(ssp, time_period)

print(agg_tbl, n = nrow(agg_tbl))

# ------------------------------------------------------------
# 10. Best and worst scenarios
# ------------------------------------------------------------
message_section("Most positive and most negative area changes")

best_case <- area_tbl %>%
  arrange(desc(delta_area_pct)) %>%
  slice(1)

worst_case <- area_tbl %>%
  arrange(delta_area_pct) %>%
  slice(1)

cat("Most positive area change:\n")
print(best_case, n = 1)

cat("\nMost negative area change:\n")
print(worst_case, n = 1)

# ------------------------------------------------------------
# 11. Compact statements
# ------------------------------------------------------------
message_section("Compact paper-ready area statements")

print_kv("Average future area change across all scenarios (km2)", fmt_int(mean(area_tbl$delta_area_km2, na.rm = TRUE)))
print_kv("Average future area change across all scenarios (%)", fmt_num(mean(area_tbl$delta_area_pct, na.rm = TRUE), 2))
print_kv("Minimum area change across scenarios (%)", fmt_num(min(area_tbl$delta_area_pct, na.rm = TRUE), 2))
print_kv("Maximum area change across scenarios (%)", fmt_num(max(area_tbl$delta_area_pct, na.rm = TRUE), 2))

cat(
  "\nCurrent-area sentence:\n",
  "Using the mean cross-validated threshold (", fmt_num4(thr),
  "), the current ensemble classified ", fmt_int(current_suitable_area_km2),
  " km2 of the valid European coastal domain as climatically suitable, corresponding to ",
  fmt_num(current_suitable_area_pct, 2), "% of the analysed coastal area.\n",
  sep = ""
)

cat("\nSSP x period summary table:\n")
print(agg_tbl, n = nrow(agg_tbl))

cat("\nSaved table at:\n", out_file, "\n", sep = "")