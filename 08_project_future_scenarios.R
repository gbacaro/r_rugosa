# 08_project_future_scenarios.R
# Project the fitted ensemble to CMIP6 CHELSA scenarios and export future suitability, uncertainty, consensus, and summary tables.

suppressPackageStartupMessages({
  library(terra)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(tibble)
  library(ranger)
})

# ------------------------------------------------------------
# 1. Local project settings
# ------------------------------------------------------------
project_dir <- getwd()

dir_models    <- file.path(project_dir, "models")
dir_fig       <- file.path(project_dir, "outputs", "figures")
dir_tab       <- file.path(project_dir, "outputs", "tables")
dir_raw_clim  <- file.path(project_dir, "data_raw", "climate_chelsa")
dir_temp      <- file.path(project_dir, "tmp_terra")

dir.create(dir_models, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_tab, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_temp, recursive = TRUE, showWarnings = FALSE)

terra::terraOptions(
  tempdir  = dir_temp,
  memfrac  = 0.60,
  progress = 1
)

# Future scenario settings: must match script 04
gcm_list  <- c("GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL")
ssp_list  <- c("ssp126", "ssp585")
time_list <- c("2041-2070", "2071-2100")

plot_max_cells     <- 250000
rf_predict_threads <- 1
cv_epsilon         <- 1e-6

# ------------------------------------------------------------
# 2. Helper functions
# ------------------------------------------------------------
message_section <- function(txt) {
  cat("\n", paste(rep("=", nchar(txt) + 4), collapse = ""), "\n", sep = "")
  cat("  ", txt, "\n", sep = "")
  cat(paste(rep("=", nchar(txt) + 4), collapse = ""), "\n\n", sep = "")
}

theme_editorial <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.25),
      axis.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 2),
      plot.subtitle = ggplot2::element_text(size = base_size),
      legend.position = "right"
    )
}

save_editorial_plot <- function(p, filename, width = 8, height = 5, dpi = 320) {
  ggplot2::ggsave(
    filename = filename,
    plot = p,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
}

extract_ranger_positive_prob <- function(pred_obj) {
  p <- pred_obj$predictions
  
  if (is.null(dim(p))) return(as.numeric(p))
  if (is.matrix(p) && ncol(p) == 1) return(as.numeric(p[, 1]))
  
  if (is.matrix(p) && ncol(p) >= 2) {
    cn <- colnames(p)
    
    if (!is.null(cn)) {
      if ("1" %in% cn)    return(as.numeric(p[, "1"]))
      if ("X1" %in% cn)   return(as.numeric(p[, "X1"]))
      if ("TRUE" %in% cn) return(as.numeric(p[, "TRUE"]))
      return(as.numeric(p[, ncol(p)]))
    }
    
    if (ncol(p) == 2) return(as.numeric(p[, 2]))
    return(as.numeric(p[, ncol(p)]))
  }
  
  stop("Unsupported ranger prediction structure.")
}

scale_block_df <- function(d, scaler_tbl, expected_names = NULL) {
  nd <- as.data.frame(d)
  
  if (!is.null(expected_names) && length(expected_names) == ncol(nd)) {
    names(nd) <- expected_names
  }
  
  for (i in seq_len(nrow(scaler_tbl))) {
    v <- scaler_tbl$variable[i]
    s <- scaler_tbl$sd[i]
    if (!is.finite(s) || s == 0) s <- 1
    nd[[v]] <- (as.numeric(nd[[v]]) - scaler_tbl$mean[i]) / s
  }
  
  nd
}

default_wopt <- function() {
  list(
    datatype = "FLT4S",
    gdal = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=IF_SAFER")
  )
}

aggregate_for_plot <- function(r, max_cells = 250000) {
  n <- terra::ncell(r)
  if (!is.finite(n) || n <= max_cells) return(r)
  fact <- ceiling(sqrt(n / max_cells))
  terra::aggregate(r, fact = fact, fun = mean, na.rm = TRUE)
}

remove_raster_artifacts <- function(filename) {
  targets <- unique(c(
    filename,
    paste0(filename, ".aux.xml"),
    paste0(filename, ".ovr"),
    paste0(filename, ".xml"),
    paste0(filename, ".vat.db")
  ))
  targets <- targets[file.exists(targets)]
  if (length(targets) > 0) {
    unlink(targets, force = TRUE)
  }
  invisible(TRUE)
}

make_unique_raster_filename <- function(dir_path, stem, ext = ".tif") {
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  file.path(
    dir_path,
    paste0(stem, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", sample.int(99999, 1), ext)
  )
}

safe_global_scalar <- function(r, fun) {
  out <- tryCatch(terra::global(r, fun = fun, na.rm = TRUE), error = function(e) NULL)
  if (is.null(out)) return(NA_real_)
  as.numeric(out[1, 1])
}

predict_algorithm_rast <- function(model, rast, scaler_tbl, algorithm, filename, rf_threads = 1) {
  vars_expected <- names(rast)
  
  fun_pred <- switch(
    algorithm,
    GLM = function(m, d) {
      nd <- scale_block_df(d, scaler_tbl, vars_expected)
      as.numeric(stats::predict(m, newdata = nd, type = "response"))
    },
    GAM = function(m, d) {
      nd <- scale_block_df(d, scaler_tbl, vars_expected)
      as.numeric(stats::predict(m, newdata = nd, type = "response"))
    },
    RF = function(m, d) {
      nd <- scale_block_df(d, scaler_tbl, vars_expected)
      pred_obj <- stats::predict(m, data = nd, num.threads = rf_threads)
      extract_ranger_positive_prob(pred_obj)
    },
    MAXNET = function(m, d) {
      nd <- scale_block_df(d, scaler_tbl, vars_expected)
      as.numeric(stats::predict(m, newdata = nd, type = "cloglog"))
    },
    stop("Unsupported algorithm: ", algorithm)
  )
  
  remove_raster_artifacts(filename)
  
  terra::predict(
    rast,
    model,
    fun = fun_pred,
    na.rm = TRUE,
    filename = filename,
    overwrite = TRUE,
    wopt = default_wopt()
  )
}

# ------------------------------------------------------------
# 3. Input files
# ------------------------------------------------------------
model_file        <- file.path(dir_models, "sdm_models_current.rds")
perf_file         <- file.path(dir_tab, "table_model_performance_cv.csv")
future_root       <- file.path(dir_raw_clim, "future")
current_mean_file <- file.path(dir_models, "current_ensemble_mean.tif")
current_sd_file   <- file.path(dir_models, "current_ensemble_sd.tif")
current_cv_file   <- file.path(dir_models, "current_ensemble_cv.tif")

if (!file.exists(model_file)) {
  stop("Missing file: ", model_file, "\nRun script 06 first.")
}
if (!file.exists(perf_file)) {
  stop("Missing file: ", perf_file, "\nRun script 06 first.")
}
if (!dir.exists(future_root)) {
  stop("Missing future climate directory: ", future_root, "\nRun script 04 first.")
}
if (!file.exists(current_mean_file) || !file.exists(current_sd_file) || !file.exists(current_cv_file)) {
  stop("Missing current ensemble rasters. Re-run script 06 first.")
}

# ------------------------------------------------------------
# 4. Read models and thresholds
# ------------------------------------------------------------
message_section("Read current SDM models and CV thresholds")

sdm_obj <- readRDS(model_file)
cv_tbl  <- readr::read_csv(perf_file, show_col_types = FALSE)

model_store <- sdm_obj$models
main_vars   <- sdm_obj$vars
scaler      <- sdm_obj$scaler

current_mean <- terra::rast(current_mean_file)
current_sd   <- terra::rast(current_sd_file)
current_cv   <- terra::rast(current_cv_file)

current_mean_global         <- safe_global_scalar(current_mean, mean)
current_suitability_sd      <- safe_global_scalar(current_mean, sd)
current_uncertainty_sd_mean <- safe_global_scalar(current_sd, mean)
current_uncertainty_cv_mean <- safe_global_scalar(current_cv, mean)

cat("Predictors used by models:\n")
print(main_vars)

thr_tbl <- cv_tbl |>
  dplyr::group_by(algorithm) |>
  dplyr::summarise(
    threshold_mean = mean(threshold, na.rm = TRUE),
    .groups = "drop"
  )

print(thr_tbl)

# ------------------------------------------------------------
# 5. Prepare output objects
# ------------------------------------------------------------
message_section("Prepare output directories")

dir_proj <- file.path(dir_models, "future_projections")
dir.create(dir_proj, recursive = TRUE, showWarnings = FALSE)

summary_rows <- list()

# ------------------------------------------------------------
# 6. Loop over future scenarios
# ------------------------------------------------------------
message_section("Project future coastal suitability")

for (gcm in gcm_list) {
  for (ssp in ssp_list) {
    for (tt in time_list) {
      
      scen_tag <- paste(gcm, ssp, tt, sep = "_")
      message_section(paste("Scenario", scen_tag))
      
      scen_dir <- file.path(future_root, scen_tag)
      fut_file <- file.path(scen_dir, paste0("future_europe_coast1km_", gcm, "_", ssp, "_", tt, ".tif"))
      
      if (!file.exists(fut_file)) {
        stop("Missing future raster: ", fut_file, "\nRun script 04 first.")
      }
      
      fut_stack <- terra::rast(fut_file)
      
      if (!all(main_vars %in% names(fut_stack))) {
        stop(
          "Future stack does not contain all model predictors.\nMissing: ",
          paste(setdiff(main_vars, names(fut_stack)), collapse = ", ")
        )
      }
      
      fut_stack <- fut_stack[[main_vars]]
      
      out_dir <- file.path(dir_proj, scen_tag)
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      
      dir_pred_tmp <- file.path(out_dir, "tmp_algorithm_predictions")
      dir.create(dir_pred_tmp, recursive = TRUE, showWarnings = FALSE)
      
      pred_files <- c(
        GLM    = make_unique_raster_filename(dir_pred_tmp, "future_glm"),
        GAM    = make_unique_raster_filename(dir_pred_tmp, "future_gam"),
        RF     = make_unique_raster_filename(dir_pred_tmp, "future_rf"),
        MAXNET = make_unique_raster_filename(dir_pred_tmp, "future_maxnet")
      )
      
      alg_file <- file.path(out_dir, paste0("future_predictions_by_algorithm_", scen_tag, ".tif"))
      ens_file <- file.path(out_dir, paste0("future_ensemble_mean_", scen_tag, ".tif"))
      sd_file  <- file.path(out_dir, paste0("future_ensemble_sd_", scen_tag, ".tif"))
      cv_file  <- file.path(out_dir, paste0("future_ensemble_cv_", scen_tag, ".tif"))
      agr_file <- file.path(out_dir, paste0("future_binary_agreement_", scen_tag, ".tif"))
      
      remove_raster_artifacts(alg_file)
      remove_raster_artifacts(ens_file)
      remove_raster_artifacts(sd_file)
      remove_raster_artifacts(cv_file)
      remove_raster_artifacts(agr_file)
      
      for (algo in names(pred_files)) {
        message_section(paste("Predict future suitability -", scen_tag, "-", algo))
        
        invisible(
          predict_algorithm_rast(
            model = model_store[[algo]],
            rast = fut_stack,
            scaler_tbl = scaler,
            algorithm = algo,
            filename = pred_files[[algo]],
            rf_threads = rf_predict_threads
          )
        )
        
        gc()
      }
      
      pred_stack <- terra::rast(unname(pred_files))
      names(pred_stack) <- names(pred_files)
      
      terra::writeRaster(
        pred_stack,
        alg_file,
        overwrite = TRUE,
        wopt = default_wopt()
      )
      
      ensemble_mean <- terra::app(
        pred_stack,
        mean,
        na.rm = TRUE,
        filename = ens_file,
        overwrite = TRUE,
        wopt = default_wopt()
      )
      
      ensemble_sd <- terra::app(
        pred_stack,
        stats::sd,
        na.rm = TRUE,
        filename = sd_file,
        overwrite = TRUE,
        wopt = default_wopt()
      )
      
      ensemble_cv <- terra::ifel(ensemble_mean > cv_epsilon, ensemble_sd / ensemble_mean, NA)
      terra::writeRaster(
        ensemble_cv,
        cv_file,
        overwrite = TRUE,
        wopt = default_wopt()
      )
      
      # Binary predictions using mean CV thresholds
      bin_layers <- list()
      for (algo in names(pred_files)) {
        thr <- thr_tbl |>
          dplyr::filter(algorithm == algo) |>
          dplyr::pull(threshold_mean)
        
        if (length(thr) == 0 || !is.finite(thr)) {
          thr <- 0.5
        }
        
        bin_layers[[algo]] <- terra::rast(pred_files[[algo]]) >= thr
      }
      
      bin_stack <- terra::rast(bin_layers)
      names(bin_stack) <- paste0(names(bin_layers), "_bin")
      
      agreement <- terra::app(bin_stack, sum, na.rm = TRUE)
      terra::writeRaster(
        agreement,
        agr_file,
        overwrite = TRUE,
        wopt = default_wopt()
      )
      
      # Scenario summary metrics
      scen_mean_suitability <- safe_global_scalar(ensemble_mean, mean)
      scen_sd_suitability   <- safe_global_scalar(ensemble_mean, sd)
      scen_unc_sd_mean      <- safe_global_scalar(ensemble_sd, mean)
      scen_unc_sd_max       <- safe_global_scalar(ensemble_sd, max)
      scen_unc_cv_mean      <- safe_global_scalar(ensemble_cv, mean)
      scen_unc_cv_max       <- safe_global_scalar(ensemble_cv, max)
      scen_agr_mean         <- safe_global_scalar(agreement, mean)
      scen_agr_max          <- safe_global_scalar(agreement, max)
      
      summary_rows[[scen_tag]] <- tibble::tibble(
        gcm                      = gcm,
        ssp                      = ssp,
        time_period              = tt,
        mean_suitability_mean    = scen_mean_suitability,
        mean_suitability_sd      = scen_sd_suitability,
        uncertainty_sd_mean      = scen_unc_sd_mean,
        uncertainty_sd_max       = scen_unc_sd_max,
        uncertainty_cv_mean      = scen_unc_cv_mean,
        uncertainty_cv_max       = scen_unc_cv_max,
        agreement_mean           = scen_agr_mean,
        agreement_max            = scen_agr_max,
        current_mean_suitability = current_mean_global,
        current_suitability_sd   = current_suitability_sd,
        current_uncertainty_sd_mean = current_uncertainty_sd_mean,
        current_uncertainty_cv_mean = current_uncertainty_cv_mean,
        delta_mean_suitability   = scen_mean_suitability - current_mean_global
      )
      
      # Editorial figures
      df_mean <- as.data.frame(aggregate_for_plot(ensemble_mean, max_cells = plot_max_cells), xy = TRUE, na.rm = TRUE)
      names(df_mean)[3] <- "suitability"
      
      p1 <- ggplot2::ggplot(df_mean, ggplot2::aes(x = x, y = y, fill = suitability)) +
        ggplot2::geom_tile() +
        ggplot2::coord_equal() +
        ggplot2::scale_fill_viridis_c(option = "magma") +
        ggplot2::labs(
          title = paste("Projected coastal suitability:", gcm, ssp, tt),
          subtitle = "Future ensemble mean",
          x = "Longitude",
          y = "Latitude",
          fill = "Suitability"
        ) +
        theme_editorial()
      
      save_editorial_plot(
        p1,
        file.path(out_dir, paste0("fig_future_ensemble_mean_", scen_tag, ".png")),
        width = 8.5,
        height = 6
      )
      
      df_sd <- as.data.frame(aggregate_for_plot(ensemble_sd, max_cells = plot_max_cells), xy = TRUE, na.rm = TRUE)
      names(df_sd)[3] <- "uncertainty"
      
      p2 <- ggplot2::ggplot(df_sd, ggplot2::aes(x = x, y = y, fill = uncertainty)) +
        ggplot2::geom_tile() +
        ggplot2::coord_equal() +
        ggplot2::scale_fill_viridis_c(option = "cividis") +
        ggplot2::labs(
          title = paste("Projection uncertainty:", gcm, ssp, tt),
          subtitle = "Future ensemble SD",
          x = "Longitude",
          y = "Latitude",
          fill = "SD"
        ) +
        theme_editorial()
      
      save_editorial_plot(
        p2,
        file.path(out_dir, paste0("fig_future_uncertainty_sd_", scen_tag, ".png")),
        width = 8.5,
        height = 6
      )
      
      df_agr <- as.data.frame(aggregate_for_plot(agreement, max_cells = plot_max_cells), xy = TRUE, na.rm = TRUE)
      names(df_agr)[3] <- "agreement"
      
      p3 <- ggplot2::ggplot(df_agr, ggplot2::aes(x = x, y = y, fill = agreement)) +
        ggplot2::geom_tile() +
        ggplot2::coord_equal() +
        ggplot2::scale_fill_viridis_c(option = "plasma", direction = -1) +
        ggplot2::labs(
          title = paste("Binary model agreement:", gcm, ssp, tt),
          subtitle = "Number of algorithms predicting suitability",
          x = "Longitude",
          y = "Latitude",
          fill = "Agreement"
        ) +
        theme_editorial()
      
      save_editorial_plot(
        p3,
        file.path(out_dir, paste0("fig_future_agreement_", scen_tag, ".png")),
        width = 8.5,
        height = 6
      )
      
      rm(pred_stack, ensemble_mean, ensemble_sd, ensemble_cv, agreement, fut_stack, bin_stack, bin_layers)
      gc()
    }
  }
}

# ------------------------------------------------------------
# 7. Save future summary table
# ------------------------------------------------------------
message_section("Save future summary table")

future_summary <- dplyr::bind_rows(summary_rows)

readr::write_csv(
  future_summary,
  file.path(dir_tab, "table_future_projection_summary.csv")
)

message_section("Future coastal projections completed")

cat("Saved future projection summary at:\n", file.path(dir_tab, "table_future_projection_summary.csv"), "\n")
cat("Saved future projections in:\n", dir_proj, "\n")