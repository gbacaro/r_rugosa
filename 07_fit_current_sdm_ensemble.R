# 07_fit_current_sdm_ensemble.R
# Fit the present-day coastal ensemble SDM, evaluate it with spatial cross-validation, and save current ensemble maps and performance tables.

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(tibble)
  library(pROC)
  library(mgcv)
  library(ranger)
  library(maxnet)
  library(blockCV)
  library(ecospat)
})

# ------------------------------------------------------------
# 1. Local project settings
# ------------------------------------------------------------
project_dir <- getwd()

dir_clean     <- file.path(project_dir, "data_clean")
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

bio_eu_file <- file.path(dir_raw_clim, "current", "current_europe_coast1km_selectedvars.tif")
occ_file    <- file.path(dir_clean, "occurrences_native_invaded_coastal.csv")
bg_file     <- file.path(dir_clean, "background_points_europe_coastal.csv")

candidate_vars          <- c("bio3", "bio8", "bio9", "bio10", "bio15", "bio18", "bio19", "fcf", "gsl", "gst")
main_cor_threshold      <- 0.7
cv_block_size_m         <- 150000
cv_folds                <- 5
max_presences_for_model <- 10000

plot_max_cells      <- 250000
rf_predict_threads  <- 1
cv_epsilon          <- 1e-6

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

thin_by_cell_df <- function(df, raster_template, lon = "decimal_longitude", lat = "decimal_latitude") {
  xy <- as.matrix(df[, c(lon, lat)])
  cell_id <- terra::cellFromXY(raster_template, xy)
  df$cell_id_tmp <- cell_id
  df <- df[!is.na(df$cell_id_tmp), , drop = FALSE]
  df <- df |>
    dplyr::group_by(cell_id_tmp) |>
    dplyr::slice(1) |>
    dplyr::ungroup()
  df$cell_id_tmp <- NULL
  df
}

fit_scaler <- function(df, vars) {
  tibble::tibble(
    variable = vars,
    mean = sapply(df[vars], mean, na.rm = TRUE),
    sd   = sapply(df[vars], sd, na.rm = TRUE)
  )
}

apply_scaler_df <- function(df, scaler_tbl) {
  out <- df
  for (i in seq_len(nrow(scaler_tbl))) {
    v <- scaler_tbl$variable[i]
    s <- scaler_tbl$sd[i]
    if (!is.finite(s) || s == 0) s <- 1
    out[[v]] <- (out[[v]] - scaler_tbl$mean[i]) / s
  }
  out
}

drop_correlated_vars <- function(df, vars, threshold = 0.7) {
  cor_mat <- stats::cor(df[, vars], use = "pairwise.complete.obs")
  cor_abs <- abs(cor_mat)
  diag(cor_abs) <- 0
  
  keep <- vars
  repeat {
    max_cor <- max(cor_abs[keep, keep], na.rm = TRUE)
    if (!is.finite(max_cor) || max_cor < threshold) break
    mean_cor <- colMeans(cor_abs[keep, keep], na.rm = TRUE)
    remove_var <- names(which.max(mean_cor))
    keep <- setdiff(keep, remove_var)
  }
  keep
}

extract_test_ids <- function(x) {
  if (is.integer(x) || is.numeric(x)) return(as.integer(x))
  if (is.list(x) && "test" %in% names(x)) return(as.integer(x$test))
  if (is.list(x) && "test_ids" %in% names(x)) return(as.integer(x$test_ids))
  if (is.list(x) && length(x) >= 2) return(as.integer(x[[2]]))
  stop("Could not extract test indices from blockCV output.")
}

calc_metrics <- function(obs, pred) {
  pred <- as.numeric(pred)
  
  if (length(unique(obs)) < 2) {
    return(c(
      auc = NA_real_,
      tss = NA_real_,
      threshold = NA_real_,
      sensitivity = NA_real_,
      specificity = NA_real_
    ))
  }
  
  roc_obj <- pROC::roc(obs, pred, quiet = TRUE)
  auc_val <- as.numeric(pROC::auc(roc_obj))
  
  best <- tryCatch(
    pROC::coords(
      roc_obj,
      x = "best",
      best.method = "youden",
      ret = c("threshold", "sensitivity", "specificity"),
      transpose = FALSE
    ),
    error = function(e) NULL
  )
  
  if (is.null(best)) {
    return(c(
      auc = auc_val,
      tss = NA_real_,
      threshold = NA_real_,
      sensitivity = NA_real_,
      specificity = NA_real_
    ))
  }
  
  if (is.data.frame(best)) {
    thr  <- as.numeric(best$threshold[1])
    sens <- as.numeric(best$sensitivity[1])
    spec <- as.numeric(best$specificity[1])
  } else {
    thr  <- as.numeric(best["threshold"])
    sens <- as.numeric(best["sensitivity"])
    spec <- as.numeric(best["specificity"])
  }
  
  c(
    auc = auc_val,
    tss = sens + spec - 1,
    threshold = thr,
    sensitivity = sens,
    specificity = spec
  )
}

calc_boyce_vector <- function(fit_all, fit_pres) {
  fit_all  <- as.numeric(fit_all)
  fit_pres <- as.numeric(fit_pres)
  
  fit_all  <- fit_all[is.finite(fit_all)]
  fit_pres <- fit_pres[is.finite(fit_pres)]
  
  if (length(fit_all) < 20 || length(fit_pres) < 10) return(NA_real_)
  if (length(unique(fit_all)) < 10 || length(unique(fit_pres)) < 5) return(NA_real_)
  
  out <- tryCatch(
    ecospat::ecospat.boyce(
      fit = fit_all,
      obs = fit_pres,
      nclass = 0,
      window.w = "default",
      res = 100,
      PEplot = FALSE,
      rm.duplicate = FALSE
    ),
    error = function(e) NULL
  )
  
  if (is.null(out)) return(NA_real_)
  
  if (!is.null(out$cor) && length(out$cor) > 0 && is.finite(out$cor[1])) {
    return(as.numeric(out$cor[1]))
  }
  
  if (!is.null(out$Spearman.cor) && length(out$Spearman.cor) > 0 && is.finite(out$Spearman.cor[1])) {
    return(as.numeric(out$Spearman.cor[1]))
  }
  
  NA_real_
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

safe_global_scalar <- function(r, fun) {
  out <- tryCatch(terra::global(r, fun = fun, na.rm = TRUE), error = function(e) NULL)
  if (is.null(out)) return(NA_real_)
  as.numeric(out[1, 1])
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
# 3. Check input files
# ------------------------------------------------------------
if (!file.exists(occ_file)) stop("Missing file: ", occ_file, "\nRun script 03 first.")
if (!file.exists(bg_file))  stop("Missing file: ", bg_file, "\nRun script 03 first.")
if (!file.exists(bio_eu_file)) stop("Missing file: ", bio_eu_file, "\nRun script 04 first.")

# ------------------------------------------------------------
# 4. Read data
# ------------------------------------------------------------
message_section("Read coastal occurrence, background and CHELSA climate data")

occ <- readr::read_csv(occ_file, show_col_types = FALSE)
bg  <- readr::read_csv(bg_file, show_col_types = FALSE)
bio_eu <- terra::rast(bio_eu_file)

cat("Europe coastal climate layers:\n")
print(names(bio_eu))

# ------------------------------------------------------------
# 5. Prepare coastal invaded presences
# ------------------------------------------------------------
message_section("Prepare invaded European coastal presences")

pres <- occ |>
  dplyr::filter(range_status == "invaded") |>
  dplyr::select(decimal_longitude, decimal_latitude) |>
  dplyr::mutate(pa = 1)

cat("Invaded European coastal presences before thinning:", nrow(pres), "\n")
pres <- thin_by_cell_df(pres, bio_eu[[1]])
cat("Invaded European coastal presences after thinning:", nrow(pres), "\n")

# ------------------------------------------------------------
# 6. Extract CHELSA values
# ------------------------------------------------------------
message_section("Extract CHELSA coastal climate values")

candidate_vars <- intersect(candidate_vars, names(bio_eu))
if (length(candidate_vars) < 2) {
  stop("Too few candidate variables found in the CHELSA coastal raster.")
}

pres_vals <- terra::extract(
  bio_eu[[candidate_vars]],
  pres[, c("decimal_longitude", "decimal_latitude")]
)[, -1, drop = FALSE]

bg_vals <- terra::extract(
  bio_eu[[candidate_vars]],
  bg[, c("decimal_longitude", "decimal_latitude")]
)[, -1, drop = FALSE]

pres_df <- cbind(pres, pres_vals) |> stats::na.omit()
bg_df   <- cbind(bg, bg_vals) |> stats::na.omit()

cat("Presence records after climate extraction:", nrow(pres_df), "\n")
cat("Background records after climate extraction:", nrow(bg_df), "\n")

if (nrow(pres_df) > max_presences_for_model) {
  set.seed(123)
  pres_df <- pres_df |>
    dplyr::slice_sample(n = max_presences_for_model)
}

cat("Presence records used for modelling:", nrow(pres_df), "\n")
cat("Background records used for modelling:", nrow(bg_df), "\n")

# ------------------------------------------------------------
# 7. Select final predictor set
# ------------------------------------------------------------
message_section("Select final climate predictor set")

pooled_env <- rbind(
  pres_df[, candidate_vars, drop = FALSE],
  bg_df[, candidate_vars, drop = FALSE]
)

main_vars <- drop_correlated_vars(
  pooled_env,
  candidate_vars,
  threshold = main_cor_threshold
)

cat("Final predictors retained:\n")
print(main_vars)

readr::write_csv(
  tibble::tibble(variable = main_vars),
  file.path(dir_tab, "table_main_predictors.csv")
)

# ------------------------------------------------------------
# 8. Build modelling table and standardize predictors
# ------------------------------------------------------------
message_section("Build modelling table")

model_df <- rbind(
  pres_df[, c("pa", "decimal_longitude", "decimal_latitude", main_vars)],
  bg_df[, c("pa", "decimal_longitude", "decimal_latitude", main_vars)]
)

scaler <- fit_scaler(model_df, main_vars)
model_df_scaled <- apply_scaler_df(model_df, scaler)

saveRDS(
  scaler,
  file.path(dir_models, "predictor_scaler_present.rds")
)

cat("Total modelling records:", nrow(model_df_scaled), "\n")

# ------------------------------------------------------------
# 9. Spatial cross-validation
# ------------------------------------------------------------
message_section("Create spatial cross-validation folds")

model_sf <- sf::st_as_sf(
  model_df_scaled,
  coords = c("decimal_longitude", "decimal_latitude"),
  crs = 4326
)

model_sf_3035 <- sf::st_transform(model_sf, 3035)

set.seed(100)
cv <- blockCV::cv_spatial(
  x = model_sf_3035,
  column = "pa",
  size = cv_block_size_m,
  k = cv_folds,
  progress = FALSE
)

fold_id <- rep(NA_integer_, nrow(model_df_scaled))
for (k in seq_along(cv$folds_list)) {
  test_ids <- extract_test_ids(cv$folds_list[[k]])
  fold_id[test_ids] <- k
}
model_df_scaled$fold <- fold_id

fold_sf <- sf::st_as_sf(
  model_df_scaled,
  coords = c("decimal_longitude", "decimal_latitude"),
  crs = 4326
)

p_fold <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = fold_sf, ggplot2::aes(color = factor(fold)), alpha = 0.6, size = 0.45) +
  ggplot2::labs(
    title = "Spatial cross-validation folds",
    x = "Longitude",
    y = "Latitude",
    color = "Fold"
  ) +
  theme_editorial()

save_editorial_plot(
  p_fold,
  file.path(dir_fig, "fig08_spatial_cv_folds.png"),
  width = 8.5,
  height = 6
)

# ------------------------------------------------------------
# 10. Model functions
# ------------------------------------------------------------
fit_glm <- function(train_df, vars) {
  f <- stats::as.formula(paste("pa ~", paste(vars, collapse = " + ")))
  stats::glm(f, data = train_df, family = stats::binomial())
}

pred_glm <- function(model, new_df) {
  stats::predict(model, newdata = new_df, type = "response")
}

fit_gam <- function(train_df, vars) {
  smooth_terms <- paste0("s(", vars, ", k = 4)", collapse = " + ")
  f <- stats::as.formula(paste("pa ~", smooth_terms))
  mgcv::gam(f, data = train_df, family = stats::binomial())
}

pred_gam <- function(model, new_df) {
  stats::predict(model, newdata = new_df, type = "response")
}

fit_rf <- function(train_df, vars) {
  train_rf <- train_df
  train_rf$pa <- factor(train_rf$pa, levels = c(0, 1))
  
  ranger::ranger(
    formula = stats::as.formula(paste("pa ~", paste(vars, collapse = " + "))),
    data = train_rf,
    probability = TRUE,
    num.trees = 500,
    importance = "impurity"
  )
}

pred_rf <- function(model, new_df) {
  pred_obj <- stats::predict(model, data = new_df, num.threads = rf_predict_threads)
  extract_ranger_positive_prob(pred_obj)
}

fit_maxnet <- function(train_df, vars) {
  maxnet::maxnet(
    p = train_df$pa,
    data = train_df[, vars, drop = FALSE],
    f = maxnet::maxnet.formula(train_df$pa, train_df[, vars, drop = FALSE], classes = "lq")
  )
}

pred_maxnet <- function(model, new_df, vars) {
  as.numeric(stats::predict(model, newdata = new_df[, vars, drop = FALSE], type = "cloglog"))
}

# ------------------------------------------------------------
# 11. Cross-validated model fitting
# ------------------------------------------------------------
message_section("Run cross-validated model fitting")

algo_names <- c("GLM", "GAM", "RF", "MAXNET")
cv_results <- list()
model_store <- list()

for (algo in algo_names) {
  cat("Running:", algo, "\n")
  fold_res <- list()
  
  for (k in sort(stats::na.omit(unique(model_df_scaled$fold)))) {
    train_df <- as.data.frame(model_df_scaled |> dplyr::filter(fold != k))
    test_df  <- as.data.frame(model_df_scaled |> dplyr::filter(fold == k))
    
    if (algo == "GLM") {
      mod <- fit_glm(train_df, main_vars)
      pr  <- pred_glm(mod, test_df)
    }
    
    if (algo == "GAM") {
      mod <- fit_gam(train_df, main_vars)
      pr  <- pred_gam(mod, test_df)
    }
    
    if (algo == "RF") {
      mod <- fit_rf(train_df, main_vars)
      pr  <- pred_rf(mod, test_df)
    }
    
    if (algo == "MAXNET") {
      mod <- fit_maxnet(train_df, main_vars)
      pr  <- pred_maxnet(mod, test_df, main_vars)
    }
    
    metrics <- calc_metrics(test_df$pa, pr)
    boyce_val <- calc_boyce_vector(
      fit_all  = pr,
      fit_pres = pr[test_df$pa == 1]
    )
    
    fold_res[[as.character(k)]] <- tibble::tibble(
      algorithm   = algo,
      fold        = k,
      auc         = as.numeric(metrics["auc"]),
      tss         = as.numeric(metrics["tss"]),
      threshold   = as.numeric(metrics["threshold"]),
      sensitivity = as.numeric(metrics["sensitivity"]),
      specificity = as.numeric(metrics["specificity"]),
      boyce       = as.numeric(boyce_val)
    )
  }
  
  cv_results[[algo]] <- dplyr::bind_rows(fold_res)
  
  full_df <- as.data.frame(model_df_scaled)
  
  if (algo == "GLM")    model_store[[algo]] <- fit_glm(full_df, main_vars)
  if (algo == "GAM")    model_store[[algo]] <- fit_gam(full_df, main_vars)
  if (algo == "RF")     model_store[[algo]] <- fit_rf(full_df, main_vars)
  if (algo == "MAXNET") model_store[[algo]] <- fit_maxnet(full_df, main_vars)
  
  rm(full_df)
  gc()
}

cv_tbl <- dplyr::bind_rows(cv_results)

readr::write_csv(
  cv_tbl,
  file.path(dir_tab, "table_model_performance_cv.csv")
)

perf_summary <- cv_tbl |>
  dplyr::group_by(algorithm) |>
  dplyr::summarise(
    auc_mean         = mean(auc, na.rm = TRUE),
    auc_sd           = sd(auc, na.rm = TRUE),
    tss_mean         = mean(tss, na.rm = TRUE),
    tss_sd           = sd(tss, na.rm = TRUE),
    threshold_mean   = mean(threshold, na.rm = TRUE),
    threshold_sd     = sd(threshold, na.rm = TRUE),
    sensitivity_mean = mean(sensitivity, na.rm = TRUE),
    sensitivity_sd   = sd(sensitivity, na.rm = TRUE),
    specificity_mean = mean(specificity, na.rm = TRUE),
    specificity_sd   = sd(specificity, na.rm = TRUE),
    boyce_mean       = mean(as.numeric(boyce), na.rm = TRUE),
    boyce_sd         = sd(as.numeric(boyce), na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(
  perf_summary,
  file.path(dir_tab, "table_model_performance_summary.csv")
)

saveRDS(
  list(
    models   = model_store,
    vars     = main_vars,
    scaler   = scaler,
    cv       = cv_tbl,
    summary  = perf_summary
  ),
  file.path(dir_models, "sdm_models_current.rds")
)

# Release objects no longer needed before raster prediction
rm(
  pres_vals, bg_vals, pooled_env,
  model_df, model_df_scaled,
  model_sf, model_sf_3035, fold_sf, cv
)
gc()

# ------------------------------------------------------------
# 12. Editorial performance figure
# ------------------------------------------------------------
message_section("Build model performance figure")

p_perf <- ggplot2::ggplot(cv_tbl, ggplot2::aes(x = algorithm, y = auc)) +
  ggplot2::geom_boxplot() +
  ggplot2::labs(
    title = "Cross-validated model performance",
    subtitle = "AUC across spatial folds",
    x = NULL,
    y = "AUC"
  ) +
  theme_editorial()

save_editorial_plot(
  p_perf,
  file.path(dir_fig, "fig09_model_performance_auc.png"),
  width = 7.5,
  height = 5
)

# ------------------------------------------------------------
# 13. Predict current coastal suitability
# ------------------------------------------------------------
message_section("Predict current coastal ensemble suitability")

current_stack <- bio_eu[[main_vars]]

dir_pred_tmp <- file.path(dir_models, "tmp_current_algorithm_predictions")
dir.create(dir_pred_tmp, recursive = TRUE, showWarnings = FALSE)

pred_files <- c(
  GLM    = make_unique_raster_filename(dir_pred_tmp, "current_glm"),
  GAM    = make_unique_raster_filename(dir_pred_tmp, "current_gam"),
  RF     = make_unique_raster_filename(dir_pred_tmp, "current_rf"),
  MAXNET = make_unique_raster_filename(dir_pred_tmp, "current_maxnet")
)

final_alg_file  <- file.path(dir_models, "current_predictions_by_algorithm.tif")
final_mean_file <- file.path(dir_models, "current_ensemble_mean.tif")
final_sd_file   <- file.path(dir_models, "current_ensemble_sd.tif")
final_cv_file   <- file.path(dir_models, "current_ensemble_cv.tif")

remove_raster_artifacts(final_alg_file)
remove_raster_artifacts(final_mean_file)
remove_raster_artifacts(final_sd_file)
remove_raster_artifacts(final_cv_file)

for (algo in names(pred_files)) {
  message_section(paste("Predict current suitability -", algo))
  
  invisible(
    predict_algorithm_rast(
      model = model_store[[algo]],
      rast = current_stack,
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
  final_alg_file,
  overwrite = TRUE,
  wopt = default_wopt()
)

ensemble_mean <- terra::app(
  pred_stack,
  mean,
  na.rm = TRUE,
  filename = final_mean_file,
  overwrite = TRUE,
  wopt = default_wopt()
)

ensemble_sd <- terra::app(
  pred_stack,
  stats::sd,
  na.rm = TRUE,
  filename = final_sd_file,
  overwrite = TRUE,
  wopt = default_wopt()
)

ensemble_cv <- terra::ifel(ensemble_mean > cv_epsilon, ensemble_sd / ensemble_mean, NA)
terra::writeRaster(
  ensemble_cv,
  final_cv_file,
  overwrite = TRUE,
  wopt = default_wopt()
)

thr_overall <- mean(cv_tbl$threshold, na.rm = TRUE)
prop_above <- safe_global_scalar(ensemble_mean >= thr_overall, mean)

current_summary <- tibble::tibble(
  mean_suitability                 = safe_global_scalar(ensemble_mean, mean),
  sd_suitability                   = safe_global_scalar(ensemble_mean, sd),
  min_suitability                  = safe_global_scalar(ensemble_mean, min),
  max_suitability                  = safe_global_scalar(ensemble_mean, max),
  uncertainty_sd_mean              = safe_global_scalar(ensemble_sd, mean),
  uncertainty_sd_max               = safe_global_scalar(ensemble_sd, max),
  uncertainty_cv_mean              = safe_global_scalar(ensemble_cv, mean),
  uncertainty_cv_max               = safe_global_scalar(ensemble_cv, max),
  threshold_mean                   = thr_overall,
  proportion_cells_above_threshold = prop_above
)

readr::write_csv(
  current_summary,
  file.path(dir_tab, "table_current_prediction_summary.csv")
)

rm(pred_stack, current_stack)
gc()

# ------------------------------------------------------------
# 14. Editorial coastal suitability figure
# ------------------------------------------------------------
message_section("Build current coastal ensemble figure")

ensemble_mean_plot <- aggregate_for_plot(ensemble_mean, max_cells = plot_max_cells)
cur_df <- as.data.frame(ensemble_mean_plot, xy = TRUE, na.rm = TRUE)
names(cur_df)[3] <- "suitability"

p_cur <- ggplot2::ggplot(cur_df, ggplot2::aes(x = x, y = y, fill = suitability)) +
  ggplot2::geom_tile() +
  ggplot2::coord_equal() +
  ggplot2::scale_fill_viridis_c(option = "magma") +
  ggplot2::labs(
    title = "Current coastal suitability for Rosa rugosa in Europe",
    x = "Longitude",
    y = "Latitude",
    fill = "Suitability"
  ) +
  theme_editorial()

save_editorial_plot(
  p_cur,
  file.path(dir_fig, "fig10_current_ensemble_suitability.png"),
  width = 8.5,
  height = 6
)

message_section("Current coastal SDM completed")

cat("Saved model performance fold-level table at:\n", file.path(dir_tab, "table_model_performance_cv.csv"), "\n")
cat("Saved model performance summary at:\n", file.path(dir_tab, "table_model_performance_summary.csv"), "\n")
cat("Saved current prediction summary at:\n", file.path(dir_tab, "table_current_prediction_summary.csv"), "\n")
cat("Saved model object at:\n", file.path(dir_models, "sdm_models_current.rds"), "\n")
cat("Saved current algorithm stack at:\n", final_alg_file, "\n")
cat("Saved current coastal ensemble raster at:\n", final_mean_file, "\n")
cat("Saved current ensemble SD raster at:\n", final_sd_file, "\n")
cat("Saved current ensemble CV raster at:\n", final_cv_file, "\n")
cat("Saved figures in:\n", dir_fig, "\n")