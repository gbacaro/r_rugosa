# 04_download_chelsa_climate_coastal.R
# Download present and future CHELSA layers, crop them remotely, and mask them to the geographic coastal domains.

library(terra)
library(sf)
library(ggplot2)
library(curl)

project_dir <- getwd()

dir_clean     <- file.path(project_dir, "data_clean")
dir_raw_clim  <- file.path(project_dir, "data_raw", "climate_chelsa")
dir_fig       <- file.path(project_dir, "outputs", "figures")

dir.create(dir_raw_clim, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_fig, recursive = TRUE, showWarnings = FALSE)

eu_ext_vec   <- c(-15, 45, 30, 72)
asia_ext_vec <- c(118, 156, 28, 58)

# Local variable names used downstream in scripts 05-07
vars_use <- c("bio3", "bio8", "bio9", "bio10", "bio15", "bio18", "bio19", "fcf", "gsl", "gst")

gcm_list  <- c("GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL")
ssp_list  <- c("ssp126", "ssp585")
time_list <- c("2041-2070", "2071-2100")

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

run_shell_checked <- function(cmd_string) {
  cat("\nRunning command:\n", cmd_string, "\n\n")
  status <- system(cmd_string, intern = FALSE, ignore.stdout = FALSE, ignore.stderr = FALSE)
  if (!identical(status, 0L) && !identical(status, 0)) {
    stop(paste("Command failed:\n", cmd_string))
  }
  invisible(status)
}

# CHELSA present bioclim variables use names like bio1, bio8, bio10
# Additional BIOCLIM+ variables such as fcf, gsl, gst are already plain names
chelsa_present_var <- function(v) {
  v
}

# For future paths, keep the same local names unless CHELSA uses a different template.
chelsa_future_var <- function(v) {
  v
}

url_exists <- function(url) {
  ok <- FALSE
  try({
    h <- curl::new_handle(nobody = TRUE, followlocation = TRUE, connecttimeout = 20, timeout = 60)
    res <- curl::curl_fetch_memory(url, handle = h)
    if (!is.null(res$status_code) && res$status_code >= 200 && res$status_code < 400) {
      ok <- TRUE
    }
  }, silent = TRUE)
  ok
}

find_first_existing_url <- function(url_candidates, label = "resource") {
  for (u in url_candidates) {
    cat("Testing URL for", label, ":\n", u, "\n")
    if (url_exists(u)) {
      cat("Resolved URL for", label, ":\n", u, "\n")
      return(u)
    }
  }
  stop("No valid remote URL found for: ", label)
}

clip_remote_raster <- function(url, ext_vec, out_file, gdal_translate_bin) {
  if (file.exists(out_file)) return(out_file)
  
  xmin <- as.character(ext_vec[1])
  xmax <- as.character(ext_vec[2])
  ymin <- as.character(ext_vec[3])
  ymax <- as.character(ext_vec[4])
  
  remote_src <- paste0("/vsicurl/", url)
  
  cmd <- paste(
    shQuote(gdal_translate_bin),
    "-projwin", xmin, ymax, xmax, ymin,
    "-co", "COMPRESS=LZW",
    "-co", "TILED=YES",
    "-of", "GTiff",
    shQuote(remote_src),
    shQuote(normalizePath(out_file, winslash = "/", mustWork = FALSE))
  )
  
  run_shell_checked(cmd)
  
  if (!file.exists(out_file)) {
    stop("Output not created: ", out_file)
  }
  
  out_file
}

stack_and_write <- function(files, var_names, out_file) {
  r <- terra::rast(files)
  names(r) <- var_names
  terra::writeRaster(r, out_file, overwrite = TRUE)
  out_file
}

mask_stack_to_belt <- function(stack_file, belt_file, out_file) {
  r <- terra::rast(stack_file)
  v <- terra::vect(readRDS(belt_file))
  r2 <- terra::mask(terra::crop(r, v), v)
  terra::writeRaster(r2, out_file, overwrite = TRUE)
  out_file
}

message_section("Check GDAL")
gdal_translate_bin <- Sys.which("gdal_translate")
if (gdal_translate_bin == "") {
  stop("gdal_translate not found in PATH.")
}

message_section("Read coastal belts")
eu_belt_file   <- file.path(dir_clean, "europe_coastal_belt_1km_geographic.rds")
asia_belt_file <- file.path(dir_clean, "asia_coastal_belt_1km_geographic.rds")

if (!file.exists(eu_belt_file) || !file.exists(asia_belt_file)) {
  stop("Run 03_build_geographic_coastal_domain.R first.")
}

# Verified present base path
base_current <- "https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010"

# Candidate future base roots / patterns.
# We resolve the real one dynamically by testing URLs.
future_base_candidates <- c(
  "https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies",
  "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies"
)

# ------------------------------------------------------------
# CURRENT
# ------------------------------------------------------------
message_section("Download current CHELSA layers")

dir_cur <- file.path(dir_raw_clim, "current")
dir.create(dir_cur, recursive = TRUE, showWarnings = FALSE)

dir_cur_eu <- file.path(dir_cur, "europe_vars")
dir_cur_as <- file.path(dir_cur, "asia_vars")
dir.create(dir_cur_eu, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_cur_as, recursive = TRUE, showWarnings = FALSE)

eu_files <- character(length(vars_use))
as_files <- character(length(vars_use))

for (i in seq_along(vars_use)) {
  v_local  <- vars_use[i]
  v_remote <- chelsa_present_var(v_local)
  
  # For CHELSA present, use the verified root and plain variable token
  url <- paste0(
    base_current,
    "/bio/",
    "CHELSA_", v_remote, "_1981-2010_V.2.1.tif"
  )
  
  if (!url_exists(url)) {
    stop("Current CHELSA file not found for variable: ", v_local, "\nURL tested:\n", url)
  }
  
  eu_files[i] <- clip_remote_raster(
    url,
    eu_ext_vec,
    file.path(dir_cur_eu, paste0("CHELSA_", v_local, "_europe_1981_2010.tif")),
    gdal_translate_bin
  )
  
  as_files[i] <- clip_remote_raster(
    url,
    asia_ext_vec,
    file.path(dir_cur_as, paste0("CHELSA_", v_local, "_asia_1981_2010.tif")),
    gdal_translate_bin
  )
}

cur_eu_bbox <- file.path(dir_cur, "current_europe_bbox_selectedvars.tif")
cur_as_bbox <- file.path(dir_cur, "current_asia_bbox_selectedvars.tif")

stack_and_write(eu_files, vars_use, cur_eu_bbox)
stack_and_write(as_files, vars_use, cur_as_bbox)

cur_eu_coast <- file.path(dir_cur, "current_europe_coast1km_selectedvars.tif")
cur_as_coast <- file.path(dir_cur, "current_asia_coast1km_selectedvars.tif")

mask_stack_to_belt(cur_eu_bbox, eu_belt_file, cur_eu_coast)
mask_stack_to_belt(cur_as_bbox, asia_belt_file, cur_as_coast)

# ------------------------------------------------------------
# FUTURE
# ------------------------------------------------------------
message_section("Download future CHELSA layers")

dir_fut <- file.path(dir_raw_clim, "future")
dir.create(dir_fut, recursive = TRUE, showWarnings = FALSE)

model_slug <- function(x) tolower(x)

build_future_candidates <- function(v_remote, gcm, ssp, tt) {
  gcm_low <- model_slug(gcm)
  
  c(
    # candidate 1
    paste0(
      future_base_candidates[1], "/", tt, "/", gcm, "/", ssp, "/bio/",
      "CHELSA_", v_remote, "_", tt, "_", gcm_low, "_", ssp, "_V.2.1.tif"
    ),
    # candidate 2
    paste0(
      future_base_candidates[1], "/", tt, "/", gcm_low, "/", ssp, "/bio/",
      "CHELSA_", v_remote, "_", tt, "_", gcm_low, "_", ssp, "_V.2.1.tif"
    ),
    # candidate 3
    paste0(
      future_base_candidates[2], "/", tt, "/", gcm, "/", ssp, "/bio/",
      "CHELSA_", v_remote, "_", tt, "_", gcm_low, "_", ssp, "_V.2.1.tif"
    ),
    # candidate 4
    paste0(
      future_base_candidates[2], "/", tt, "/", gcm_low, "/", ssp, "/bio/",
      "CHELSA_", v_remote, "_", tt, "_", gcm_low, "_", ssp, "_V.2.1.tif"
    )
  )
}

for (gcm in gcm_list) {
  for (ssp in ssp_list) {
    for (tt in time_list) {
      
      message_section(paste("Future stack", gcm, ssp, tt))
      
      dir_combo <- file.path(dir_fut, paste(gcm, ssp, tt, sep = "_"))
      dir.create(dir_combo, recursive = TRUE, showWarnings = FALSE)
      
      fut_files <- character(length(vars_use))
      
      for (i in seq_along(vars_use)) {
        v_local  <- vars_use[i]
        v_remote <- chelsa_future_var(v_local)
        
        url <- find_first_existing_url(
          build_future_candidates(v_remote, gcm, ssp, tt),
          label = paste(v_local, gcm, ssp, tt)
        )
        
        fut_files[i] <- clip_remote_raster(
          url,
          eu_ext_vec,
          file.path(dir_combo, paste0("CHELSA_", v_local, "_", gcm, "_", ssp, "_", tt, "_europe.tif")),
          gdal_translate_bin
        )
      }
      
      fut_bbox <- file.path(dir_combo, paste0("future_europe_bbox_", gcm, "_", ssp, "_", tt, ".tif"))
      fut_coast <- file.path(dir_combo, paste0("future_europe_coast1km_", gcm, "_", ssp, "_", tt, ".tif"))
      
      stack_and_write(fut_files, vars_use, fut_bbox)
      mask_stack_to_belt(fut_bbox, eu_belt_file, fut_coast)
    }
  }
}

# ------------------------------------------------------------
# EDITORIAL DIAGNOSTICS
# ------------------------------------------------------------
message_section("Build editorial diagnostics")

r_eu <- terra::rast(cur_eu_coast)

plot_vars <- intersect(c("bio10", "bio18", "gsl"), names(r_eu))
df <- as.data.frame(r_eu[[plot_vars]], xy = TRUE, na.rm = TRUE)

make_raster_plot <- function(df, value_col, title_txt) {
  ggplot(df, aes(x = x, y = y, fill = .data[[value_col]])) +
    geom_raster() +
    coord_equal() +
    scale_fill_viridis_c(option = "magma") +
    labs(
      title = title_txt,
      x = "Longitude",
      y = "Latitude",
      fill = value_col
    ) +
    theme_editorial()
}

if ("bio10" %in% plot_vars) {
  save_editorial_plot(
    make_raster_plot(df, "bio10", "CHELSA current coastal climate: BIO10"),
    file.path(dir_fig, "fig11_current_climate_bio10_coast.png"),
    width = 7.5,
    height = 5
  )
}

if ("bio18" %in% plot_vars) {
  save_editorial_plot(
    make_raster_plot(df, "bio18", "CHELSA current coastal climate: BIO18"),
    file.path(dir_fig, "fig12_current_climate_bio18_coast.png"),
    width = 7.5,
    height = 5
  )
}

if ("gsl" %in% plot_vars) {
  save_editorial_plot(
    make_raster_plot(df, "gsl", "CHELSA current coastal climate: GSL"),
    file.path(dir_fig, "fig13_current_climate_gsl_coast.png"),
    width = 7.5,
    height = 5
  )
}

message_section("CHELSA coastal climate preparation completed")
cat("Current Europe coastal stack:\n", cur_eu_coast, "\n")
cat("Current Asia coastal stack:\n", cur_as_coast, "\n")
cat("Future Europe coastal stacks in:\n", dir_fut, "\n")