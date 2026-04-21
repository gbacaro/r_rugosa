# 03_build_geographic_coastal_domain.R
# Build the 1-km geographic coastal domains used to clip present and future CHELSA rasters.

suppressPackageStartupMessages({
  library(sf)
  library(readr)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(rnaturalearth)
  library(rnaturalearthdata)
})

project_dir <- getwd()

dir_clean <- file.path(project_dir, "data_clean")
dir_fig   <- file.path(project_dir, "outputs", "figures")
dir_tab   <- file.path(project_dir, "outputs", "tables")

dir.create(dir_clean, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_tab, recursive = TRUE, showWarnings = FALSE)

eu_ext_vec   <- c(-15, 45, 30, 72)
asia_ext_vec <- c(118, 156, 28, 58)
coastal_buffer_km <- 1

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

bbox_to_sf <- function(ext_vec) {
  sf::st_as_sfc(sf::st_bbox(c(
    xmin = ext_vec[1], xmax = ext_vec[2],
    ymin = ext_vec[3], ymax = ext_vec[4]
  ), crs = 4326))
}

build_coastal_belt <- function(ext_vec, buffer_km) {
  bb <- bbox_to_sf(ext_vec)

  coast <- rnaturalearth::ne_coastline(scale = "medium", returnclass = "sf") |>
    suppressWarnings(sf::st_intersection(bb))
  land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") |>
    suppressWarnings(sf::st_intersection(bb))

  coast_m <- sf::st_transform(coast, 3857)
  land_m  <- sf::st_transform(land, 3857)
  coast_buf <- sf::st_buffer(coast_m, dist = buffer_km * 1000)

  belt <- suppressWarnings(sf::st_intersection(sf::st_union(land_m), sf::st_union(coast_buf)))
  belt <- sf::st_make_valid(belt)
  belt <- suppressWarnings(sf::st_collection_extract(belt, "POLYGON"))
  sf::st_transform(belt, 4326)
}

belt_area_km2 <- function(x) {
  as.numeric(sum(sf::st_area(sf::st_transform(x, 6933)))) / 1e6
}

message_section("Build geographic coastal domains")

eu_belt_1km   <- build_coastal_belt(eu_ext_vec, coastal_buffer_km)
asia_belt_1km <- build_coastal_belt(asia_ext_vec, coastal_buffer_km)

saveRDS(eu_belt_1km,   file.path(dir_clean, "europe_coastal_belt_1km_geographic.rds"))
saveRDS(asia_belt_1km, file.path(dir_clean, "asia_coastal_belt_1km_geographic.rds"))

domain_tbl <- tibble::tibble(
  region = c("Europe", "East Asia"),
  coastal_buffer_km = coastal_buffer_km,
  xmin = c(eu_ext_vec[1], asia_ext_vec[1]),
  xmax = c(eu_ext_vec[2], asia_ext_vec[2]),
  ymin = c(eu_ext_vec[3], asia_ext_vec[3]),
  ymax = c(eu_ext_vec[4], asia_ext_vec[4]),
  coastal_belt_area_km2 = c(belt_area_km2(eu_belt_1km), belt_area_km2(asia_belt_1km))
)

readr::write_csv(
  domain_tbl,
  file.path(dir_tab, "table_geographic_coastal_domains.csv")
)

message_section("Build geographic-domain map")

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

p <- ggplot() +
  geom_sf(data = world, fill = "grey96", color = "grey75", linewidth = 0.2) +
  geom_sf(data = eu_belt_1km, fill = NA, color = "#2C7FB8", linewidth = 0.5) +
  geom_sf(data = asia_belt_1km, fill = NA, color = "#1A9641", linewidth = 0.5) +
  coord_sf(xlim = c(-20, 180), ylim = c(15, 75), expand = FALSE) +
  labs(
    title = "Geographic 1-km coastal domains",
    subtitle = "Domains used to clip CHELSA rasters before raster-valid filtering",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_editorial()

save_editorial_plot(
  p,
  file.path(dir_fig, "fig03_geographic_coastal_domains.png"),
  width = 10.5,
  height = 5.8
)

message_section("Geographic coastal domains saved")
print(domain_tbl)
