# Rosa rugosa coastal-climate workflow
# Run this script from the repository root.

species_name <- "Rosa rugosa"
project_dir  <- getwd()

dir_raw_gbif <- file.path(project_dir, "data_raw", "gbif")
dir_raw_clim <- file.path(project_dir, "data_raw", "climate_chelsa")
dir_clean    <- file.path(project_dir, "data_clean")
dir_models   <- file.path(project_dir, "models")
dir_fig      <- file.path(project_dir, "outputs", "figures")
dir_tab      <- file.path(project_dir, "outputs", "tables")
dir_paper    <- file.path(project_dir, "outputs", "paper_revision1")
dir_logs     <- file.path(project_dir, "logs")
dir_temp     <- file.path(project_dir, "tmp_terra")

all_dirs <- c(
  dir_raw_gbif, dir_raw_clim, dir_clean, dir_models,
  dir_fig, dir_tab, dir_paper, dir_logs, dir_temp
)
invisible(lapply(all_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

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
      legend.position = "right",
      strip.text = ggplot2::element_text(face = "bold")
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

project_settings <- list(
  species_name = species_name,
  project_dir = project_dir,
  dir_raw_gbif = dir_raw_gbif,
  dir_raw_clim = dir_raw_clim,
  dir_clean = dir_clean,
  dir_models = dir_models,
  dir_fig = dir_fig,
  dir_tab = dir_tab,
  dir_paper = dir_paper,
  dir_logs = dir_logs,
  dir_temp = dir_temp
)

saveRDS(project_settings, file.path(project_dir, "project_settings.rds"))

message_section("Project directories initialized")
print(project_settings)
