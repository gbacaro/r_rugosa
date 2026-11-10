# ============================================================
# ROSA RUGOSA - REVISION 2
# 00_setup_REV2.R
#
# Create a minimal, clean REV2 workspace.
# REV1 is treated as read-only.
# ============================================================

rm(list = ls())
gc()

# ------------------------------------------------------------
# 1. Locate REV2 project
# ------------------------------------------------------------

rev2_dir <- NA_character_

if (requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) {
  rev2_dir <- tryCatch(
    rstudioapi::getActiveProject(),
    error = function(e) NA_character_
  )
}

if (is.na(rev2_dir) || !nzchar(rev2_dir)) {
  rev2_dir <- "C:/Users/works/Desktop/ROSA RUGOSA REV 2"
}

rev2_dir <- normalizePath(rev2_dir, winslash = "/", mustWork = TRUE)

# ------------------------------------------------------------
# 2. REV1 archive paths - READ ONLY
# ------------------------------------------------------------

rev1_dir <- normalizePath(
  "C:/Users/works/Desktop/PAPER ROSA RUGOSA/Rosa rugosa bis/REVISIONE_1_ROSA",
  winslash = "/",
  mustWork = TRUE
)

rev1_rdata_general <- file.path(rev1_dir, "dati_per_revisione.RData")
rev1_rdata_7pred   <- file.path(rev1_dir, "dati_revisione_7_predittori.RData")

# ------------------------------------------------------------
# 3. Minimal REV2 folder structure
# ------------------------------------------------------------

folders <- c(
  "01_input",
  "02_scripts",
  "03_audit",
  "04_work",
  "05_output",
  "06_submission"
)

for (f in folders) {
  dir.create(file.path(rev2_dir, f), recursive = TRUE, showWarnings = FALSE)
}

input_dir      <- file.path(rev2_dir, "01_input")
scripts_dir    <- file.path(rev2_dir, "02_scripts")
audit_dir      <- file.path(rev2_dir, "03_audit")
work_dir       <- file.path(rev2_dir, "04_work")
output_dir     <- file.path(rev2_dir, "05_output")
submission_dir <- file.path(rev2_dir, "06_submission")

# ------------------------------------------------------------
# 4. Find and extract Sara/Ribotta archive
# ------------------------------------------------------------

zip_candidates <- list.files(
  rev2_dir,
  pattern = "^R__Major_Revisions_requested_BINV-D-26-00324R1.*\\.zip$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(zip_candidates) == 0) {
  stop(
    "Sara/Ribotta ZIP not found in REV2 root. Expected a file named like ",
    "R__Major_Revisions_requested_BINV-D-26-00324R1.zip"
  )
}

if (length(zip_candidates) > 1) {
  warning("More than one Sara/Ribotta ZIP found. Using: ", basename(zip_candidates[1]))
}

sara_zip <- normalizePath(zip_candidates[1], winslash = "/", mustWork = TRUE)
sara_dir <- file.path(input_dir, "sara_ribotta")

dir.create(sara_dir, recursive = TRUE, showWarnings = FALSE)

old_sara <- list.files(sara_dir, full.names = TRUE, all.files = TRUE, no.. = TRUE)
if (length(old_sara) > 0) {
  unlink(old_sara, recursive = TRUE, force = TRUE)
}

unzip(sara_zip, exdir = sara_dir, overwrite = TRUE)

# ------------------------------------------------------------
# 5. Write central path configuration
# ------------------------------------------------------------

paths_file <- file.path(scripts_dir, "REV2_paths.R")

cfg <- c(
  "# Generated automatically by 00_setup_REV2.R",
  sprintf('rev2_dir <- "%s"', rev2_dir),
  sprintf('rev1_dir <- "%s"', rev1_dir),
  sprintf('rev1_rdata_general <- "%s"', normalizePath(rev1_rdata_general, winslash = "/", mustWork = FALSE)),
  sprintf('rev1_rdata_7pred <- "%s"', normalizePath(rev1_rdata_7pred, winslash = "/", mustWork = FALSE)),
  sprintf('sara_zip <- "%s"', sara_zip),
  sprintf('sara_dir <- "%s"', normalizePath(sara_dir, winslash = "/", mustWork = FALSE)),
  'input_dir <- file.path(rev2_dir, "01_input")',
  'scripts_dir <- file.path(rev2_dir, "02_scripts")',
  'audit_dir <- file.path(rev2_dir, "03_audit")',
  'work_dir <- file.path(rev2_dir, "04_work")',
  'output_dir <- file.path(rev2_dir, "05_output")',
  'submission_dir <- file.path(rev2_dir, "06_submission")'
)

writeLines(cfg, paths_file)

# ------------------------------------------------------------
# 6. Input and structure checks
# ------------------------------------------------------------

shp_files <- list.files(
  sara_dir,
  pattern = "\\.shp$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

input_check <- data.frame(
  item = c(
    "REV2 project",
    "REV1 archive",
    "REV1 general RData",
    "REV1 7-predictor RData",
    "Sara/Ribotta ZIP",
    "Sara/Ribotta shapefile"
  ),
  path = c(
    rev2_dir,
    rev1_dir,
    rev1_rdata_general,
    rev1_rdata_7pred,
    sara_zip,
    if (length(shp_files) > 0) shp_files[1] else NA_character_
  ),
  exists = c(
    dir.exists(rev2_dir),
    dir.exists(rev1_dir),
    file.exists(rev1_rdata_general),
    file.exists(rev1_rdata_7pred),
    file.exists(sara_zip),
    length(shp_files) > 0
  ),
  stringsAsFactors = FALSE
)

write.csv(
  input_check,
  file.path(audit_dir, "00_input_check.csv"),
  row.names = FALSE
)

structure <- data.frame(
  folder = folders,
  path = file.path(rev2_dir, folders),
  exists = dir.exists(file.path(rev2_dir, folders)),
  stringsAsFactors = FALSE
)

write.csv(
  structure,
  file.path(audit_dir, "00_folder_structure.csv"),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(audit_dir, "00_sessionInfo.txt")
)

# ------------------------------------------------------------
# 7. Report
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("ROSA RUGOSA REV2 - CLEAN WORKSPACE INITIALIZED\n")
cat("============================================================\n\n")
cat("REV2 root:\n", rev2_dir, "\n\n")
cat("REV1 archive (READ ONLY):\n", rev1_dir, "\n\n")
cat("Input check:\n")
print(input_check, row.names = FALSE)
cat("\nFolders created:\n")
print(structure[, c("folder", "exists")], row.names = FALSE)
cat("\nCentral path file:\n", paths_file, "\n")
cat("\nNo REV1 file has been modified. No model has been fitted.\n")
cat("============================================================\n")
