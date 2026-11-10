# ============================================================
# ROSA RUGOSA - REVISION 2
# 02_model_input_audit_REV2.R
#
# Audit how presences/background were encoded and used by the
# four REV1 modelling algorithms.
#
# No model is fitted. REV1 remains read only.
# ============================================================

rm(list = ls())
gc()

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

paths_file <- file.path(getwd(), "02_scripts", "REV2_paths.R")
if (!file.exists(paths_file)) {
  stop("Missing 02_scripts/REV2_paths.R. Run 00_setup_REV2.R first.")
}
source(paths_file)

if (!requireNamespace("readr", quietly = TRUE) ||
    !requireNamespace("dplyr", quietly = TRUE) ||
    !requireNamespace("tibble", quietly = TRUE)) {
  stop("Required packages: readr, dplyr, tibble")
}

dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
fun_dir <- file.path(audit_dir, "02_model_functions")
dir.create(fun_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. Load REV1 workspace privately
# ------------------------------------------------------------

if (!file.exists(rev1_rdata_general)) {
  stop("REV1 RData not found: ", rev1_rdata_general)
}

rev1 <- new.env(parent = baseenv())
load(rev1_rdata_general, envir = rev1)

get_rev1 <- function(nm, required = FALSE) {
  if (!exists(nm, envir = rev1, inherits = FALSE)) {
    if (required) stop("Required REV1 object not found: ", nm)
    return(NULL)
  }
  get(nm, envir = rev1, inherits = FALSE)
}

# ------------------------------------------------------------
# 3. Model-table inventory
# ------------------------------------------------------------

object_names <- c(
  "bg",
  "bg_revision",
  "raw",
  "full_df",
  "model_df",
  "train_df",
  "tr",
  "te",
  "test_df"
)

inventory <- lapply(object_names, function(nm) {
  x <- get_rev1(nm)

  tibble::tibble(
    object = nm,
    exists = !is.null(x),
    class = if (is.null(x)) NA_character_ else paste(class(x), collapse = "; "),
    nrow = if (is.data.frame(x) || is.matrix(x)) nrow(x) else NA_integer_,
    ncol = if (is.data.frame(x) || is.matrix(x)) ncol(x) else NA_integer_,
    fields = if (is.null(x)) NA_character_ else paste(names(x), collapse = "; ")
  )
}) |>
  dplyr::bind_rows()

readr::write_csv(
  inventory,
  file.path(audit_dir, "02A_model_object_inventory.csv")
)

# ------------------------------------------------------------
# 4. pa / class counts
# ------------------------------------------------------------

pa_rows <- list()

for (nm in object_names) {
  x <- get_rev1(nm)

  if (is.null(x) || !is.data.frame(x) || !"pa" %in% names(x)) next

  tab <- table(x$pa, useNA = "ifany")

  pa_rows[[nm]] <- tibble::tibble(
    object = nm,
    pa_value = names(tab),
    n = as.integer(tab)
  )
}

pa_counts <- if (length(pa_rows) > 0) {
  dplyr::bind_rows(pa_rows)
} else {
  tibble::tibble(
    object = character(),
    pa_value = character(),
    n = integer()
  )
}

readr::write_csv(
  pa_counts,
  file.path(audit_dir, "02B_pa_class_counts.csv")
)

# ------------------------------------------------------------
# 5. Search model tables for explicit weight fields
# ------------------------------------------------------------

weight_rows <- list()

for (nm in object_names) {
  x <- get_rev1(nm)

  if (is.null(x) || !is.data.frame(x)) next

  hits <- grep(
    "weight|case_weight|class_weight|sample_weight",
    names(x),
    ignore.case = TRUE,
    value = TRUE
  )

  weight_rows[[nm]] <- tibble::tibble(
    object = nm,
    explicit_weight_field = if (length(hits) == 0) NA_character_ else paste(hits, collapse = "; ")
  )
}

weight_fields <- if (length(weight_rows) > 0) {
  dplyr::bind_rows(weight_rows)
} else {
  tibble::tibble(object = character(), explicit_weight_field = character())
}

readr::write_csv(
  weight_fields,
  file.path(audit_dir, "02C_weight_fields_in_model_tables.csv")
)

# ------------------------------------------------------------
# 6. Export actual REV1 model functions
# ------------------------------------------------------------

function_names <- c(
  "prepare_revision_model_table",
  "fit_glm_revision",
  "fit_gam_revision",
  "fit_rf_revision",
  "fit_maxnet_revision",
  "fit_all_revision_models",
  "predict_revision_model",
  "predict_all_revision_models",
  "evaluate_fold_models"
)

function_audit <- list()

for (fn in function_names) {
  f <- get_rev1(fn)

  if (is.function(f)) {
    txt <- paste(deparse(f, width.cutoff = 500L), collapse = "\n")

    writeLines(
      txt,
      file.path(fun_dir, paste0(fn, ".R.txt"))
    )

    function_audit[[fn]] <- tibble::tibble(
      function_name = fn,
      found = TRUE,
      contains_pa = grepl("\\bpa\\b", txt),
      contains_weights_argument = grepl(
        "weights\\s*=|case\\.weights|class\\.weights|sample\\.fraction",
        txt,
        ignore.case = TRUE
      ),
      contains_binomial = grepl("binomial", txt, ignore.case = TRUE),
      contains_ranger = grepl("ranger", txt, ignore.case = TRUE),
      contains_maxnet = grepl("maxnet", txt, ignore.case = TRUE)
    )
  } else {
    function_audit[[fn]] <- tibble::tibble(
      function_name = fn,
      found = FALSE,
      contains_pa = NA,
      contains_weights_argument = NA,
      contains_binomial = NA,
      contains_ranger = NA,
      contains_maxnet = NA
    )
  }
}

function_audit <- dplyr::bind_rows(function_audit)

readr::write_csv(
  function_audit,
  file.path(audit_dir, "02D_model_function_audit.csv")
)

# ------------------------------------------------------------
# 7. Export REV1 settings and parameter trace
# ------------------------------------------------------------

settings <- get_rev1("revision_settings")
if (!is.null(settings)) {
  capture.output(
    str(settings, max.level = 6),
    file = file.path(audit_dir, "02E_revision_settings_structure.txt")
  )
  capture.output(
    print(settings),
    file = file.path(audit_dir, "02F_revision_settings_values.txt")
  )
}

parameter_trace <- get_rev1("parameter_trace")
if (is.data.frame(parameter_trace)) {
  readr::write_csv(
    parameter_trace,
    file.path(audit_dir, "02G_parameter_trace.csv")
  )
}

balance_tbl <- get_rev1("balance_tbl")
if (is.data.frame(balance_tbl)) {
  readr::write_csv(
    balance_tbl,
    file.path(audit_dir, "02H_REV1_balance_table.csv")
  )
}

# ------------------------------------------------------------
# 8. Determine what the REV1 model matrix actually represents
# ------------------------------------------------------------

bg <- get_rev1("bg_revision")
if (is.null(bg)) bg <- get_rev1("bg")

full_df <- get_rev1("full_df")
model_df <- get_rev1("model_df")

bg_zero_encoded <- FALSE
if (is.data.frame(bg) && "pa" %in% names(bg)) {
  vals <- unique(stats::na.omit(bg$pa))
  bg_zero_encoded <- length(vals) > 0 && all(as.character(vals) == "0")
}

model_has_both_classes <- FALSE
if (is.data.frame(full_df) && "pa" %in% names(full_df)) {
  vals <- sort(unique(stats::na.omit(full_df$pa)))
  model_has_both_classes <- all(c(0, 1) %in% vals)
}

any_explicit_weight_field <- any(
  !is.na(weight_fields$explicit_weight_field) &
  nzchar(weight_fields$explicit_weight_field)
)

fit_rows <- function_audit |>
  dplyr::filter(
    function_name %in% c(
      "fit_glm_revision",
      "fit_gam_revision",
      "fit_rf_revision",
      "fit_maxnet_revision"
    ),
    found
  )

any_fit_weighting <- any(
  fit_rows$contains_weights_argument %in% TRUE,
  na.rm = TRUE
)

diagnosis <- if (bg_zero_encoded && model_has_both_classes) {
  paste(
    "REV1 background records are encoded as pa = 0 and presences as pa = 1.",
    "For Maxnet this is a standard presence-background formulation.",
    "For GLM, GAM and random forest, the same background rows operate as the zero class",
    "(i.e. background/pseudo-absence surrogates) unless explicit case/class weighting is applied."
  )
} else {
  paste(
    "The saved REV1 objects do not support an automatic diagnosis of a 1/0 model matrix.",
    "Inspect the exported model functions and pa-count table."
  )
}

weight_diagnosis <- if (!any_explicit_weight_field && !any_fit_weighting) {
  paste(
    "No explicit observation/class weight field or weighting argument was detected",
    "in the saved model tables/functions audited here."
  )
} else {
  paste(
    "At least one explicit weighting field or weighting argument was detected.",
    "Inspect 02C and 02D before interpreting model balance."
  )
}

diagnostic_summary <- tibble::tibble(
  item = c(
    "Background encoded exclusively as pa=0",
    "Saved full modelling table contains pa=0 and pa=1",
    "Explicit weight fields found in model tables",
    "Explicit weighting arguments detected in fitting functions",
    "Interpretation",
    "Weighting interpretation"
  ),
  value = c(
    as.character(bg_zero_encoded),
    as.character(model_has_both_classes),
    as.character(any_explicit_weight_field),
    as.character(any_fit_weighting),
    diagnosis,
    weight_diagnosis
  )
)

readr::write_csv(
  diagnostic_summary,
  file.path(audit_dir, "02_model_input_diagnosis.csv")
)

# ------------------------------------------------------------
# 9. Report
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("REV2 MODEL-INPUT AUDIT COMPLETE\n")
cat("============================================================\n\n")

cat("Object inventory:\n")
print(inventory, n = Inf)

cat("\npa class counts:\n")
print(pa_counts, n = Inf)

cat("\nModel-function audit:\n")
print(function_audit, n = Inf)

cat("\nDIAGNOSIS:\n")
cat(diagnosis, "\n\n")
cat(weight_diagnosis, "\n")

cat("\nNo model has been refitted and no REV1 file has been modified.\n")
cat("============================================================\n")
