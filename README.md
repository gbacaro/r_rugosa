# Rosa rugosa coastal climate workflow

This folder contains a cleaned, execution-ordered version of the pipeline used for the manuscript on the current and future climatic suitability of *Rosa rugosa* along European coasts. The workflow is organized so that the analytical core can be rerun from raw GBIF download to present and future projections, with a separate manuscript-output section for the final tables and figures.

## Execution order

Run all scripts from the repository root, in the order shown below.

### Core analytical workflow
1. `00_setup_project.R`  
   Creates the project directories and writes `project_settings.rds`.

2. `01_download_gbif_occurrences.R`  
   Downloads the GBIF archive and saves the raw occurrence table and metadata.

3. `02_clean_occurrences.R`  
   Cleans the raw occurrence table and exports the modelling dataset.

4. `03_build_geographic_coastal_domain.R`  
   Builds the 1-km geographic coastal domains for Europe and East Asia.

5. `04_download_chelsa_climate_coastal.R`  
   Downloads CHELSA present and future layers, crops them remotely, and masks them to the geographic coastal domains.

6. `05_define_raster_valid_coastal_occurrences_and_background.R`  
   Classifies occurrences as native or invaded, snaps records to raster-valid coastal cells, and samples the European coastal background.

7. `06_run_pca_env_niche_analysis.R`  
   Runs the native-versus-invaded coastal niche comparison.

8. `07_fit_current_sdm_ensemble.R`  
   Fits the present-day coastal SDM ensemble and exports the current ensemble products.

9. `08_project_future_scenarios.R`  
   Projects the ensemble to future CHELSA CMIP6 scenarios and exports scenario summaries.

10. `09_calculate_threshold_area_change.R`  
    Calculates current and future area above the mean cross-validated threshold.

### Manuscript-output workflow
11. `10_build_paper_outputs_revision1.R`  
    Rebuilds the paper tables and the main/supplementary figure set in `outputs/paper_revision1`.

12. `11_refresh_figure_2_climatic_space.R`  
    Rebuilds the final main-text Figure 2.

13. `12_refresh_main_figures_3_5_final.R`  
    Rebuilds the final main-text Figures 3–5.

## Main requirements

The workflow assumes:
- R with the packages used inside each script (`R_packages.txt` lists the package set referenced across the workflow)
- internet access for the GBIF and CHELSA download steps
- command-line GDAL available in `PATH`, especially `gdal_translate`
- GBIF credentials available as environment variables:
  - `GBIF_USER`
  - `GBIF_PWD`
  - `GBIF_EMAIL`

## Main output folders

After execution, the main outputs are written to:
- `data_raw/gbif`
- `data_raw/climate_chelsa`
- `data_clean`
- `models`
- `outputs/tables`
- `outputs/figures`
- `outputs/paper_revision1`
