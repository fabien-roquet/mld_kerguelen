# MLD Kerguelen

Reproducible analysis workflow for Kerguelen mixed-layer-depth variability.

The project is organized as a three-stage pipeline:

1. `1_data_processing/`: process the raw CMA and GLORYS NetCDF inputs into gridded products and R-ready tables.
2. `2_compute_fPCA_R/`: run the functional PCA analysis in R.
3. `3_generate_figures/`: build figures from processed gridded products and fPCA outputs.

Large NetCDF inputs and generated outputs are intentionally ignored by Git. Place the required local inputs under `data/` before running the pipeline.
