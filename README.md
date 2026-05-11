# MLD Kerguelen

Reproducible analysis workflow for Kerguelen mixed-layer-depth variability.

The code is organized as the three-step analysis used in the notebooks:

1. `1_data_processing/`: one Python script per processed dataset.
2. `2_compute_fPCA_R/`: one R script per fPCA dataset.
3. `3_generate_figures/`: one Python script per figure notebook.

Large NetCDF inputs and generated outputs are intentionally ignored by Git. Place the required local inputs under `data/` before running the pipeline.

## Inputs

Expected local files:

- `data/CORA_MEOP_ARGO_2026.nc`
- `data/GLORYS_2026.nc`
- `data/GEBCO_ker_large.nc`
- `data/fronts_62985.nc`
- `data/GLORYS_1000m_section_timemean.nc` for Figure 9
- `data/kerfix.csv` for Figure 10

The Figure 9 time-mean section file is generated once from the large temporary
source `data/GLORYS_1000m_section_plot.nc`:

```bash
uv run python 1_data_processing/process_GLORYS_section.py
```

After `data/GLORYS_1000m_section_timemean.nc` exists, the large source file is
not required by the figure script.

## Running

Install Python dependencies:

```bash
uv sync
```

Run the whole modular pipeline:

```bash
uv run python run_analysis.py
```

Run only one stage:

```bash
uv run python run_analysis.py --stage data
uv run python run_analysis.py --stage fpca
uv run python run_analysis.py --stage figures --figures 2 3 4
```

Run individual scripts directly:

```bash
uv run python 1_data_processing/process_CMA.py
Rscript 2_compute_fPCA_R/script_PCA_CM_2026.R
uv run python 3_generate_figures/Figure_4_MU_modes.py
```

## Reference Outputs

`processed_reference/` is a local, ignored snapshot of the current `processed/` folder. Use it to check whether regenerated outputs still mostly match the reference:

```bash
uv run python scripts/compare_processed_reference.py
```

Figure 1 is represented by a script, but it requires source files that are not part of the reduced two-input workflow.
