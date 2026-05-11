# Data Inputs

Expected local input files:

- `CORA_MEOP_ARGO_2026.nc`
- `GLORYS_2026.nc`
- `GEBCO_ker_large.nc`
- `fronts_62985.nc`
- `kerfix.csv` for the KERFIX figure workflow

The NetCDF files are ignored by Git because they are large. Keep them locally in this directory.

Some original notebooks referenced additional raw or section-specific files. Those figures are kept as explicit scripts, but they are not part of the default reproducible two-input pipeline unless the extra files are added here.
