# Data Inputs

Expected local input files:

- `CORA_MEOP_ARGO_2026.nc`
- `GLORYS_2026.nc`
- `GEBCO_ker_large.nc`
- `fronts_62985.nc`
- `kerfix.csv` for the KERFIX figure workflow
- `GLORYS_1000m_section_timemean.nc` for Figure 9

The large raw NetCDF files are ignored by Git. Keep them locally in this directory.
`GLORYS_1000m_section_timemean.nc` is generated from `GLORYS_1000m_section_plot.nc`
with `1_data_processing/process_GLORYS_section.py`; the large source file is not
needed after the time-mean file has been created.

Some original notebooks referenced additional raw or section-specific files. Those
figures are kept as explicit scripts, but they are not part of the reduced
two-input pipeline unless the extra files are added here.
