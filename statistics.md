# Statistics Audit

This file collects the main numbers currently used to update `latex/manuscript_doriot_v3.tex`. Values are from the local `data/` and `processed/` outputs unless stated otherwise.

## Source Data and Gridding

| Quantity | Value | Source |
|---|---:|---|
| CMA source MLD profiles | 97,815 | `data/CORA_MEOP_ARGO_2026.nc` |
| Time span | Jan 2007-Dec 2023 | source data and processed grids |
| Monthly fields | 204 | `processed/1_gridded_data/*.nc` |
| Horizontal grid | 39 x 39 | processing scripts |
| Grid cells | 1,521 | processing scripts |
| Total monthly grid-cell slots | 310,284 | 204 x 39 x 39 |
| Kerguelen island mask | 20 cells | lon 68.25-70.75 E, lat 50-48 S |
| CMA finite monthly grid values | 27,149 | `CMA_gridded.nc` |
| CMA gridded coverage | 8.75% | 27,149 / 310,284 |
| CMA R-input non-missing values | 31,229 | includes island-mask zeros |
| CMA R-input coverage | 10.06% | 31,229 / 310,284 |

The merged CMA file does not contain a source/platform label, so the MEOP/CORA-Argo split cannot be recomputed directly from this reduced input file.

Profile counts by year in `data/CORA_MEOP_ARGO_2026.nc`:

| Year | Profiles |
|---:|---:|
| 2007 | 616 |
| 2008 | 2,323 |
| 2009 | 3,246 |
| 2010 | 2,335 |
| 2011 | 6,629 |
| 2012 | 6,710 |
| 2013 | 8,293 |
| 2014 | 11,187 |
| 2015 | 6,261 |
| 2016 | 3,982 |
| 2017 | 6,415 |
| 2018 | 5,370 |
| 2019 | 10,758 |
| 2020 | 5,866 |
| 2021 | 7,302 |
| 2022 | 7,075 |
| 2023 | 3,447 |

The average over 2011-2021 is 7,161 profiles per year, rounded to about 7,160 in the manuscript.

Profile counts by calendar month:

| Month | Profiles |
|---:|---:|
| 1 | 14,000 |
| 2 | 11,801 |
| 3 | 8,817 |
| 4 | 7,150 |
| 5 | 8,070 |
| 6 | 7,512 |
| 7 | 7,386 |
| 8 | 6,789 |
| 9 | 7,463 |
| 10 | 6,032 |
| 11 | 6,547 |
| 12 | 6,248 |

## Gridded MLD Statistics

Mean and standard deviation are over all finite monthly grid cells.

| Product | Valid values | Mean MLD (m) | Std MLD (m) |
|---|---:|---:|---:|
| GLORYS | 307,470 | 98.45 | 60.43 |
| GLORYS_CL | 27,411 | 93.99 | 54.55 |
| CMA | 27,149 | 109.21 | 68.94 |

Pointwise GLORYS-CMA difference at observed monthly grid cells:

| Difference | n | Mean (m) | Std (m) | RMSE (m) |
|---|---:|---:|---:|---:|
| GLORYS - CMA | 27,041 | -13.67 | 40.47 | 42.72 |
| GLORYS_CL - CMA | 27,041 | -13.67 | 40.47 | 42.72 |

Pointwise correlations at observed monthly grid cells:

| Pair | Correlation | n |
|---|---:|---:|
| GLORYS vs CMA | 0.786 | 27,041 |
| GLORYS_CL vs CMA | 0.786 | 27,041 |
| GLORYS vs GLORYS_CL | 1.000 | 27,411 |

Spatial climatology correlations:

| Pair | Correlation | n |
|---|---:|---:|
| GLORYS vs CMA | 0.668 | 1,501 |
| GLORYS_CL vs CMA | 0.840 | 1,501 |
| GLORYS vs GLORYS_CL | 0.741 | 1,521 |

## Seasonal Cycle and Domain-Mean Anomalies

Statistics below use the products plotted in Figure 3: monthly climatology files and domain-mean anomaly files.

| Product | Climatology min (m) | Climatology max (m) | Seasonal amplitude (m) | Anomaly mean (m) | Anomaly std (m) |
|---|---:|---:|---:|---:|---:|
| GLORYS | 49.12 | 176.23 | 127.11 | -0.96 | 8.57 |
| GLORYS_CL | 46.82 | 162.67 | 115.86 | -0.04 | 10.23 |
| CMA | 61.71 | 182.19 | 120.47 | -0.49 | 12.77 |

Domain-mean anomaly differences:

| Difference | Mean (m) | Std (m) | Min (m) | Max (m) |
|---|---:|---:|---:|---:|
| GLORYS_CL - GLORYS | 0.92 | 7.68 | -24.57 | 23.44 |
| CMA - GLORYS | 0.47 | 11.01 | -37.23 | 68.71 |

## fPCA/PACE Statistics

| Product | Retained modes | Mode 1 (%) | Mode 2 (%) | First 2 modes (%) | First 5 modes (%) |
|---|---:|---:|---:|---:|---:|
| GLORYS | 125 | 52.78 | 14.81 | 67.58 | 75.53 |
| GLORYS_CL | 86 | 48.44 | 10.62 | 59.06 | 76.64 |
| CMA | 88 | 42.75 | 12.08 | 54.84 | 72.36 |

Temporal-mode summary:

| Product | mean xi1 | std xi1 | mean xi2 | std xi2 |
|---|---:|---:|---:|---:|
| GLORYS | 0.179 | 0.174 | -0.139 | 0.209 |
| GLORYS_CL | 0.217 | 0.124 | -0.026 | 0.249 |
| CMA | 0.207 | 0.141 | -0.087 | 0.235 |

## Reconstruction Error

Figure 6 RMSE values are evaluated on the co-located `GLORYS_CL` sampling mask.

| Product | Mean RMSE (m) | Std RMSE (m) |
|---|---:|---:|
| GLORYS | 9.32 | 7.53 |
| GLORYS_CL | 16.11 | 8.48 |
| CMA | 29.11 | 17.59 |

## Domain-Mean Trends

Slopes are in `m yr-1`, with p-values in parentheses.

| Period | GLORYS | GLORYS_CL | CMA |
|---|---:|---:|---:|
| Annual | -0.20 (0.21) | -0.39 (0.05) | -0.30 (0.14) |
| Summer (JFM) | 0.17 (0.33) | 0.23 (0.36) | 0.23 (0.43) |
| Winter (JAS) | 0.06 (0.85) | -0.89 (0.05) | -0.63 (0.37) |

## KERFIX / Figure 10

KERFIX station position used in Figure 10: 50.6667 S, 68.4167 E. The GLORYS grid point selected by nearest-neighbor interpolation is 50.5125 S, 68.4615 E.

KERFIX MLD is recomputed from `data/kerfix.csv` with the same `0.03 kg m-3` density threshold relative to 10 m used in Figure 10.

| Product | Period | n | Mean MLD (m) | Std MLD (m) | Min (m) | Max (m) |
|---|---|---:|---:|---:|---:|---:|
| KERFIX | Dec 1990-Dec 1994 | 45 | 139.62 | 53.60 | 49.56 | 247.69 |
| GLORYS at KERFIX | Jan 2007-Dec 2023 | 204 | 107.07 | 45.91 | 34.05 | 229.94 |

Seasonal means:

| Product | Annual mean (m) | Summer JFM (m) | Winter JAS (m) |
|---|---:|---:|---:|
| KERFIX | 139.62 | 92.51 | 192.40 |
| GLORYS at KERFIX | 107.07 | 64.11 | 159.59 |

Simple monthly linear fits at the station:

| Product | Slope (m yr-1) | p-value |
|---|---:|---:|
| KERFIX, 1990-1994 | 15.23 | 0.037 |
| GLORYS at KERFIX, 2007-2023 | 0.13 | 0.840 |

Using the `0.02 kg m-3` threshold reported by Park et al. (1998), the same `kerfix.csv` gives a KERFIX mean MLD of `124.43 +/- 55.40 m`, with January and August monthly means of `69.4 m` and `188.3 m`, close to the Park et al. abstract values of about `60 m` in January and `185 m` in August.

## Park et al. (1998) Context

ScienceDirect metadata/abstract for Park et al. (1998) states that monthly hydrographic data were collected at KERFIX from May 1991 to December 1994, at about 50 deg 40 min S, 68 deg 25 min E. MLD was determined using a `0.02 sigma_theta` density-difference criterion. The monthly mean MLD varied from about 60 m in January to 185 m in August. The study reported warm/low-salinity anomalies in 1992, cold/high-salinity anomalies in 1994, and significantly deeper mixed layers in the later years of the short record, with possible links to ENSO.

Bibliographic check: the DOI `10.1016/S0924-7963(98)00065-7` corresponds to Journal of Marine Systems volume 17, issues 1-4, pages 571-586, not pages 233-247.
