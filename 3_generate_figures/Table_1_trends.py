#!/usr/bin/env python3
"""Generate the manuscript trend table from the reconstructed products."""

from __future__ import annotations

import argparse

import numpy as np
import xarray as xr
from scipy.stats import linregress

from figure_common import load_fpca, parse_project_root_arg, paths


def seasonal_series_full_region(ds_in, months):
    ts = ds_in["mld"].mean(dim=("long", "lat"), skipna=True)
    sub = ts.where(ts["time"].dt.month.isin(months), drop=True)
    if set(months) == {12, 1, 2}:
        season_year = xr.where(sub["time"].dt.month == 12, sub["time"].dt.year + 1, sub["time"].dt.year)
    else:
        season_year = sub["time"].dt.year
    return sub.groupby(season_year.rename("season_year")).mean("time", skipna=True)


def trend_stats(ds_in, months) -> tuple[float, float]:
    ts_y = seasonal_series_full_region(ds_in, months)
    x = ts_y["season_year"].values.astype(float)
    y = ts_y.values
    valid = np.isfinite(x) & np.isfinite(y)
    if valid.sum() < 8:
        return np.nan, np.nan
    lr = linregress(x[valid], y[valid])
    return float(lr.slope), float(lr.pvalue)


def fmt_number(value: float) -> str:
    if not np.isfinite(value):
        return "--"
    text = f"{value:.2f}".rstrip("0").rstrip(".")
    return "0" if text == "-0" else text


def fmt_cell(slope: float, pvalue: float) -> str:
    return rf"${fmt_number(slope)}$ ({fmt_number(pvalue)})"


def main() -> None:
    parser = parse_project_root_arg(argparse.ArgumentParser(description=__doc__))
    args = parser.parse_args()

    fpca = load_fpca(args.project_root)
    periods = {"Annual": list(range(1, 13)), "Summer": [1, 2, 3], "Winter": [7, 8, 9]}
    dataset_order = ["GLORYS", "GLORYS_CL", "CMA"]
    ds_map = {"GLORYS": fpca["GLORYS"][0], "GLORYS_CL": fpca["GLORYS_CL"][0], "CMA": fpca["CMA"][0]}

    rows = []
    for period_name, months in periods.items():
        cells = [fmt_cell(*trend_stats(ds_map[dataset_name], months)) for dataset_name in dataset_order]
        rows.append(f"{period_name} & " + " & ".join(cells) + r" \\")

    table = "\n".join(
        [
            r"\begin{table}",
            r"\centering",
            r"\caption{Summary of MLD annual and seasonal trend slope coefficients (in $\mathrm{m\,yr^{-1}}$) with their associated $p$-value in parenthesis.}",
            r"\label{table2}",
            r"\small",
            r"\begin{tabular}{c|c|c|c}",
            r"\hline",
            r"Period & GLORYS & $\textrm{GLORYS}_{\textrm{CL}}$ & CMA \\",
            r"\hline",
            *rows,
            r"\hline",
            r"\end{tabular}",
            r"\end{table}",
            "",
        ]
    )

    out_file = paths(args.project_root)["figures"] / "Table_1_trends.tex"
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text(table)
    print(f"Wrote {out_file}")


if __name__ == "__main__":
    main()
