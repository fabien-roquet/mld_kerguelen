#!/usr/bin/env python3
"""Python script version of Figure_7_1D_trends.ipynb."""

from __future__ import annotations

import argparse

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from scipy.stats import linregress

from figure_common import CL_COLOR, CMA_COLOR, G_COLOR, LEGEND_FS, load_fpca, parse_project_root_arg, paths, save_figure


def seasonal_series_full_region(ds_in, months):
    ts = ds_in["mld"].mean(dim=("long", "lat"), skipna=True)
    sub = ts.where(ts["time"].dt.month.isin(months), drop=True)
    if set(months) == {12, 1, 2}:
        season_year = xr.where(sub["time"].dt.month == 12, sub["time"].dt.year + 1, sub["time"].dt.year)
    else:
        season_year = sub["time"].dt.year
    return sub.groupby(season_year.rename("season_year")).mean("time", skipna=True)


def main() -> None:
    parser = parse_project_root_arg(argparse.ArgumentParser(description=__doc__))
    args = parser.parse_args()
    plt.rcParams.update({"font.size": 20})

    fpca = load_fpca(args.project_root)
    season_months = {"Annual": list(range(1, 13)), "Summer (JFM)": [1, 2, 3], "Winter (JAS)": [7, 8, 9]}
    ds_map = {"GLORYS": fpca["GLORYS"][0], "GLORYS_CL": fpca["GLORYS_CL"][0], "CMA": fpca["CMA"][0]}
    colors = {"GLORYS": G_COLOR, "GLORYS_CL": CL_COLOR, "CMA": CMA_COLOR}

    fig, axes = plt.subplots(3, 1, figsize=(15, 15), sharex=True, tight_layout=True, gridspec_kw={"hspace": 0.05})
    for ax, (sname, months), letter in zip(axes, season_months.items(), ["a", "b", "c"]):
        for name, ds_in in ds_map.items():
            ts_y = seasonal_series_full_region(ds_in, months)
            x = ts_y["season_year"].values.astype(float)
            y = ts_y.values
            valid = np.isfinite(x) & np.isfinite(y)
            ax.plot(x[valid], y[valid], marker="o", ms=4, lw=1.3, color=colors[name], alpha=0.5, label="_nolegend_")
            if valid.sum() >= 8:
                lr = linregress(x[valid], y[valid])
                yhat = lr.intercept + lr.slope * x[valid]
                ax.plot(x[valid], yhat, lw=2.2, color=colors[name], label=f"{lr.slope:.2f} m yr$^{{-1}}$ (p={lr.pvalue:.2f})")
        ax.set_ylabel("MLD [m]")
        ax.axhline(0, color="k", linestyle="--", linewidth=1)
        ax.grid(alpha=0.3)
        ax.text(0.01, 0.98, f"({letter}) {sname}", transform=ax.transAxes, ha="left", va="top", fontsize=20, fontweight="bold")

    axes[0].set_ylim(-8, 8)
    axes[0].legend(fontsize=LEGEND_FS, loc="upper right", ncol=1)
    axes[1].set_ylim(-15, 15)
    axes[1].legend(fontsize=LEGEND_FS, loc="lower right", ncol=1)
    axes[2].set_ylim(-30, 30)
    axes[2].legend(fontsize=LEGEND_FS, loc="lower left", ncol=1)
    axes[-1].set_xlabel("Time")
    fig.align_ylabels()
    save_figure(fig, paths(args.project_root)["figures"] / "Figure_7_1D_trends.png")


if __name__ == "__main__":
    main()
