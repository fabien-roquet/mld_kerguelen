#!/usr/bin/env python3
"""Python script version of Figure_8_trend_maps.ipynb."""

from __future__ import annotations

import argparse

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from figure_common import add_common_map_layers, cmo, load_fpca, parse_project_root_arg, paths, save_figure, topo_fronts


def seasonal_trend_map(ds_in, months, min_years=8):
    da = ds_in["mld"]
    sub = da.where(da["time"].dt.month.isin(months), drop=True)
    if set(months) == {12, 1, 2}:
        season_year = xr.where(sub["time"].dt.month == 12, sub["time"].dt.year + 1, sub["time"].dt.year)
    else:
        season_year = sub["time"].dt.year
    da_y = sub.groupby(season_year.rename("season_year")).mean("time", skipna=True)
    x = da_y["season_year"].astype(float)

    def _slope(y, x_values, min_points=8):
        valid = np.isfinite(y) & np.isfinite(x_values)
        if valid.sum() < min_points:
            return np.nan
        return np.polyfit(x_values[valid], y[valid], 1)[0]

    return xr.apply_ufunc(
        _slope,
        da_y,
        x,
        input_core_dims=[["season_year"], ["season_year"]],
        output_core_dims=[[]],
        vectorize=True,
        dask="allowed",
        kwargs={"min_points": min_years},
        output_dtypes=[float],
    )


def main() -> None:
    parser = parse_project_root_arg(argparse.ArgumentParser(description=__doc__))
    args = parser.parse_args()
    plt.rcParams.update({"font.size": 20})

    elevation, ds_front = topo_fronts(args.project_root)
    fpca = load_fpca(args.project_root)
    ds_map = {"GLORYS": fpca["GLORYS"][0], "GLORYS_CL": fpca["GLORYS_CL"][0], "CMA": fpca["CMA"][0]}
    periods = {"Annual": list(range(1, 13)), "Summer (JFM)": [1, 2, 3], "Winter (JAS)": [7, 8, 9]}
    dataset_order = ["GLORYS", "GLORYS_CL", "CMA"]
    trend_maps_period_dataset = {
        p_name: {d_name: seasonal_trend_map(ds_map[d_name], months) for d_name in dataset_order}
        for p_name, months in periods.items()
    }

    fig, axes = plt.subplots(
        nrows=len(periods),
        ncols=len(dataset_order),
        figsize=(18, 15),
        sharex=True,
        sharey=True,
        constrained_layout=True,
        squeeze=False,
    )
    cbar_ticks_uniform = np.arange(-2.0, 2.5, 1)
    for i, period_name in enumerate(periods.keys()):
        for j, d_name in enumerate(dataset_order):
            ax = axes[i, j]
            tr = trend_maps_period_dataset[period_name][d_name]
            im = ax.pcolormesh(tr["long"], tr["lat"], tr.transpose("lat", "long"), shading="auto", cmap=cmo.balance, vmin=-2, vmax=2)
            add_common_map_layers(ax, elevation, ds_front)
            ax.set_xlabel("Longitude [deg E]" if i == len(periods) - 1 else "")
            ax.set_ylabel("Latitude [deg N]" if j == 0 else "")
        cbar = fig.colorbar(im, ax=axes[i, :], orientation="vertical", shrink=0.95, pad=0.02, ticks=cbar_ticks_uniform)
        cbar.set_label("MLD trend [m yr$^{-1}$]")

    for i, ax_i in enumerate(axes.ravel()):
        period_idx = i // len(dataset_order)
        dataset_idx = i % len(dataset_order)
        d_name = dataset_order[dataset_idx]
        label = r"GLORYS$_{\mathregular{CL}}$" if d_name == "GLORYS_CL" else d_name
        ax_i.text(0.01, 0.98, f"({chr(97 + i)}) {label} - {list(periods.keys())[period_idx]}", transform=ax_i.transAxes, ha="left", va="top", fontsize=20, fontweight="bold")

    save_figure(fig, paths(args.project_root)["figures"] / "Figure_8_2D_trends.png")


if __name__ == "__main__":
    main()
