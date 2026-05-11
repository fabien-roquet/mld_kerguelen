#!/usr/bin/env python3
"""Build the gridded MLD products and R-input tables for the 2026 analysis."""

from __future__ import annotations

import argparse
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from pandas.tseries.offsets import DateOffset


NB_BINS = 39
END_DATE = pd.Timestamp("2023-12-31")
KERGUELEN_BOX = {
    "lon_min": 68.25,
    "lon_max": 70.75,
    "lat_min": -50.0,
    "lat_max": -48.0,
}


@dataclass(frozen=True)
class ProjectPaths:
    root: Path
    data: Path
    gridded: Path
    r_input: Path


def project_paths(project_root: Path) -> ProjectPaths:
    root = project_root.resolve()
    gridded = root / "processed" / "1_gridded_data"
    paths = ProjectPaths(
        root=root,
        data=root / "data",
        gridded=gridded,
        r_input=gridded / "r_input",
    )
    paths.gridded.mkdir(parents=True, exist_ok=True)
    paths.r_input.mkdir(parents=True, exist_ok=True)
    return paths


def rename_obs_coords(ds: xr.Dataset) -> xr.Dataset:
    rename = {}
    if "LATITUDE" in ds:
        rename["LATITUDE"] = "latitude"
    if "LONGITUDE" in ds:
        rename["LONGITUDE"] = "longitude"
    return ds.rename(rename) if rename else ds


def make_ds_cut(
    df: pd.DataFrame,
    nb_bins: int = NB_BINS,
    start_date: pd.Timestamp | None = None,
    end_date: pd.Timestamp = END_DATE,
    freq: str = "ME",
) -> tuple[xr.Dataset, pd.DataFrame]:
    """Bin point MLD values onto the 39 x 39 monthly grid used by the notebooks."""
    lon_col = "longitude"
    lat_col = "latitude"
    time_col = "time"

    bins_dt = pd.date_range(
        start=start_date if start_date is not None else df[time_col].min() + DateOffset(months=-1),
        end=end_date + DateOffset(months=1),
        freq=freq,
    )

    cut_lat_label = pd.cut(df[lat_col], nb_bins)
    cut_lon_label = pd.cut(df[lon_col], nb_bins)
    cut_time_label = pd.cut(df[time_col], bins=bins_dt)

    df_cut_label = df.drop([lat_col, lon_col, time_col], axis=1)
    df_cut_label = df_cut_label.groupby(
        [cut_time_label, cut_lon_label, cut_lat_label],
        observed=False,
    ).mean()

    lat_mid = pd.IntervalIndex(df_cut_label.index.get_level_values(2)).mid.unique()
    lon_mid = pd.IntervalIndex(df_cut_label.index.get_level_values(1)).mid.unique()
    time_mid = pd.IntervalIndex(df_cut_label.index.get_level_values(0)).mid.unique()

    df_cut_label.index = df_cut_label.index.set_levels(time_mid.values, level=0)
    df_cut_label.index = df_cut_label.index.set_levels(lon_mid, level=1)
    df_cut_label.index = df_cut_label.index.set_levels(lat_mid.values, level=2)

    df_cut = df_cut_label.copy()
    df_cut.replace(0, np.nan, inplace=True)
    df_cut = df_cut.drop(columns=["uid"], errors="ignore")

    ds_cut = df_cut.to_xarray()
    ds_cut["latitude"] = sorted(lat_mid)
    ds_cut["longitude"] = sorted(lon_mid)
    ds_cut["time"] = time_mid
    return ds_cut, df_cut


def kerguelen_mask(ds: xr.Dataset, lon_name: str = "longitude", lat_name: str = "latitude") -> xr.DataArray:
    return (
        (ds[lon_name] >= KERGUELEN_BOX["lon_min"])
        & (ds[lon_name] <= KERGUELEN_BOX["lon_max"])
        & (ds[lat_name] >= KERGUELEN_BOX["lat_min"])
        & (ds[lat_name] <= KERGUELEN_BOX["lat_max"])
    )


def write_r_input(ds: xr.Dataset, out_file: Path) -> None:
    df = ds.where(ds.time.dt.year < 2024, drop=True).to_dataframe().reset_index()
    df["year"] = df.time.dt.year
    df["mth"] = df.time.dt.month
    df = (
        df.drop(columns="time")
        .rename(columns={"latitude": "lat", "longitude": "long"})
        [["long", "lat", "mth", "year", "mld"]]
        .sort_values(["year", "mth", "long", "lat"])
        .reset_index(drop=True)
    )
    out_file.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_file, sep=" ", index=False, na_rep="NA")
    print(f"Wrote {out_file.relative_to(out_file.parents[3]) if len(out_file.parents) > 3 else out_file}")


def write_anomaly_products(ds: xr.Dataset, mask: xr.DataArray, anom_file: Path, clim_file: Path) -> xr.Dataset:
    masked = ds.where(~mask)
    clim = masked.groupby("time.month").mean("time").mean(["latitude", "longitude"])
    anom = masked.groupby("time.month") - clim
    anom = anom.where(~mask, other=0)
    clim = clim.where(~mask, other=0)
    anom.to_netcdf(anom_file)
    clim.to_netcdf(clim_file)
    return anom


def build_observation_grid(paths: ProjectPaths, nb_bins: int) -> xr.Dataset:
    ds_obs = xr.open_dataset(paths.data / "CORA_MEOP_ARGO_2026.nc")
    ds_obs = rename_obs_coords(ds_obs)
    df_obs = ds_obs.mld.to_dataframe().reset_index()
    unmasked, _ = make_ds_cut(df_obs, nb_bins=nb_bins)
    return unmasked


def load_observation_grid(paths: ProjectPaths, force: bool, nb_bins: int) -> tuple[xr.Dataset, xr.Dataset]:
    """Return unmasked and Kerguelen-masked CMA grids."""
    cma_file = paths.gridded / "CMA_gridded.nc"
    anom_file = paths.gridded / "CMA_anom.nc"
    clim_file = paths.gridded / "CMA_clim.nc"

    if not force and cma_file.exists() and anom_file.exists() and clim_file.exists():
        masked = xr.open_dataset(cma_file)
        write_r_input(masked, paths.r_input / "CMA_masked.txt")
        return masked, masked

    unmasked = build_observation_grid(paths, nb_bins)
    mask = kerguelen_mask(unmasked)
    masked = unmasked.where(~mask)
    write_anomaly_products(masked, mask, anom_file, clim_file)
    masked.to_netcdf(cma_file)
    write_r_input(masked, paths.r_input / "CMA_masked.txt")
    return unmasked, masked


def sigma0_unesco(salinity: xr.DataArray, temperature: xr.DataArray) -> xr.DataArray:
    """EOS-80 density anomaly at atmospheric pressure, used if gsw is unavailable."""
    t = temperature
    s = salinity
    rho_w = (
        999.842594
        + 6.793952e-2 * t
        - 9.095290e-3 * t**2
        + 1.001685e-4 * t**3
        - 1.120083e-6 * t**4
        + 6.536332e-9 * t**5
    )
    a = 0.824493 - 4.0899e-3 * t + 7.6438e-5 * t**2 - 8.2467e-7 * t**3 + 5.3875e-9 * t**4
    b = -5.72466e-3 + 1.0227e-4 * t - 1.6546e-6 * t**2
    c = 4.8314e-4
    return (rho_w + a * s + b * s**1.5 + c * s**2 - 1000.0).rename("sigma0")


def compute_sigma0(ds: xr.Dataset) -> xr.DataArray:
    try:
        import gsw  # type: ignore
    except ModuleNotFoundError:
        warnings.warn(
            "gsw is not installed; using an EOS-80 sigma0 approximation for GLORYS MLD.",
            RuntimeWarning,
            stacklevel=2,
        )
        return sigma0_unesco(ds["so"], ds["thetao"])

    return xr.apply_ufunc(
        gsw.sigma0,
        ds["so"],
        ds["thetao"],
        dask="allowed",
        output_dtypes=[float],
    ).rename("sigma0")


def mld_from_delta(delta: xr.DataArray, depth: xr.DataArray, threshold: float) -> xr.DataArray:
    hit = delta >= threshold
    depth_index = hit.argmax("depth")
    has_mld = hit.any("depth")
    mld = xr.where(has_mld, depth.isel(depth=depth_index), np.nan)
    mld = mld.where(mld < depth.max())
    mld.name = "mld"
    mld.attrs["standard_name"] = "ocean_mixed_layer_thickness"
    mld.attrs["units"] = depth.attrs.get("units", "m")
    return mld


def compute_glorys_mlds(
    ds: xr.Dataset,
    thresholds: dict[str, float],
    z_ref: float = 10.0,
    block_size: int = 12,
) -> dict[str, xr.DataArray]:
    """Compute GLORYS MLD for one or more density thresholds in time blocks."""
    chunks: dict[str, list[xr.DataArray]] = {name: [] for name in thresholds}
    depth = ds["depth"]

    for start in range(0, ds.sizes["time"], block_size):
        stop = min(start + block_size, ds.sizes["time"])
        part = ds.isel(time=slice(start, stop))
        sigma0 = compute_sigma0(part)
        sigma0_ref = sigma0.interp(depth=z_ref)
        delta_sigma0 = sigma0 - sigma0_ref
        for name, threshold in thresholds.items():
            chunks[name].append(mld_from_delta(delta_sigma0, depth, threshold).load())
        print(f"Computed GLORYS MLD block {start + 1}-{stop} of {ds.sizes['time']}")

    return {name: xr.concat(parts, dim=ds["time"]) for name, parts in chunks.items()}


def grid_glorys_mld(mld: xr.DataArray, obs_grid: xr.Dataset, nb_bins: int) -> xr.Dataset:
    ds_mld = mld.to_dataset(name="mld")
    df = ds_mld.mld.to_dataframe().reset_index()
    gridded, _ = make_ds_cut(df, nb_bins=nb_bins)
    gridded = gridded.interp(
        latitude=obs_grid.latitude,
        longitude=obs_grid.longitude,
        time=obs_grid.time,
        method="nearest",
    )
    return gridded.sel(time=slice(obs_grid.time.min(), obs_grid.time.max()))


def process_existing_products(paths: ProjectPaths, selected: set[str]) -> set[str]:
    processed = set()
    if "cma" in selected and (paths.gridded / "CMA_gridded.nc").exists():
        write_r_input(xr.open_dataset(paths.gridded / "CMA_gridded.nc"), paths.r_input / "CMA_masked.txt")
        processed.add("cma")
    if "glorys" in selected and (paths.gridded / "GLORYS_anom.nc").exists():
        write_r_input(xr.open_dataset(paths.gridded / "GLORYS_anom.nc"), paths.r_input / "GLORYS_masked.txt")
        processed.add("glorys")
    if "glorys-cl" in selected and (paths.gridded / "GLORYS_CL_anom.nc").exists():
        write_r_input(xr.open_dataset(paths.gridded / "GLORYS_CL_anom.nc"), paths.r_input / "GLORYS_CL_masked.txt")
        processed.add("glorys-cl")
    return processed


def process_glorys_products(paths: ProjectPaths, selected: set[str], force: bool, nb_bins: int) -> None:
    required = {
        "glorys": [
            paths.gridded / "GLORYS_gridded.nc",
            paths.gridded / "GLORYS_anom.nc",
            paths.gridded / "GLORYS_clim.nc",
        ],
        "glorys-cl": [
            paths.gridded / "GLORYS_CL_gridded.nc",
            paths.gridded / "GLORYS_CL_anom.nc",
            paths.gridded / "GLORYS_CL_clim.nc",
        ],
    }
    needs_compute = {
        name
        for name in selected.intersection(required)
        if force or not all(path.exists() for path in required[name])
    }
    if not needs_compute:
        return

    if "glorys-cl" in needs_compute:
        obs_unmasked = build_observation_grid(paths, nb_bins)
    else:
        obs_unmasked, _ = load_observation_grid(paths, force=False, nb_bins=nb_bins)
    raw_obs = xr.open_dataset(paths.data / "CORA_MEOP_ARGO_2026.nc")
    raw_obs = rename_obs_coords(raw_obs)
    ds_glorys = xr.open_dataset(paths.data / "GLORYS_2026.nc")
    ds_glorys = ds_glorys.sel(
        latitude=slice(raw_obs.latitude.min() - 0.5, raw_obs.latitude.max() + 0.5),
        longitude=slice(raw_obs.longitude.min() - 0.5, raw_obs.longitude.max() + 0.5),
    )

    thresholds = {}
    if "glorys" in needs_compute:
        thresholds["glorys"] = 0.03
    if "glorys-cl" in needs_compute:
        thresholds["glorys-cl"] = 0.02

    mlds = compute_glorys_mlds(ds_glorys, thresholds)

    if "glorys" in needs_compute:
        ds_g = grid_glorys_mld(mlds["glorys"], obs_unmasked, nb_bins=nb_bins)
        ds_g.to_netcdf(paths.gridded / "GLORYS_gridded.nc")
        mask = kerguelen_mask(ds_g)
        anom = write_anomaly_products(ds_g, mask, paths.gridded / "GLORYS_anom.nc", paths.gridded / "GLORYS_clim.nc")
        write_r_input(anom, paths.r_input / "GLORYS_masked.txt")

    if "glorys-cl" in needs_compute:
        ds_cl = grid_glorys_mld(mlds["glorys-cl"], obs_unmasked, nb_bins=nb_bins)
        ds_cl["mld"] = ds_cl.mld.where(~np.isnan(obs_unmasked.mld))
        ds_cl.to_netcdf(paths.gridded / "GLORYS_CL_gridded.nc")
        mask = kerguelen_mask(ds_cl)
        anom = write_anomaly_products(ds_cl, mask, paths.gridded / "GLORYS_CL_anom.nc", paths.gridded / "GLORYS_CL_clim.nc")
        write_r_input(anom, paths.r_input / "GLORYS_CL_masked.txt")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="Repository root containing data/ and processed/.",
    )
    parser.add_argument(
        "--only",
        nargs="+",
        choices=["cma", "glorys", "glorys-cl"],
        default=["cma", "glorys", "glorys-cl"],
        help="Subset of products to build.",
    )
    parser.add_argument("--force", action="store_true", help="Recompute NetCDF products even when they exist.")
    parser.add_argument("--nb-bins", type=int, default=NB_BINS, help="Number of latitude/longitude bins.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    paths = project_paths(args.project_root)
    selected = set(args.only)

    reused = set() if args.force else process_existing_products(paths, selected)
    if "cma" in selected and "cma" not in reused:
        load_observation_grid(paths, force=True, nb_bins=args.nb_bins)
    process_glorys_products(paths, selected - reused, force=args.force, nb_bins=args.nb_bins)


if __name__ == "__main__":
    main()
