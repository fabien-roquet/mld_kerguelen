#!/usr/bin/env python3
"""Build the compact GLORYS 1000 m time-mean product used by Figure 9.

The large ``data/GLORYS_1000m_section_plot.nc`` source is only needed to create
``data/GLORYS_1000m_section_timemean.nc``. Figure 9 reads the compact time-mean
file so the large source can be removed afterwards.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import xarray as xr

from processing_common import PROJECT_ROOT, mld_from_density_threshold, project_paths


SOURCE_NAME = "GLORYS_1000m_section_plot.nc"
TIMEMEAN_NAME = "GLORYS_1000m_section_timemean.nc"


def _accumulate_time_mean(
    total: xr.Dataset | None,
    count: xr.Dataset | None,
    block: xr.Dataset,
    variables: list[str],
) -> tuple[xr.Dataset, xr.Dataset]:
    block_vars = block[variables]
    block_total = block_vars.sum("time", skipna=True).load()
    block_count = block_vars.notnull().sum("time").load()
    if total is None or count is None:
        return block_total, block_count
    return total + block_total, count + block_count


def build_glorys_section_timemean(
    project_root: str | Path = PROJECT_ROOT,
    force: bool = False,
    block_size: int = 12,
    density_threshold: float = 0.03,
    z_ref: float = 10.0,
    end_date: str = "2023-12-31",
) -> Path:
    paths = project_paths(project_root)
    source_file = paths.data / SOURCE_NAME
    out_file = paths.data / TIMEMEAN_NAME
    if out_file.exists() and not force:
        print(f"Reused existing {out_file}")
        return out_file
    if not source_file.exists():
        raise FileNotFoundError(
            f"Missing {source_file}. It is only needed to create {out_file}; "
            "after that, Figure 9 uses the compact time-mean file."
        )

    ds = xr.open_dataset(source_file).sel(time=slice(None, pd.Timestamp(end_date)))
    variables = [name for name in ("thetao", "so") if name in ds]
    if "thetao" not in variables or "so" not in variables:
        raise ValueError(f"{source_file} must contain thetao and so to compute the Figure 9 section and MLD.")

    field_total: xr.Dataset | None = None
    field_count: xr.Dataset | None = None
    mld_total = None
    mld_count = None

    for start in range(0, ds.sizes["time"], block_size):
        stop = min(start + block_size, ds.sizes["time"])
        block = ds.isel(time=slice(start, stop))
        field_total, field_count = _accumulate_time_mean(field_total, field_count, block, variables)

        mld = mld_from_density_threshold(block, density_threshold=density_threshold, z_ref=z_ref)
        block_mld_total = mld.sum("time", skipna=True).load()
        block_mld_count = mld.notnull().sum("time").load()
        mld_total = block_mld_total if mld_total is None else mld_total + block_mld_total
        mld_count = block_mld_count if mld_count is None else mld_count + block_mld_count
        print(f"Processed GLORYS 1000 m section block {start + 1}-{stop} of {ds.sizes['time']}")

    if field_total is None or field_count is None or mld_total is None or mld_count is None:
        raise ValueError(f"{source_file} has no time samples.")

    ds_mean = field_total / field_count.where(field_count > 0)
    ds_mean["mld"] = (mld_total / mld_count.where(mld_count > 0)).rename("mld")
    for var_name in ds_mean.data_vars:
        ds_mean[var_name] = ds_mean[var_name].astype("float32")
    ds_mean["mld"].attrs.update(
        {
            "standard_name": "ocean_mixed_layer_thickness",
            "long_name": "time-mean mixed layer depth from density threshold",
            "units": ds.depth.attrs.get("units", "m"),
            "density_threshold_kg_m3": density_threshold,
            "reference_depth_m": z_ref,
        }
    )
    ds_mean.attrs.update(ds.attrs)
    ds_mean.attrs.update(
        {
            "source_file": SOURCE_NAME,
            "time_mean_source_start": str(ds.time.min().values),
            "time_mean_source_end": str(ds.time.max().values),
            "time_mean_source_count": int(ds.sizes["time"]),
            "mld_density_threshold_kg_m3": density_threshold,
            "mld_reference_depth_m": z_ref,
        }
    )
    encoding = {
        var_name: {"zlib": True, "complevel": 4}
        for var_name in ds_mean.data_vars
    }
    ds_mean.to_netcdf(out_file, encoding=encoding)
    print(f"Wrote {out_file}")
    return out_file


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT)
    parser.add_argument("--force", action="store_true", help="Recompute even if the compact file exists.")
    parser.add_argument("--block-size", type=int, default=12, help="Number of source months to process at once.")
    parser.add_argument("--density-threshold", type=float, default=0.03)
    parser.add_argument("--z-ref", type=float, default=10.0)
    parser.add_argument("--end-date", default="2023-12-31", help="Last source date included in the time mean.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    build_glorys_section_timemean(
        project_root=args.project_root,
        force=args.force,
        block_size=args.block_size,
        density_threshold=args.density_threshold,
        z_ref=args.z_ref,
        end_date=args.end_date,
    )


if __name__ == "__main__":
    main()
