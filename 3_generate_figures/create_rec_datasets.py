"""Read R fPCA outputs and reconstruct gridded xarray datasets."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import xarray as xr


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _fpca_dir(project_root: Path, dataset_name: str, n_modes: str, dense: bool) -> Path:
    suffix = "_dense" if dense and not dataset_name.endswith("_dense") else ""
    if n_modes == "original":
        return project_root / "processed" / "2_fPCA" / f"{dataset_name}{suffix}"
    return project_root / f"processed_{n_modes}_modes" / "2_fPCA" / f"{dataset_name}{suffix}"


def r_analysis_df(
    dataset_name: str,
    n_modes: str = "original",
    dense: bool = False,
    project_root: str | Path | None = None,
):
    """Load fPCA CSV outputs and rebuild the reconstructed MLD dataset.

    Parameters
    ----------
    dataset_name:
        Base R-input/fPCA dataset name, e.g. ``"CMA_masked"``,
        ``"GLORYS_masked"``, or ``"GLORYS_CL_masked"``.
    n_modes:
        Name of an alternate fPCA output set. ``"original"`` maps to
        ``processed/2_fPCA``.
    dense:
        Append ``"_dense"`` to the fPCA output directory. This is used for the
        GLORYS dense run.
    project_root:
        Repository root. Defaults to the parent of ``3_generate_figures``.
    """
    root = Path(project_root).resolve() if project_root is not None else PROJECT_ROOT
    base_raw = root / "processed" / "1_gridded_data" / "r_input"
    base_processed = _fpca_dir(root, dataset_name, n_modes, dense)

    raw_file = base_raw / f"{dataset_name}.txt"
    if not raw_file.exists():
        raise FileNotFoundError(f"Missing R input table: {raw_file}")
    if not base_processed.exists():
        raise FileNotFoundError(f"Missing fPCA output directory: {base_processed}")

    df = (
        pd.read_csv(raw_file, sep=r"\s+")
        .assign(time=lambda d: pd.to_datetime(d["year"] * 100 + d["mth"], format="%Y%m"))
        .drop(columns=["year", "mth"])
    )

    ds_coord = (
        df.loc[df["time"].eq(df["time"].min()), ["long", "lat"]]
        .reset_index(drop=True)
        .to_xarray()
    )

    def _read_processed(filename: str, rename: dict[str, str] | None = None) -> pd.DataFrame:
        out = pd.read_csv(base_processed / filename).drop(columns="Unnamed: 0", errors="ignore")
        if rename is not None:
            out = out.rename(columns=rename)
        return out

    df_XIEST = _read_processed("PCA_XIEST.csv")
    df_PHI = _read_processed("PCA_PHI.csv")
    df_GRID = _read_processed("PCA_GRID.csv", rename={"x": "time"})
    df_MU = _read_processed("PCA_MU.csv").rename(columns={"x": "MU"})
    df_lambda = _read_processed("PCA_LAMBDA.csv")

    coords = ds_coord[["long", "lat"]].to_dataframe().reset_index()
    df_XIEST = coords.join(df_XIEST.reset_index(drop=True))

    time_index = pd.date_range(start="2007-01-01", periods=len(df_PHI), freq="ME")
    df_PHI = df_PHI.set_axis(time_index, axis=0).rename_axis("time")
    df_MU = df_MU.set_axis(time_index, axis=0).rename_axis("time")

    n_modes_loaded = len(df_PHI.columns)
    mode_numbers = range(1, n_modes_loaded + 1)

    ds_MU = df_MU.to_xarray()
    ds_PHI = df_PHI.to_xarray().rename({f"V{i}": f"phi{i}" for i in mode_numbers})
    ds_XIEST = (
        df_XIEST
        .set_index(["long", "lat"])
        .to_xarray()
        .rename({f"V{i}": f"xi{i}" for i in mode_numbers})
    )

    ds_PHI["phi1"] = -ds_PHI["phi1"]
    ds_XIEST["xi1"] = -ds_XIEST["xi1"]
    if dataset_name == "GLORYS_CL_masked" and "phi2" in ds_PHI and "xi2" in ds_XIEST:
        ds_PHI["phi2"] = -ds_PHI["phi2"]
        ds_XIEST["xi2"] = -ds_XIEST["xi2"]

    ds = xr.merge((ds_PHI, ds_MU, ds_XIEST))
    xi_phi_terms = []
    for i in mode_numbers:
        term_name = f"xi_phi{i}"
        ds[term_name] = ds[f"xi{i}"] * ds[f"phi{i}"]
        xi_phi_terms.append(ds[term_name])

    ds["xi_phi_tot"] = sum(xi_phi_terms[1:], xi_phi_terms[0])
    ds["mld"] = ds["MU"] + ds["xi_phi_tot"]

    return ds, df, df_XIEST, df_PHI, df_GRID, df_MU, df_lambda
