import xarray as xr
import pandas as pd

def r_analysis_df(dataset_name, n_modes="original",dense=False):
    ### Retrieve lat/lon from original df
    base_raw = "/home/jupyter-vincent2/vincent/process_profiles/data/R_analysis_2026/raw"

    if n_modes == "original":
        base_processed = f"/home/jupyter-vincent2/vincent/process_profiles/data/R_analysis_2026/processed/{dataset_name}"
    else:
        base_processed = f"/home/jupyter-vincent2/vincent/process_profiles/data/R_analysis_2026/processed_20_modes/{dataset_name}"

    if dense:
        base_processed += "_dense"
    
    df = (
        pd.read_csv(f"{base_raw}/{dataset_name}.txt", sep=" ")
        .assign(time=lambda d: pd.to_datetime(d["year"] * 100 + d["mth"], format="%Y%m"))
        .drop(columns=["year", "mth"])
    )

    ds_coord = (
        df.loc[df["time"].eq(df["time"].min()), ["long", "lat"]]
        .reset_index(drop=True)
        .to_xarray()
    )

    def _read_processed(filename, rename=None):
        out = pd.read_csv(f"{base_processed}/{filename}").drop(columns="Unnamed: 0", errors="ignore")
        if rename is not None:
            out = out.rename(columns=rename)
        return out

    ### Parameters from PACE method
    ### Xiest
    df_XIEST = _read_processed("PCA_XIEST.csv")
    ### PHI
    df_PHI = _read_processed("PCA_PHI.csv")
    ### Grid
    df_GRID = _read_processed("PCA_GRID.csv", rename={"x": "time"})
    print(len(df_GRID))
    ### MU 
    df_MU = _read_processed("PCA_MU.csv").rename(columns={"x": "MU"})
    ### Variance explained by each mode
    df_lambda = _read_processed("PCA_LAMBDA.csv")

    ### Replace row index with lon/lat coordinates
    coords = ds_coord[["long", "lat"]].to_dataframe().reset_index()
    df_XIEST = coords.join(df_XIEST.reset_index(drop=True))
    # df_XIEST["V1"] = -df_XIEST["V1"]

    ### Time index
    time_index = pd.date_range(start="2007-01", end="2024-01", freq="ME")
    df_PHI = df_PHI.set_axis(time_index, axis=0).rename_axis("time")
    # df_PHI["V1"] = -df_PHI["V1"]
    df_MU = df_MU.set_axis(time_index, axis=0).rename_axis("time")

    ### To xarray
    n_modes = len(df_PHI.columns)
    mode_numbers = range(1, n_modes + 1)

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

    ds = xr.merge((ds_PHI, ds_MU, ds_XIEST))

    ### Check on the number of modes
    xi_phi_terms = []
    for i in mode_numbers:
        term_name = f"xi_phi{i}"
        ds[term_name] = ds[f"xi{i}"] * ds[f"phi{i}"]
        xi_phi_terms.append(ds[term_name])

    ds["xi_phi_tot"] = sum(xi_phi_terms[1:], xi_phi_terms[0])
    ds["mld"] = ds["MU"] + ds["xi_phi_tot"]

    return ds,df,df_XIEST,df_PHI,df_GRID,df_MU,df_lambda