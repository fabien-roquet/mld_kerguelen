#!/usr/bin/env python3
"""Compare generated processed outputs against processed_reference."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _relative_files(root: Path, suffixes: tuple[str, ...]) -> set[Path]:
    return {
        path.relative_to(root)
        for path in root.rglob("*")
        if path.is_file() and path.suffix.lower() in suffixes and path.name != ".DS_Store"
    }


def _max_abs_diff(a: np.ndarray, b: np.ndarray) -> float:
    if a.size == 0 and b.size == 0:
        return 0.0
    diff = np.abs(a - b)
    if np.all(np.isnan(diff)):
        return 0.0
    return float(np.nanmax(diff))


def compare_netcdf(generated: Path, reference: Path, rtol: float, atol: float) -> list[str]:
    issues = []
    with xr.open_dataset(generated) as ds_gen, xr.open_dataset(reference) as ds_ref:
        if dict(ds_gen.sizes) != dict(ds_ref.sizes):
            issues.append(f"sizes differ: generated={dict(ds_gen.sizes)} reference={dict(ds_ref.sizes)}")

        all_coords = sorted(set(ds_gen.coords) | set(ds_ref.coords))
        for name in all_coords:
            if name not in ds_gen.coords or name not in ds_ref.coords:
                issues.append(f"coordinate {name!r} missing in one dataset")
                continue
            left, right = xr.align(ds_gen[name], ds_ref[name], join="outer")
            if left.shape != right.shape:
                issues.append(f"coordinate {name!r} shape differs: {left.shape} vs {right.shape}")
                continue
            if np.issubdtype(left.dtype, np.number):
                max_diff = _max_abs_diff(left.values, right.values)
                if not np.allclose(left.values, right.values, rtol=rtol, atol=atol, equal_nan=True):
                    issues.append(f"coordinate {name!r} differs; max abs diff={max_diff:.6g}")
            elif not np.array_equal(left.values, right.values):
                issues.append(f"coordinate {name!r} differs")

        all_vars = sorted(set(ds_gen.data_vars) | set(ds_ref.data_vars))
        for name in all_vars:
            if name not in ds_gen.data_vars or name not in ds_ref.data_vars:
                issues.append(f"variable {name!r} missing in one dataset")
                continue
            left, right = xr.align(ds_gen[name], ds_ref[name], join="outer")
            if left.dims != right.dims:
                issues.append(f"variable {name!r} dims differ: {left.dims} vs {right.dims}")
                continue
            if left.shape != right.shape:
                issues.append(f"variable {name!r} shape differs: {left.shape} vs {right.shape}")
                continue
            if np.issubdtype(left.dtype, np.number):
                max_diff = _max_abs_diff(left.values, right.values)
                if not np.allclose(left.values, right.values, rtol=rtol, atol=atol, equal_nan=True):
                    issues.append(f"variable {name!r} differs; max abs diff={max_diff:.6g}")
            elif not np.array_equal(left.values, right.values):
                issues.append(f"variable {name!r} differs")
    return issues


def compare_csv(generated: Path, reference: Path, rtol: float, atol: float) -> list[str]:
    issues = []
    gen = pd.read_csv(generated)
    ref = pd.read_csv(reference)
    gen = gen.drop(columns="Unnamed: 0", errors="ignore")
    ref = ref.drop(columns="Unnamed: 0", errors="ignore")

    if gen.shape != ref.shape:
        issues.append(f"shape differs: generated={gen.shape} reference={ref.shape}")
        return issues
    if list(gen.columns) != list(ref.columns):
        issues.append(f"columns differ: generated={list(gen.columns)} reference={list(ref.columns)}")
        return issues

    for column in gen.columns:
        if pd.api.types.is_numeric_dtype(gen[column]) and pd.api.types.is_numeric_dtype(ref[column]):
            left = gen[column].to_numpy()
            right = ref[column].to_numpy()
            max_diff = _max_abs_diff(left, right)
            if not np.allclose(left, right, rtol=rtol, atol=atol, equal_nan=True):
                issues.append(f"column {column!r} differs; max abs diff={max_diff:.6g}")
        elif not gen[column].equals(ref[column]):
            issues.append(f"column {column!r} differs")
    return issues


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT)
    parser.add_argument("--generated", type=Path, default=None, help="Generated processed directory.")
    parser.add_argument("--reference", type=Path, default=None, help="Reference processed directory.")
    parser.add_argument("--rtol", type=float, default=1e-5)
    parser.add_argument("--atol", type=float, default=1e-4)
    parser.add_argument("--fail-missing", action="store_true", help="Fail when one side has extra files.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    project_root = args.project_root.resolve()
    generated = args.generated.resolve() if args.generated else project_root / "processed"
    reference = args.reference.resolve() if args.reference else project_root / "processed_reference"

    if not generated.exists():
        raise FileNotFoundError(f"Generated processed directory not found: {generated}")
    if not reference.exists():
        raise FileNotFoundError(f"Reference processed directory not found: {reference}")

    generated_files = _relative_files(generated, (".nc", ".csv"))
    reference_files = _relative_files(reference, (".nc", ".csv"))
    common_files = sorted(generated_files & reference_files)
    only_generated = sorted(generated_files - reference_files)
    only_reference = sorted(reference_files - generated_files)

    failed = False
    for rel_path in common_files:
        gen_file = generated / rel_path
        ref_file = reference / rel_path
        if rel_path.suffix.lower() == ".nc":
            issues = compare_netcdf(gen_file, ref_file, rtol=args.rtol, atol=args.atol)
        else:
            issues = compare_csv(gen_file, ref_file, rtol=args.rtol, atol=args.atol)

        if issues:
            failed = True
            print(f"[FAIL] {rel_path}")
            for issue in issues:
                print(f"  - {issue}")
        else:
            print(f"[OK] {rel_path}")

    if only_generated:
        print("Only in generated:")
        for path in only_generated:
            print(f"  - {path}")
        failed = failed or args.fail_missing

    if only_reference:
        print("Only in reference:")
        for path in only_reference:
            print(f"  - {path}")
        failed = failed or args.fail_missing

    if failed:
        raise SystemExit(1)

    print(f"Compared {len(common_files)} files successfully.")


if __name__ == "__main__":
    main()
