#!/usr/bin/env python3
"""Run the modular three-stage MLD Kerguelen analysis."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parent

DATA_SCRIPTS = {
    "CMA": PROJECT_ROOT / "1_data_processing" / "process_CMA.py",
    "GLORYS": PROJECT_ROOT / "1_data_processing" / "process_GLORYS.py",
    "GLORYS_CL": PROJECT_ROOT / "1_data_processing" / "process_GLORYS_CL.py",
    "GLORYS_SECTION": PROJECT_ROOT / "1_data_processing" / "process_GLORYS_section.py",
}

FPCA_SCRIPTS = {
    "CMA": PROJECT_ROOT / "2_compute_fPCA_R" / "script_PCA_CM_2026.R",
    "GLORYS": PROJECT_ROOT / "2_compute_fPCA_R" / "script_PCA_GLORYS_2026.R",
    "GLORYS_CL": PROJECT_ROOT / "2_compute_fPCA_R" / "script_PCA_GLORYS_CL_2026.R",
}

FIGURE_SCRIPTS = {
    "1": PROJECT_ROOT / "3_generate_figures" / "Figure_1_prof_dist_year.py",
    "2": PROJECT_ROOT / "3_generate_figures" / "Figure_2_map_mld.py",
    "3": PROJECT_ROOT / "3_generate_figures" / "Figure_3_seasonal_cycle.py",
    "4": PROJECT_ROOT / "3_generate_figures" / "Figure_4_MU_modes.py",
    "5": PROJECT_ROOT / "3_generate_figures" / "Figure_5_maps_cp.py",
    "6": PROJECT_ROOT / "3_generate_figures" / "Figure_6_quadratic_error.py",
    "7": PROJECT_ROOT / "3_generate_figures" / "Figure_7_1D_trends.py",
    "8": PROJECT_ROOT / "3_generate_figures" / "Figure_8_trend_maps.py",
    "9": PROJECT_ROOT / "3_generate_figures" / "Figure_9_sections.py",
    "10": PROJECT_ROOT / "3_generate_figures" / "Figure_10_KERFIX.py",
}

DEFAULT_FIGURES = ["2", "3", "4", "5", "6", "7", "8", "9", "10"]


def run_command(command: list[str], cwd: Path = PROJECT_ROOT) -> None:
    print("+ " + " ".join(command), flush=True)
    subprocess.run(command, cwd=cwd, check=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stage",
        nargs="+",
        choices=["data", "fpca", "figures", "compare"],
        default=["data", "fpca", "figures", "compare"],
        help="Pipeline stages to run.",
    )
    parser.add_argument("--datasets", nargs="+", choices=DATA_SCRIPTS.keys(), default=list(DATA_SCRIPTS.keys()))
    parser.add_argument("--figures", nargs="+", choices=FIGURE_SCRIPTS.keys(), default=DEFAULT_FIGURES)
    parser.add_argument("--force-data", action="store_true", help="Recompute data-processing outputs.")
    parser.add_argument("--skip-reference-compare", action="store_true", help="Alias for omitting the compare stage.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    stages = set(args.stage)
    if args.skip_reference_compare:
        stages.discard("compare")

    if "data" in stages:
        for dataset in args.datasets:
            command = [sys.executable, str(DATA_SCRIPTS[dataset]), "--project-root", str(PROJECT_ROOT)]
            if args.force_data:
                command.append("--force")
            run_command(command)

    if "fpca" in stages:
        for dataset in args.datasets:
            script = FPCA_SCRIPTS.get(dataset)
            if script is not None:
                run_command(["Rscript", str(script)])

    if "figures" in stages:
        for figure in args.figures:
            run_command([sys.executable, str(FIGURE_SCRIPTS[figure]), "--project-root", str(PROJECT_ROOT)])

    if "compare" in stages:
        run_command([sys.executable, str(PROJECT_ROOT / "scripts" / "compare_processed_reference.py")])


if __name__ == "__main__":
    main()
