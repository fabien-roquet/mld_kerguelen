#!/usr/bin/env python3
"""Figure 1 conversion note.

Figure_1_prof_dist_year.ipynb uses individual raw CORA, MEOP, and ARGO profile
files that are not part of the reorganized two-input pipeline. The combined
``data/CORA_MEOP_ARGO_2026.nc`` file is enough for CMA processing, but it does
not preserve the source labels needed to reproduce the stacked CORA/MEOP profile
histograms exactly.
"""

from __future__ import annotations

import argparse

from figure_common import parse_project_root_arg, paths, require_file


def main() -> None:
    parser = parse_project_root_arg(argparse.ArgumentParser(description=__doc__))
    args = parser.parse_args()
    data_dir = paths(args.project_root)["data"]
    required = [
        data_dir / "new_all_CORA_2024.nc",
        data_dir / "all_MEOP_2024.nc",
        data_dir / "all_ARGO_2024.nc",
    ]
    for file in required:
        require_file(file, "Needed to reproduce Figure 1 exactly from the notebook.")
    raise SystemExit("Figure 1 script has the input checks in place, but plotting is not wired for the reduced two-input repo.")


if __name__ == "__main__":
    main()
