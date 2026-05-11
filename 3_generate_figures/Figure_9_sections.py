#!/usr/bin/env python3
"""Figure 9 conversion note.

Figure_9_sections.ipynb depends on section-specific GLORYS/CMA files that are
not included in the reorganized two-input pipeline. This script records the
required inputs and fails early with a clear message if they are absent.
"""

from __future__ import annotations

import argparse

from figure_common import parse_project_root_arg, paths, require_file


def main() -> None:
    parser = parse_project_root_arg(argparse.ArgumentParser(description=__doc__))
    args = parser.parse_args()
    data_dir = paths(args.project_root)["data"]
    required = [
        data_dir / "GLORYS_1000m_section_plot.nc",
        data_dir / "GLORYS_2026_with_mld.nc",
        data_dir / "all_section_CMA.nc",
    ]
    for file in required:
        require_file(file, "Needed to reproduce Figure 9 exactly from the notebook.")
    raise SystemExit("Figure 9 script has the input checks in place, but plotting is not wired for missing section inputs.")


if __name__ == "__main__":
    main()
