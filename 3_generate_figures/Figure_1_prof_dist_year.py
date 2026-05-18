#!/usr/bin/env python3
"""Add the KERFIX station marker to Figure 1.

Figure_1_prof_dist_year.ipynb uses individual raw CORA, MEOP, and ARGO profile
files that are not part of the reorganized two-input pipeline. The combined
``data/CORA_MEOP_ARGO_2026.nc`` file is enough for CMA processing, but it does
not preserve the source labels needed to reproduce the stacked CORA/MEOP profile
histograms exactly.

This script therefore applies the reproducible KERFIX overlay to the rendered
Figure 1 image that is distributed with the processed outputs.
"""

from __future__ import annotations

import argparse

import matplotlib.image as mpimg
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
import scipy.ndimage as ndi

from figure_common import parse_project_root_arg, paths


KERFIX_LON = 68.4167
KERFIX_LAT = -50.6667

# Pixel bounds of panel (a)'s map area in the rendered Figure 1 PNG. They map to
# the plotted domain 60-80 deg E and 60-40 deg S.
MAP_LEFT_PX = 332
MAP_RIGHT_PX = 3568
MAP_TOP_PX = 61
MAP_BOTTOM_PX = 2371
MAP_LON_MIN = 60.0
MAP_LON_MAX = 80.0
MAP_LAT_MIN = -60.0
MAP_LAT_MAX = -40.0
OLD_SECTION_A_BBOX = (1700, 1440, 1820, 1555)
NEW_SECTION_A_X = 2325
NEW_SECTION_A_Y = 1295


def data_to_pixel(lon: float, lat: float) -> tuple[float, float]:
    x = MAP_LEFT_PX + (lon - MAP_LON_MIN) / (MAP_LON_MAX - MAP_LON_MIN) * (MAP_RIGHT_PX - MAP_LEFT_PX)
    y = MAP_BOTTOM_PX - (lat - MAP_LAT_MIN) / (MAP_LAT_MAX - MAP_LAT_MIN) * (MAP_BOTTOM_PX - MAP_TOP_PX)
    return x, y


def remove_old_section_label(image):
    x0, y0, x1, y1 = OLD_SECTION_A_BBOX
    patch = image[y0:y1, x0:x1].copy()
    rgb = patch[..., :3]
    yellow_mask = (rgb[..., 0] > 0.75) & (rgb[..., 1] > 0.45) & (rgb[..., 1] < 0.9) & (rgb[..., 2] < 0.35)
    if yellow_mask.any():
        nearest = ndi.distance_transform_edt(yellow_mask, return_distances=False, return_indices=True)
        patch[yellow_mask] = patch[tuple(nearest[:, yellow_mask])]
        image[y0:y1, x0:x1] = patch
    return image


def main() -> None:
    parser = parse_project_root_arg(argparse.ArgumentParser(description=__doc__))
    args = parser.parse_args()

    out_file = paths(args.project_root)["figures"] / "Figure_1_prof_dist.png"
    reference_file = paths(args.project_root)["root"] / "processed_reference" / "3_figures" / "Figure_1_prof_dist.png"
    base_file = reference_file if reference_file.exists() else out_file
    if not base_file.exists():
        raise FileNotFoundError(
            f"Missing rendered Figure 1 image: {base_file}. The reduced repo cannot rebuild Figure 1 from raw source files."
        )

    image = remove_old_section_label(mpimg.imread(base_file))
    height, width = image.shape[:2]
    fig = plt.figure(figsize=(width / 100, height / 100), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(image)
    ax.axis("off")
    ax.set_xlim(0, width)
    ax.set_ylim(height, 0)

    kerfix_x, kerfix_y = data_to_pixel(KERFIX_LON, KERFIX_LAT)
    ax.scatter(
        kerfix_x,
        kerfix_y,
        marker="*",
        s=25000,
        facecolor="#FFD21F",
        edgecolor="black",
        linewidth=3.5,
        zorder=10,
    )
    ax.text(
        kerfix_x - 150,
        kerfix_y,
        "KERFIX",
        ha="right",
        va="center",
        fontsize=68,
        fontweight="bold",
        color="black",
        bbox={"facecolor": "white", "edgecolor": "black", "alpha": 0.9, "pad": 6},
        zorder=11,
    )
    ax.text(
        NEW_SECTION_A_X,
        NEW_SECTION_A_Y,
        "A",
        ha="center",
        va="center",
        fontsize=78,
        fontweight="bold",
        color="#FFD21F",
        path_effects=[path_effects.withStroke(linewidth=5, foreground="black")],
        zorder=11,
    )
    out_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_file, dpi=100, pad_inches=0)
    plt.close(fig)
    print(f"Wrote {out_file}")


if __name__ == "__main__":
    main()
