"""Pipeline crop sizes in µm (Visium pitch 100 µm).

spot / spot_margin match the original Visium FOVs.
non_overlap is the largest axis-aligned square with no interior overlap on a
flat-top hex lattice (√3/2 × pitch).
hest matches HEST-1k / H-optimus-1 / Virchow2@20× / UNI2 resized FOV
(224 px × 0.5 µm/px).
"""

from __future__ import annotations

import math

PITCH_UM = 100.0
TARGET_MPP = 0.5  # harmonized read resolution for lazyslide TileSpec

CROP_UM: dict[str, float] = {
    "spot": 55.0,
    "spot_margin": 75.0,
    "non_overlap": PITCH_UM * math.sqrt(3) / 2,
    "hest": 224 * TARGET_MPP,
}

CROP_ORDER: list[str] = ["spot", "spot_margin", "non_overlap", "hest"]


def tile_px(crop_name: str, mpp: float = TARGET_MPP) -> int:
    if crop_name not in CROP_UM:
        raise KeyError(f"unknown crop_size {crop_name!r}; choose from {CROP_ORDER}")
    return max(1, int(round(CROP_UM[crop_name] / mpp)))
