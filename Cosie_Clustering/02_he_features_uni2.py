"""
UNI2-h features for every Visium spot of the COSIE mouse brain sections.

Only the Space Ranger hires H&E image is available (~3.5 um/px, one 55 um spot ~ 16 px),
no full-resolution image. UNI2 is trained on 224 px tiles at ~0.5 um/px (~112 um field
of view), so for every spot we crop the same field of view (~32 px) around the spot
centre and upsample it to 224 px. This runs, but the input is far blurrier than what
UNI2 was trained on; with full-resolution images the features would be much better.

Requirements: access to https://huggingface.co/MahmoodLab/UNI2-h, `huggingface-cli login`,
and an env with torch, timm and huggingface_hub. A GPU is recommended (CPU works, slower).

Run from the SpatialLeiden2-Study folder:
    python Cosie_Clustering/02_he_features_uni2.py
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import timm
import torch
from PIL import Image
from timm.data import resolve_data_config
from timm.data.transforms_factory import create_transform

data_dir = Path("data/01_mouse_brain_MALDI_visium")
slide_dir = data_dir / "sma" / "V11L12-109"
out_dir = data_dir / "processed" / "he_uni2"
out_dir.mkdir(parents=True, exist_ok=True)

SECTIONS = {
    "V11L12-109_A1": slide_dir / "V11L12-109_A1/output_data/V11L12-109_A1_RNA/outs/spatial",
    "V11L12-109_B1": slide_dir / "V11L12-109_B1/output_data/V11L12-109_B1_RNA/outs/spatial",
    "V11L12-109_C1": slide_dir / "V11L12-109_C1/output_data/V11L12-109_C1_RNA/outs/spatial",
    "V11L12-109_D1": slide_dir / "V11L12-109_D1/output_data/V11L12-109_D1/outs/spatial",
}

# Space Ranger's spot_diameter_fullres corresponds to 65 um (55 um spot plus padding)
spot_diameter_um = 65
field_of_view_um = 112  # what a 224 px UNI2 tile covers at 0.5 um/px
batch_size = 64
device = "cuda" if torch.cuda.is_available() else "cpu"

# model definition from the UNI2-h model card
timm_kwargs = {
    "img_size": 224, "patch_size": 14, "depth": 24, "num_heads": 24, "init_values": 1e-5,
    "embed_dim": 1536, "mlp_ratio": 2.66667 * 2, "num_classes": 0, "no_embed_class": True,
    "mlp_layer": timm.layers.SwiGLUPacked, "act_layer": torch.nn.SiLU, "reg_tokens": 8,
    "dynamic_img_size": True,
}
model = timm.create_model("hf-hub:MahmoodLab/UNI2-h", pretrained=True, **timm_kwargs).eval().to(device)
transform = create_transform(**resolve_data_config(model.pretrained_cfg, model=model))


def spot_patches(spatial_dir):
    img = Image.open(spatial_dir / "tissue_hires_image.png").convert("RGB")
    scalef = json.loads((spatial_dir / "scalefactors_json.json").read_text())
    um_per_px = spot_diameter_um / (scalef["spot_diameter_fullres"] * scalef["tissue_hires_scalef"])
    half = int(round(field_of_view_um / um_per_px / 2))

    pos = pd.read_csv(
        spatial_dir / "tissue_positions_list.csv", header=None,
        names=["barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"],
    )
    pos = pos[pos["in_tissue"] == 1]
    rows = (pos["pxl_row"] * scalef["tissue_hires_scalef"]).round().astype(int)
    cols = (pos["pxl_col"] * scalef["tissue_hires_scalef"]).round().astype(int)

    # pad with white so crops at the image border keep their size
    arr = np.pad(np.asarray(img), ((half, half), (half, half), (0, 0)), constant_values=255)
    patches = [
        Image.fromarray(arr[r : r + 2 * half, c : c + 2 * half]).resize((224, 224), Image.BICUBIC)
        for r, c in zip(rows, cols)
    ]
    print(f"  {um_per_px:.2f} um/px, crop {2 * half} px -> 224 px, {len(patches):,} spots")
    return pos["barcode"].to_numpy(), patches


@torch.inference_mode()
def embed(patches):
    feats = []
    for i in range(0, len(patches), batch_size):
        x = torch.stack([transform(p) for p in patches[i : i + batch_size]]).to(device)
        with torch.autocast(device_type=device, dtype=torch.float16, enabled=device == "cuda"):
            feats.append(model(x).float().cpu().numpy())
    return np.concatenate(feats)


for section, spatial_dir in SECTIONS.items():
    print(section)
    barcodes, patches = spot_patches(spatial_dir)
    feats = embed(patches)
    df = pd.DataFrame(feats, index=[f"{section}_{bc}" for bc in barcodes])
    df.columns = [f"uni2_{i}" for i in range(df.shape[1])]
    df.to_csv(out_dir / f"{section}.csv.gz")
    print(f"  -> {out_dir / f'{section}.csv.gz'} ({df.shape[0]:,} x {df.shape[1]})")
