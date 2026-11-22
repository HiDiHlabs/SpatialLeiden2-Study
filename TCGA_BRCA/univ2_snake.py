import argparse
import os
import sys
from pathlib import Path

import numpy as np
import scanpy as sc
import torch
import timm
from PIL import Image
from tqdm import tqdm
from torchvision import transforms

sys.path.insert(0, str(Path(__file__).resolve().parent))
from thin_obsm import write_thin_obsm

parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True)
parser.add_argument('--output', required=True)
parser.add_argument('--device', default='auto', choices=('auto', 'cuda', 'cpu'),
                    help='Device: auto (cuda if available else cpu), cuda, or cpu')
parser.add_argument('--batch-size', type=int, default=32,
                    help='Tiles per batch (default: 32; use smaller on CPU, e.g. 8)')
parser.add_argument('--model-dir', default=None,
                    help='Directory containing pytorch_model.bin (UNI2 weights). Override with --model-dir or UNIV2_MODEL_DIR.')
args = parser.parse_args()

if args.device == 'auto':
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
else:
    device = torch.device(args.device)
print(f"Using device: {device}")

adata = sc.read_h5ad(args.input)
tile_paths = adata.obs['tile_path'].tolist()
batch_size = max(1, args.batch_size)

# UNI2 model setup: directory containing pytorch_model.bin (set via --model-dir or UNIV2_MODEL_DIR)
model_dir = args.model_dir or os.environ.get("UNIV2_MODEL_DIR", "models")
model_path = os.path.join(model_dir, "pytorch_model.bin")
if not os.path.isfile(model_path):
    raise FileNotFoundError(f"UNI2 weights not found at {model_path}. Set --model-dir or UNIV2_MODEL_DIR to the directory containing pytorch_model.bin.")
timm_kwargs = {
    'model_name': 'vit_giant_patch14_224',
    'img_size': 224,
    'patch_size': 14,
    'depth': 24,
    'num_heads': 24,
    'init_values': 1e-5,
    'embed_dim': 1536,
    'mlp_ratio': 2.66667*2,
    'num_classes': 0,
    'no_embed_class': True,
    'mlp_layer': timm.layers.SwiGLUPacked,
    'act_layer': torch.nn.SiLU,
    'reg_tokens': 8,
    'dynamic_img_size': True
}
model = timm.create_model(pretrained=False, **timm_kwargs)
model.load_state_dict(torch.load(model_path, map_location="cpu"), strict=True)
# Resize to (224, 224) so all tensors match when stacking (Resize(224) keeps aspect ratio -> variable size)
transform = transforms.Compose([
    transforms.Resize((224, 224)),
    transforms.ToTensor(),
    transforms.Normalize(mean=(0.485, 0.456, 0.406), std=(0.229, 0.224, 0.225)),
])
model.to(device)
model.eval()

# Batched feature extraction
all_features = []
for i in tqdm(range(0, len(tile_paths), batch_size), desc="univ2"):
    batch_paths = tile_paths[i : i + batch_size]
    imgs = []
    for path in batch_paths:
        img = Image.open(path)
        img = transform(img)
        imgs.append(img)
    batch = torch.stack(imgs).to(device)
    with torch.inference_mode():
        feat = model(batch)
    all_features.append(feat.cpu())
features = torch.cat(all_features, dim=0).numpy()
write_thin_obsm(adata, features, 'univ2', args.output) 