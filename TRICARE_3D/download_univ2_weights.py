#!/usr/bin/env python3
"""Download gated UNI2-h weights (pytorch_model.bin) from Hugging Face.

Repo: https://huggingface.co/MahmoodLab/UNI2-h  (request access first;
institutional email required). Then huggingface-cli login.

Usage
-----
  conda activate tricare_univ2   # or any env with huggingface_hub
  # login (cli entrypoint varies by huggingface_hub version):
  python -m huggingface_hub.commands.huggingface_cli login
  # or:  hf auth login

  export MODEL_DIR=./models
  python TRICARE_3D/download_univ2_weights.py --model-dir "$MODEL_DIR"

  # then TRICARE univ2 array:
  export MODEL_DIR=./models
  bash TRICARE_3D/submit_tricare_array.sh univ2
"""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        "--model-dir",
        type=Path,
        default=Path("models/uni2-h"),
        help="Directory to place pytorch_model.bin (default: models/uni2-h)",
    )
    p.add_argument("--repo-id", default="MahmoodLab/UNI2-h")
    p.add_argument("--filename", default="pytorch_model.bin")
    args = p.parse_args()

    try:
        from huggingface_hub import hf_hub_download
        from huggingface_hub.utils import HfHubHTTPError
    except ImportError as e:
        raise SystemExit("Install huggingface_hub: pip install huggingface_hub") from e

    out = args.model_dir.expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)
    dest = out / args.filename
    if dest.is_file() and dest.stat().st_size > 0:
        print(f"Already present: {dest}")
        print(f"export MODEL_DIR={out}")
        return

    print(f"Downloading {args.repo_id}/{args.filename} → {out}")
    try:
        path = hf_hub_download(
            repo_id=args.repo_id,
            filename=args.filename,
            local_dir=str(out),
        )
    except HfHubHTTPError as e:
        raise SystemExit(
            f"{e}\n"
            "Gated model: request access at https://huggingface.co/MahmoodLab/UNI2-h\n"
            "Then login:\n"
            "  python -m huggingface_hub.commands.huggingface_cli login\n"
            "  # or: hf auth login"
        ) from e
    print(f"Wrote {path}")
    print(f"export MODEL_DIR={out}")


if __name__ == "__main__":
    main()
