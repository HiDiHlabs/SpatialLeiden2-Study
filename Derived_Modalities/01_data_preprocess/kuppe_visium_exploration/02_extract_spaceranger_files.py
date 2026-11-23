"""Copy small Space Ranger files out of each Kuppe 2022 Visium archive.

For every 10X_Visium_<id>.tar.gz writes, into the tables folder:
<id>_lowres.png, <id>_scalefactors.json, <id>_metrics_summary.csv.
Usage: python 02_extract_spaceranger_files.py <visium_raw dir> <tables dir>
"""

import io
import sys
import tarfile
import zipfile
from pathlib import Path

raw_dir, out_dir = Path(sys.argv[1]), Path(sys.argv[2])
out_dir.mkdir(parents=True, exist_ok=True)
for tar_path in sorted(raw_dir.glob("10X_Visium_*.tar.gz")):
    sid = tar_path.name.removeprefix("10X_Visium_").removesuffix(".tar.gz")
    if (out_dir / f"{sid}_lowres.png").exists():
        continue
    with tarfile.open(tar_path) as t:
        outs = f"./{sid}/outs/"
        metrics = t.extractfile(f"{outs}spatial_{sid}_metrics_summary.csv").read()
        (out_dir / f"{sid}_metrics_summary.csv").write_bytes(metrics)
        z = zipfile.ZipFile(
            io.BytesIO(t.extractfile(f"{outs}spatial_{sid}_spatial.zip").read())
        )
        for name in z.namelist():
            if name.endswith("tissue_lowres_image.png"):
                (out_dir / f"{sid}_lowres.png").write_bytes(z.read(name))
            elif name.endswith("scalefactors_json.json"):
                (out_dir / f"{sid}_scalefactors.json").write_bytes(z.read(name))
    print("done", sid)
