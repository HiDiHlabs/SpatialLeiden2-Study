"""Download Kuppe et al. 2022 (Nature 608, 766-777) data from Zenodo.

The study is split over several Zenodo records. Choose what to download with
--sets (default: the "relevant" sets below, about 82 GB):

    processed1   6578047  27 Visium h5ad + snRNA-seq + snATAC-seq h5ad    25.4 GB
    visium_h5ad  6578047  only the 27 Visium h5ad of processed1            4.4 GB
    processed2   7098004  per-cell-type snRNA+snATAC objects (R .Rds)     34.1 GB
    visium_raw   6580069  Visium Space Ranger outputs + metadata          17.5 GB
    he_images    6778620  H&E images of the Visium sections + metadata     5.1 GB
    metadata     6578617  RNA and ATAC sample metadata sheets             <1 MB
    atac_frag    6578553 + 6578617  snATAC fragment files                63 GB   (not default)
    rna_raw      6578617  snRNA Cell Ranger outputs                       21 GB   (not default)

Run from this folder:

    python kuppe_2022.py                         # relevant sets
    python kuppe_2022.py --sets processed1       # one set
    python kuppe_2022.py --sets atac_frag
    python kuppe_2022.py --sets visium_h5ad visium_raw metadata   # all Visium, ~22 GB

Files go to ../data/kuppe_2022/<set>/. Files already present with a matching
checksum are skipped. The script stops if less than 50 GB would remain free.
"""

import argparse
import hashlib
import json
import shutil
import time
import urllib.error
import urllib.request
from pathlib import Path

OUT_ROOT = Path(__file__).resolve().parent.parent / "data" / "kuppe_2022"
MIN_FREE_GB = 50

# set name -> list of (zenodo record, filter on file name)
SETS = {
    "processed1": [("6578047", lambda n: True)],
    "visium_h5ad": [("6578047", lambda n: n.startswith("Visium_"))],
    "processed2": [("7098004", lambda n: True)],
    "visium_raw": [("6580069", lambda n: True)],
    "he_images": [("6778620", lambda n: True)],
    "metadata": [("6578617", lambda n: n.startswith("metadata"))],
    "atac_frag": [
        ("6578553", lambda n: "_ATAC_" in n),
        ("6578617", lambda n: "_ATAC_" in n),
    ],
    "rna_raw": [("6578617", lambda n: "_RNA_" in n)],
}
DEFAULT = ["processed1", "processed2", "visium_raw", "he_images", "metadata"]


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 24), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(url: str, dest: Path, size: int, retries: int = 20) -> None:
    """Download url to dest, resuming after dropped connections."""
    for attempt in range(1, retries + 1):
        have = dest.stat().st_size if dest.exists() else 0
        if have >= size:
            return
        req = urllib.request.Request(url, headers={"Range": f"bytes={have}-"})
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                mode = "ab" if r.status == 206 else "wb"  # 200 = server ignored Range
                with open(dest, mode) as f:
                    shutil.copyfileobj(r, f, 1 << 24)
        except (OSError, urllib.error.URLError) as e:
            print(f"       retry {attempt}/{retries} after: {e}", flush=True)
            time.sleep(min(60, 5 * attempt))
    if not dest.exists() or dest.stat().st_size < size:
        raise SystemExit(f"could not download {url} after {retries} attempts")


def list_files(record: str) -> list[dict]:
    with urllib.request.urlopen(f"https://zenodo.org/api/records/{record}") as r:
        return json.load(r)["files"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--sets", nargs="+", choices=SETS, default=DEFAULT)
    args = parser.parse_args()

    jobs = []
    for name in args.sets:
        for record, keep in SETS[name]:
            jobs += [(name, f) for f in list_files(record) if keep(f["key"])]

    total = sum(f["size"] for _, f in jobs)
    print(f"{len(jobs)} files, {total / 1e9:.1f} GB in sets {args.sets} -> {OUT_ROOT}")

    for name, f in sorted(jobs, key=lambda j: j[1]["size"]):
        out_dir = OUT_ROOT / name
        out_dir.mkdir(parents=True, exist_ok=True)
        key, expected = f["key"], f["checksum"].removeprefix("md5:")
        path = out_dir / key
        if path.exists() and md5(path) == expected:
            print(f"skip   {name}/{key}")
            continue
        free_after = (shutil.disk_usage(OUT_ROOT).free - f["size"]) / 1e9
        if free_after < MIN_FREE_GB:
            raise SystemExit(f"stop: {name}/{key} would leave {free_after:.0f} GB free")
        print(f"get    {name}/{key} ({f['size'] / 1e9:.2f} GB)", flush=True)
        tmp = path.with_suffix(path.suffix + ".part")
        fetch(f["links"]["self"], tmp, f["size"])
        if md5(tmp) != expected:
            tmp.unlink()
            raise SystemExit(f"checksum mismatch for {name}/{key}")
        tmp.rename(path)

    print("done")


if __name__ == "__main__":
    main()
