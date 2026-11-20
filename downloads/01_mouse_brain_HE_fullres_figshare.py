#!/usr/bin/env python3
"""
COSIE-Datensatz 1, Ergänzung: Vollauflösungs-H&E-Bilder (rohe Mikroskopbilder) zu
Vicari et al., Nat Biotechnol 2024. SciLifeLab Figshare, DOI 10.17044/scilifelab.22778920
  https://figshare.scilifelab.se/articles/dataset/Spatial_Multimodal_Analysis_SMA_-_Spatial_Transcriptomics/22778920

Der Mendeley-Datensatz (01_mouse_brain_MALDI_visium_mendeley.py) enthält nur die
verkleinerten Space-Ranger-Bilder (~3.5 um/px). Hier liegen die Originalbilder, auf die
sich pxl_row/pxl_col in tissue_positions_list.csv beziehen (~0.15 um/px).

Standard: die vier Striatum-Schnitte V11L12-109 A1-D1 (Bilder + Loupe-Alignment-JSON).
Mit --list nur anzeigen. Mit --match REGEX andere Dateien wählen (z.B. "Spot.*\\.jpg$"
für alle Mikroskopbilder). Die Dateien werden per MD5 gegen die Figshare-Angabe geprüft.
"""
import argparse
import hashlib
import json
import os
import re
import sys
import urllib.request

ARTICLE = 22778920
API = f"https://api.figshare.com/v2/articles/{ARTICLE}"
DEFAULT_MATCH = r"V11L12-109[._-][A-D]1-Spot\d+\.jpg$|^V11L12-109-[A-D]1\.json$"


def get(url):
    req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "figshare_dl/1.0"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.load(r)


def list_all():
    return [
        {"name": f["name"], "size": f["size"], "url": f["download_url"], "md5": f.get("computed_md5")}
        for f in get(API)["files"]
    ]


def md5sum(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        while chunk := fh.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


def download(f, dest):
    if os.path.exists(dest) and os.path.getsize(dest) == f["size"]:
        print(f"  [da]  {dest}")
        return
    os.makedirs(os.path.dirname(dest) or ".", exist_ok=True)
    tmp = dest + ".part"
    req = urllib.request.Request(f["url"], headers={"User-Agent": "figshare_dl/1.0"})
    with urllib.request.urlopen(req, timeout=600) as r, open(tmp, "wb") as fh:
        while chunk := r.read(1 << 20):
            fh.write(chunk)
    if f["md5"] and md5sum(tmp) != f["md5"]:
        raise RuntimeError("MD5 stimmt nicht, Datei bleibt als .part liegen")
    os.replace(tmp, dest)
    print(f"  [ok]  {dest}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--out", default="data/01_mouse_brain_MALDI_visium/he_fullres")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--match", default=DEFAULT_MATCH, help="Regex auf Dateiname")
    a = ap.parse_args()

    rx = re.compile(a.match, re.I)
    files = [f for f in list_all() if rx.search(f["name"])]
    total = sum(f["size"] for f in files)
    print(f"{len(files)} Dateien, {total/1e9:.2f} GB")
    for f in files:
        print(f"  {f['size']/1e6:9.1f} MB  {f['name']}")
    if a.list:
        return
    for f in files:
        try:
            download(f, os.path.join(a.out, f["name"]))
        except Exception as e:
            print(f"  [FEHLER] {f['name']}: {e}", file=sys.stderr)


if __name__ == "__main__":
    main()
