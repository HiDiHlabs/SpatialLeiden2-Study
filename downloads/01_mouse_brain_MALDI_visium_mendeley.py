#!/usr/bin/env python3
"""
COSIE-Datensatz 1: Mouse Brain, Spatial MALDI-Visium (SMA), Parkinson-Modell (6-OHDA)
H&E + Metabolom (MALDI-MSI) + RNA (Visium). Vicari et al., Nat Biotechnol 2024.
Mendeley Data, DOI 10.17632/w7nw4km7xd.1

Lädt alle Dateien (inkl. Unterordner) über die öffentliche Mendeley-Data-API.
Mit --list nur anzeigen. Mit --match REGEX nur passende Pfade laden.
Falls die API nicht antwortet: Seite im Browser öffnen und "Download All" klicken.
  https://data.mendeley.com/datasets/w7nw4km7xd/1

Hinweis: Der Datensatz enthält mehr als COSIE nutzt (auch Humanproben, Kontrollen).
COSIE: 3 Sections (Section 1 von Maus A, Sections 2+3 von Maus B; Section 3 ohne MALDI).
"""
import argparse
import json
import os
import re
import sys
import urllib.request

DS, VER = "w7nw4km7xd", 1
API = "https://data.mendeley.com/public-api/datasets/{ds}"


def get(url):
    req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "mendeley_dl/1.0"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.load(r)


def list_all(ds, ver):
    base = API.format(ds=ds)
    folders = {"root": ""}
    try:
        for f in get(f"{base}/folders/{ver}"):
            folders[f["id"]] = f["name"]
        # vollständige Pfade aus parent_id bauen
        raw = {f["id"]: f for f in get(f"{base}/folders/{ver}")}

        def path(fid):
            f = raw.get(fid)
            if not f:
                return ""
            p = f.get("parent_id")
            return os.path.join(path(p), f["name"]) if p and p != "root" else f["name"]

        folders.update({fid: path(fid) for fid in raw})
    except Exception as e:
        print(f"(Ordnerliste nicht abrufbar: {e}; lade nur Root-Ebene)", file=sys.stderr)
    files = []
    for fid, fpath in folders.items():
        for f in get(f"{base}/files?folder_id={fid}&version={ver}"):
            cd = f.get("content_details", {})
            files.append({
                "path": os.path.join(fpath, f["filename"]),
                "size": cd.get("size") or f.get("size") or 0,
                "url": cd.get("download_url"),
            })
    return files


def download(url, dest):
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        print(f"  [da]  {dest}")
        return
    os.makedirs(os.path.dirname(dest) or ".", exist_ok=True)
    tmp = dest + ".part"
    req = urllib.request.Request(url, headers={"User-Agent": "mendeley_dl/1.0"})
    with urllib.request.urlopen(req, timeout=600) as r, open(tmp, "wb") as fh:
        while chunk := r.read(1 << 20):
            fh.write(chunk)
    os.replace(tmp, dest)
    print(f"  [ok]  {dest}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--out", default="data/01_mouse_brain_MALDI_visium")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--match", help="Regex auf Dateipfad")
    a = ap.parse_args()

    files = list_all(DS, VER)
    if a.match:
        rx = re.compile(a.match, re.I)
        files = [f for f in files if rx.search(f["path"])]
    total = sum(f["size"] for f in files)
    print(f"{len(files)} Dateien, {total/1e9:.2f} GB")
    for f in files:
        print(f"  {f['size']/1e6:9.1f} MB  {f['path']}")
    if a.list:
        return
    for f in files:
        if not f["url"]:
            print(f"  [keine URL] {f['path']}", file=sys.stderr)
            continue
        try:
            download(f["url"], os.path.join(a.out, f["path"]))
        except Exception as e:
            print(f"  [FEHLER] {f['path']}: {e}", file=sys.stderr)


if __name__ == "__main__":
    main()
