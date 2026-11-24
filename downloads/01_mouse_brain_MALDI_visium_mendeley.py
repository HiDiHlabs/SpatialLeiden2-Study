#!/usr/bin/env python3
"""
COSIE-Datensatz 1: Mouse Brain, Spatial MALDI-Visium (SMA), Parkinson-Modell (6-OHDA).
H&E + Metabolom (MALDI-MSI) + RNA (Visium). Vicari et al., Nat Biotechnol 2024.
Mendeley Data, DOI 10.17632/w7nw4km7xd.1 (6 Dateien, ca. 6.2 GB).

"""
import argparse
import json
import os
import re
import sys
import time
import urllib.request
import zipfile
from pathlib import Path

DS, VER = "w7nw4km7xd", 1
API = "https://data.mendeley.com/public-api/datasets/{ds}"
HEADERS = {"User-Agent": "spatialleiden-study-download/1.0"}
DEFAULT_OUT = Path(__file__).resolve().parent.parent / "data" / "01_mouse_brain_MALDI_visium"


def get_json(url, retries=3):
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers={**HEADERS, "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.load(r)
        except Exception:
            if attempt == retries:
                raise
            time.sleep(5 * attempt)


def list_files(ds, ver):
    """Alle Dateien des Datensatzes mit Pfad, Groesse und Download-URL (inkl. Unterordner)."""
    base = API.format(ds=ds)
    try:
        raw = {f["id"]: f for f in get_json(f"{base}/folders/{ver}")}
    except Exception as e:
        print(f"(Ordnerliste nicht abrufbar: {e}; lade nur die oberste Ebene)", file=sys.stderr)
        raw = {}

    def folder_path(fid):
        f = raw.get(fid)
        if not f:
            return ""
        parent = f.get("parent_id")
        return os.path.join(folder_path(parent), f["name"]) if parent and parent != "root" else f["name"]

    folders = {"root": "", **{fid: folder_path(fid) for fid in raw}}
    files = []
    for fid, fpath in folders.items():
        for f in get_json(f"{base}/files?folder_id={fid}&version={ver}"):
            cd = f.get("content_details", {})
            files.append({
                "path": os.path.join(fpath, f["filename"]),
                "size": int(cd.get("size") or f.get("size") or 0),
                "url": cd.get("download_url"),
            })
    return sorted(files, key=lambda f: f["path"])


def download(url, dest, size, retries=3):
    """Laedt url nach dest (ueber dest.part); ueberspringt, wenn dest schon mit der erwarteten Groesse da ist."""
    if dest.exists() and (size == 0 or dest.stat().st_size == size):
        print(f"  [da]     {dest.name}")
        return
    if dest.exists():
        print(f"  [neu]    {dest.name}: Groesse {dest.stat().st_size:,} statt {size:,} Bytes, lade neu")
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(dest.name + ".part")
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers=HEADERS)
            with urllib.request.urlopen(req, timeout=600) as r, open(tmp, "wb") as fh:
                done, next_report = 0, 0.1
                while chunk := r.read(1 << 20):
                    fh.write(chunk)
                    done += len(chunk)
                    if size > 100e6 and done / size >= next_report:  # progress only for large files
                        print(f"           {dest.name}: {done / size:4.0%} ({done / 1e9:.2f} / {size / 1e9:.2f} GB)", flush=True)
                        next_report += 0.1
            if size and tmp.stat().st_size != size:
                raise IOError(f"unvollstaendig: {tmp.stat().st_size:,} statt {size:,} Bytes")
            os.replace(tmp, dest)
            print(f"  [ok]     {dest.name}")
            return
        except Exception as e:
            tmp.unlink(missing_ok=True)
            if attempt == retries:
                raise
            print(f"  [retry]  {dest.name}: {e} (Versuch {attempt}/{retries})", file=sys.stderr)
            time.sleep(10 * attempt)


def extract(zip_path):
    """Entpackt ein Zip neben sich, ohne macOS-Metadaten (__MACOSX, ._*); ueberspringt bereits entpackte Archive."""
    marker = zip_path.with_name(f".{zip_path.stem}.extracted")
    with zipfile.ZipFile(zip_path) as zf:
        members = [m for m in zf.namelist() if not m.startswith("__MACOSX/") and not Path(m).name.startswith("._")]
        tops = {m.split("/")[0] for m in members}
        if marker.exists() or all((zip_path.parent / t).exists() for t in tops):
            print(f"  [da]     {zip_path.stem}/ (bereits entpackt)")
            marker.touch()
            return
        zf.extractall(zip_path.parent, members=members)
    marker.touch()
    print(f"  [entp.]  {zip_path.name} -> {', '.join(sorted(t + '/' for t in tops))}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--out", type=Path, default=DEFAULT_OUT, help=f"Zielordner (Standard: {DEFAULT_OUT})")
    ap.add_argument("--list", action="store_true", help="nur die Dateien anzeigen, nichts laden")
    ap.add_argument("--match", help='nur Dateien, deren Pfad zu dieser Regex passt, z. B. "metadata|visium|misc"')
    ap.add_argument("--no-extract", action="store_true", help="Zip-Archive nicht entpacken")
    ap.add_argument("--delete-zip", action="store_true", help="Zip-Archive nach dem Entpacken loeschen")
    a = ap.parse_args()

    files = list_files(DS, VER)
    if a.match:
        rx = re.compile(a.match, re.I)
        files = [f for f in files if rx.search(f["path"])]
    print(f"{len(files)} Dateien, {sum(f['size'] for f in files) / 1e9:.2f} GB -> {a.out}")
    for f in files:
        print(f"  {f['size'] / 1e6:9.1f} MB  {f['path']}")
    if a.list or not files:
        return

    errors = []
    for f in files:
        dest = a.out / f["path"]
        if not f["url"]:
            errors.append(f"{f['path']}: keine Download-URL")
            continue
        try:
            download(f["url"], dest, f["size"])
            if dest.suffix == ".zip" and not a.no_extract:
                extract(dest)
                if a.delete_zip:
                    dest.unlink()
                    print(f"  [del]    {dest.name}")
        except Exception as e:
            errors.append(f"{f['path']}: {e}")

    if errors:
        print("\nFehler:", *errors, sep="\n  ", file=sys.stderr)
        sys.exit(1)
    print("\nFertig.")


if __name__ == "__main__":
    main()
