"""Our own CCC layer for all 27 Kuppe 2022 Visium sections: LIANA+ local bivariate
ligand-receptor scores, calibrated to the colleague's tables of 4 slides.

The colleague's scores (data/kuppe_2022/modalities/ccc_<run>.csv for 10X001, ACH0012, ACH0013,
ACH0023) exist only for those 4 slides; their exact settings are not documented. Calibration on
P9 IZ (ACH0012) against their table, LIANA+ 1.8.1 li.mt.bivariate, local cosine score,
Gaussian kernel, cutoff 0.1, log-normalised X of the Zenodo h5ad:
    bandwidth (x spot spacing)  0.4 / 0.6 / 0.75 / 0.855 / 1 / 2 / 3
    median per-pair Pearson r   0.20 / 0.77 / 0.84 / 0.83 / 0.82 / 0.67 / 0.58
-> bandwidth = 0.75 x the slide's median nearest-spot distance (in full-resolution pixels).
r = 0.84 is close but not identical: this is a reconstruction of the colleague's layer, not
their layer. The per-slide agreement on the 4 shared slides is written to compare_colleague.csv.

Per section:
  X = log-normalised expression shipped in processed1/*.h5ad (Seurat LogNormalize),
  spatial = pxl_col_fullres / pxl_row_fullres from visium_tables/<run>_spots.csv,
  li.ut.spatial_neighbors(bandwidth, cutoff=0.1, kernel="gaussian", set_diag=True),
  li.mt.bivariate(resource_name="consensus", local_name="cosine", global_name="morans",
                  n_perms=None, use_raw=False)  (LIANA+ default expression filters)
Outputs:
  data/kuppe_2022/modalities/ccc_all27_ours/ccc_<run>.csv  spots x "LIGAND^RECEPTOR" (0-1)
  results/kuppe_2022/ccc_all27_ours/liana_qc.csv, compare_colleague.csv
These per-slide tables are the input of 07_integration/ccc_harmony.py (folder mode).

Usage (env with liana >= 1.8, e.g. masld_atlas; from Derived_Modalities/):
    python 02_modalities/liana/liana_local_all27.py data/kuppe_2022 results/kuppe_2022/ccc_all27_ours [workers]
"""

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys  # noqa: E402
import time  # noqa: E402
from multiprocessing import get_context  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

DATA, RES = Path(sys.argv[1]), Path(sys.argv[2])
N_WORKERS = int(sys.argv[3]) if len(sys.argv) > 3 else 6
MOD = DATA / "modalities" / "ccc_all27_ours"
BW_FACTOR = 0.75


def run_one(job):
    import anndata as ad
    import liana as li
    from scipy.spatial import cKDTree

    f, run = job
    t0 = time.time()
    a = ad.read_h5ad(f)
    sp = pd.read_csv(DATA / "visium_tables" / f"{run}_spots.csv").set_index("barcode").loc[a.obs_names]
    a.obsm["spatial"] = sp[["pxl_col_fullres", "pxl_row_fullres"]].to_numpy(float)
    pitch = float(np.median(cKDTree(a.obsm["spatial"]).query(a.obsm["spatial"], k=2)[0][:, 1]))
    li.ut.spatial_neighbors(a, bandwidth=BW_FACTOR * pitch, cutoff=0.1, kernel="gaussian", set_diag=True)
    lr = li.mt.bivariate(a, resource_name="consensus", local_name="cosine", global_name="morans",
                         n_perms=None, use_raw=False, verbose=False)
    X = lr.X.toarray() if hasattr(lr.X, "toarray") else np.asarray(lr.X)
    df = pd.DataFrame(X.astype(np.float32), index=lr.obs_names, columns=lr.var_names)
    df.to_csv(MOD / f"ccc_{run}.csv", float_format="%.5g")
    return {"run": run, "spots": df.shape[0], "pairs": df.shape[1], "pitch_px": pitch,
            "bandwidth_px": BW_FACTOR * pitch, "seconds": round(time.time() - t0)}


def main():
    import anndata as ad

    MOD.mkdir(parents=True, exist_ok=True)
    RES.mkdir(parents=True, exist_ok=True)
    tables = {}
    for t in (DATA / "visium_tables").glob("*_spots.csv"):
        tables[t.name.removesuffix("_spots.csv")] = set(pd.read_csv(t, usecols=["barcode"]).barcode)
    jobs = []
    for f in sorted((DATA / "processed1").glob("Visium_*.h5ad")):
        names = set(ad.read_h5ad(f, backed="r").obs_names)
        run = [r for r, b in tables.items() if b == names]
        assert len(run) == 1, f"{f.name}: matches {run}"
        jobs.append((f, run[0]))
    print(f"{len(jobs)} sections; bandwidth {BW_FACTOR} x spot spacing", flush=True)
    rows = []
    with get_context("spawn").Pool(N_WORKERS) as pool:
        for r in pool.imap_unordered(run_one, jobs):
            print(r, flush=True)
            rows.append(r)
    pd.DataFrame(rows).sort_values("run").to_csv(RES / "liana_qc.csv", index=False)
    comp = []
    for run in ["10X001", "ACH0012", "ACH0013", "ACH0023"]:
        theirs = DATA / "modalities" / f"ccc_{run}.csv"
        if not theirs.exists():
            continue
        t = pd.read_csv(theirs, index_col=0)
        o = pd.read_csv(MOD / f"ccc_{run}.csv", index_col=0).loc[t.index]
        common = o.columns.intersection(t.columns)
        r = np.array([np.corrcoef(o[p], t[p])[0, 1] for p in common])
        comp.append({"run": run, "pairs_ours": o.shape[1], "pairs_colleague": t.shape[1], "shared": len(common),
                     "median_r": np.nanmedian(r), "p5_r": np.nanpercentile(r, 5),
                     "median_abs_diff": float(np.median(np.abs(o[common].values - t[common].values)))})
    pd.DataFrame(comp).to_csv(RES / "compare_colleague.csv", index=False)
    print(pd.DataFrame(comp).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
