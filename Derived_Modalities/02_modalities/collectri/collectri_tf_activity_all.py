"""CollecTRI TF activity (decoupler ULM) on all 27 Visium sections of Kuppe et al. 2022.

Same recipe as collectri_tf_activity.py (which covers the three experiment-0 slides): raw counts
from the Zenodo h5ad, normalize_total 1e4 + log1p on all genes, dc.mt.ulm with tmin = 5, CollecTRI
human from dc.op.collectri. All 27 sections are scored with one download of the network, saved
next to the results, so every slide uses the same network. The run of each section is found by
matching its barcodes to the authors' spot tables (the h5ad name does not identify RZ_P6's run).

Outputs (separate folder, so the 23 Nov files of the three slides are not overwritten):
  data/kuppe_2022/modalities/collectri_all27/collectri_ulm_<run>.h5ad   spots x TFs, layers["padj"]
  results/kuppe_2022/collectri_all27/collectri_qc.csv, collectri_network_used.csv.gz
  results/kuppe_2022/collectri_all27/compare_23nov.csv  agreement with the 23 Nov files (3 slides)

Usage (env decoupler, from Derived_Modalities/):
    python 02_modalities/collectri/collectri_tf_activity_all.py data/kuppe_2022 results/kuppe_2022/collectri_all27
"""
import sys
from pathlib import Path

import anndata as ad
import decoupler as dc
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

DATA, RES = Path(sys.argv[1]), Path(sys.argv[2])
MOD = DATA / "modalities" / "collectri_all27"
for d in (MOD, RES):
    d.mkdir(parents=True, exist_ok=True)
TMIN = 5

print("decoupler", dc.__version__, flush=True)
net = dc.op.collectri(organism="human")
net.to_csv(RES / "collectri_network_used.csv.gz", index=False)
print(f"CollecTRI: {net.source.nunique()} TFs, {len(net)} edges", flush=True)

meta = pd.read_csv(DATA / "visium_raw" / "metadata-Visium.csv")
tables = {r: set(pd.read_csv(DATA / "visium_tables" / f"{r}_spots.csv", usecols=["barcode"]).barcode)
          for r in meta.hca_sample_id if (DATA / "visium_tables" / f"{r}_spots.csv").exists()}

qc = []
for f in sorted((DATA / "processed1").glob("Visium_*.h5ad")):
    h5 = ad.read_h5ad(f)
    run = [r for r, bc in tables.items() if bc == set(h5.obs_names)]
    assert len(run) == 1, f"{f.name}: {run}"
    run = run[0]
    spots = pd.read_csv(DATA / "visium_tables" / f"{run}_spots.csv").set_index("barcode").loc[h5.obs_names]
    a = ad.AnnData(X=sparse.csr_matrix(h5.raw.X), obs=spots.copy(), var=pd.DataFrame(index=h5.raw.var_names))
    a.obsm["spatial"] = spots[["pxl_col_fullres", "pxl_row_fullres"]].to_numpy(float)
    ngenes = np.asarray((a.X > 0).sum(1)).ravel()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    dc.mt.ulm(data=a, net=net, tmin=TMIN)
    scores = dc.pp.get_obsm(a, key="score_ulm")
    out = ad.AnnData(X=scores.X.astype(np.float32), obs=a.obs.copy(), var=pd.DataFrame(index=scores.var_names))
    out.layers["padj"] = np.asarray(a.obsm["padj_ulm"].loc[out.obs_names, out.var_names], dtype=np.float32)
    out.obsm["spatial"] = a.obsm["spatial"]
    out.uns["collectri_ulm"] = {"decoupler_version": dc.__version__, "method": "ulm", "tmin": TMIN,
                                "network": "CollecTRI human (dc.op.collectri), collectri_network_used.csv.gz",
                                "normalization": "normalize_total 1e4 + log1p, all genes in the file",
                                "source_h5ad": f.name, "n_genes_in_slide": int(a.n_vars)}
    out.write_h5ad(MOD / f"collectri_ulm_{run}.h5ad", compression="gzip")
    qc.append({"run": run, "h5ad": f.name, "spots": out.n_obs, "genes_in_file": a.n_vars,
               "median_genes_per_spot": float(np.median(ngenes)), "tfs_scored": out.n_vars,
               "frac_scores_padj_lt_0.05": float((out.layers["padj"] < 0.05).mean())})
    print(qc[-1], flush=True)

qc = pd.DataFrame(qc)
common = set.intersection(*(set(ad.read_h5ad(MOD / f"collectri_ulm_{r}.h5ad", backed="r").var_names) for r in qc.run))
qc["tfs_scored_in_all_27"] = len(common)
qc.to_csv(RES / "collectri_qc.csv", index=False)
print(f"TFs scored in all 27 sections: {len(common)}")

# agreement with the 23 Nov outputs of the three experiment-0 slides (network may have been updated)
rows = []
for run in ["10X001", "ACH0012", "ACH0013"]:
    old_f = DATA / "modalities" / f"collectri_ulm_{run}.h5ad"
    if not old_f.exists():
        continue
    old, new = ad.read_h5ad(old_f), ad.read_h5ad(MOD / f"collectri_ulm_{run}.h5ad")
    tf = old.var_names.intersection(new.var_names)
    o, n = old[new.obs_names, tf].X, new[:, tf].X
    r = [np.corrcoef(o[:, j], n[:, j])[0, 1] for j in range(len(tf)) if o[:, j].std() > 0 and n[:, j].std() > 0]
    rows.append({"run": run, "tfs_23nov": old.n_vars, "tfs_now": new.n_vars, "tfs_shared": len(tf),
                 "max_abs_diff": float(np.abs(o - n).max()), "median_per_tf_pearson": float(np.median(r)),
                 "min_per_tf_pearson": float(np.min(r))})
pd.DataFrame(rows).to_csv(RES / "compare_23nov.csv", index=False)
print(pd.DataFrame(rows).to_string(index=False))
print("done")
