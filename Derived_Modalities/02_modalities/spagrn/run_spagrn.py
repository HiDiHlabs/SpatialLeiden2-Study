"""Run SpaGRN (spatial gene regulatory networks) on one Kuppe 2022 Visium slide.

SpaGRN infers TF -> target modules from spatial co-expression (bivariate Moran's I on a
spatial kNN graph), prunes them with cisTarget motif enrichment (pySCENIC), and scores
each regulon per spot with AUCell. Settings follow the authors' CLI defaults
(spagrn/cli/spagrn_parser.py, v1.1.4): model "danb" on raw counts, 30 spatial
neighbours, mode "moran", gene/spot filter from InferNetwork.preprocess.

Deviation from the defaults, on purpose:
  * preprocess(max_gene_num=...) is raised from 4,000 to effectively no limit. The
    default silently drops every spot with >= 4,000 detected genes, which would remove
    much of the deep P9 rep2 section and break the one-to-one spot match with the RNA
    and CollecTRI layers. All spots are kept.

Inputs: Zenodo h5ad raw counts (.raw) joined by barcode to the authors' spot table
(true coordinates, niches). cisTarget hg38 v10 databases and the TF list on the vault.

Usage (env spagrn, from Derived_Modalities/):
    python 02_modalities/spagrn/run_spagrn.py <run_id> data/kuppe_2022 results/kuppe_2022/spagrn --workers 8
Output: results/kuppe_2022/spagrn/<run_id>/<run_id>_spagrn.h5ad (+ tmp_files/ with the
regulons, motifs, AUC matrix, selected genes).
"""

import argparse
import os
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

CISTARGET = os.environ.get("CISTARGET_DIR") or sys.exit("set CISTARGET_DIR to the folder with the cisTarget hg38 v10 files (resources.aertslab.org/cistarget)")
DATABASES = f"{CISTARGET}/hg38_*_full_tx_v10_clust.genes_vs_motifs.rankings.feather"  # glob, 2 databases
MOTIF_ANNO = f"{CISTARGET}/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
TFS = f"{CISTARGET}/allTFs_hg38.txt"
H5AD = {"10X001": "Visium_control_P1.h5ad", "ACH0013": "Visium_GT_IZ_P9_rep2.h5ad",
        "ACH0012": "Visium_GT_IZ_P9.h5ad"}


def build(data: Path, run: str) -> ad.AnnData:
    h5 = ad.read_h5ad(data / "processed1" / H5AD[run])
    spots = pd.read_csv(data / "visium_tables" / f"{run}_spots.csv").set_index("barcode")
    assert set(spots.index) == set(h5.obs_names), "barcodes differ"
    spots = spots.loc[h5.obs_names]
    a = ad.AnnData(X=h5.raw.X.copy(), obs=spots.copy(), var=pd.DataFrame(index=h5.raw.var_names))
    assert np.allclose(np.asarray(a.X.sum(1)).ravel(), spots.nCount_Spatial.to_numpy()), "counts mismatch"
    a.obsm["spatial"] = spots[["pxl_col_fullres", "pxl_row_fullres"]].to_numpy(float)
    a.obs["celltype_niche"] = a.obs["celltype_niche"].astype(str)
    return a


if __name__ == "__main__":  # required by SpaGRN/dask multiprocessing
    p = argparse.ArgumentParser()
    p.add_argument("run", choices=list(H5AD))
    p.add_argument("data", type=Path)
    p.add_argument("out", type=Path)
    p.add_argument("--workers", type=int, default=8)
    args = p.parse_args()

    from spagrn.regulatory_network import InferNetwork as irn

    t0 = time.time()
    out = args.out / args.run
    out.mkdir(parents=True, exist_ok=True)
    a = build(args.data, args.run)
    n_spots = a.n_obs
    a = irn.preprocess(a, max_gene_num=10**9)  # keep every spot (see module docstring)
    assert a.n_obs == n_spots, f"preprocess dropped spots: {n_spots} -> {a.n_obs}"
    a = a.copy()
    print(f"{args.run}: {a.n_obs} spots, {a.n_vars} genes after SpaGRN gene filter", flush=True)

    grn = irn(a, project_name=args.run)
    grn.infer(DATABASES, MOTIF_ANNO, TFS,
              niche_df=None,
              num_workers=args.workers,
              cache=False,
              output_dir=str(out),
              save_tmp=True,  # required: infer() builds paths from tmp_dir
              layers=None,    # X = raw counts
              latent_obsm_key="spatial",
              model="danb",
              n_neighbors=30,
              mode="moran",
              cluster_label="celltype_niche",
              noweights=False)
    print(f"{args.run}: done in {(time.time() - t0) / 60:.1f} min", flush=True)
