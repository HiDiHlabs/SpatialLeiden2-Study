"""One MuData per slide with every matrix we inferred, so later work starts from the objects
instead of re-running the pipelines.

data/kuppe_2022/mudata/kuppe_<run>_derived.h5mu, modalities (same spots, same order):
    rna      Visium counts (layers["counts"]; X = log1p CP10k), HVG flags, PCA 30, 15-NN graph,
             as built by experiment 0 (exp0_sweep.preprocess, but X kept log-normalised)
    tf       CollecTRI ULM TF activity (X = score, layers["padj"])
    ccc      team CCC local ligand-receptor scores (0-1; spatially smoothed upstream)
    met      scCellFie metabolic task activity
    spagrn   SpaGRN regulon AUCell scores; uns["regulons"] = {regulon: [targets]} (JSON string)
mdata.obs: the authors' spot table (niches, cell2location proportions, QC) and our cluster labels
(seed 0, K = authors' cell-type niche count):
    sl_rna, sl_rna_space            RNA, spatial 0 / 0.75              (exp0_k_celltype_niche)
    sl_tf, sl_tf_space              TF,  spatial 0 / 0.75              (exp0tf_k_celltype_niche)
    sl_ccc, sl_ccc_space025/075     CCC, spatial 0 / 0.25 / 0.75       (exp0ccc_k_celltype_niche)
    sl_met, sl_met_space            metabolic, spatial 0 / 0.75        (exp0met_k_celltype_niche)
    sl_rna_tf_space                 RNA 1, TF 0.5, spatial 0.75         (mm_proto_rna_tf)
    sl_rna_ccc_space025/075         RNA 1, CCC 0.5, spatial 0.25/0.75   (mm_configs_sp0.25 / mm_configs)
    sl_rna_met_space, sl_rna_all_space                                  (mm_configs)
    vessel_k2, vessel_k3            P9 rep2 only: nested sub-clusters of the artery band
mdata.obsm["spatial"]: full-resolution pixel coordinates. All five seeds stay in results/.

Usage (env spatialleiden, from Derived_Modalities/): python 02_modalities/build_mudata.py data/kuppe_2022 results/kuppe_2022
"""

import importlib.util
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
spec = importlib.util.spec_from_file_location(
    "exp0_sweep", HERE.parents[1] / "03_single_modality_clustering" / "exp0_layer_ratio_sweep" / "exp0_sweep.py")
e0 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(e0)

import anndata as ad  # noqa: E402
import mudata as md  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402

DATA, RES = Path(sys.argv[1]), Path(sys.argv[2])
OUT = DATA / "mudata"
OUT.mkdir(parents=True, exist_ok=True)
LABELS = {"sl_rna": ("exp0_k_celltype_niche", "r0_s0"), "sl_rna_space": ("exp0_k_celltype_niche", "r0.75_s0"),
          "sl_tf": ("exp0tf_k_celltype_niche", "r0_s0"), "sl_tf_space": ("exp0tf_k_celltype_niche", "r0.75_s0"),
          "sl_ccc": ("exp0ccc_k_celltype_niche", "r0_s0"), "sl_ccc_space025": ("exp0ccc_k_celltype_niche", "r0.25_s0"),
          "sl_ccc_space075": ("exp0ccc_k_celltype_niche", "r0.75_s0"),
          "sl_met": ("exp0met_k_celltype_niche", "r0_s0"), "sl_met_space": ("exp0met_k_celltype_niche", "r0.75_s0"),
          "sl_rna_tf_space": ("mm_proto_rna_tf", "sp0.75_tf0.5_s0"),
          "sl_rna_ccc_space025": ("mm_configs_sp0.25", "rna_ccc_s0"), "sl_rna_ccc_space075": ("mm_configs", "rna_ccc_s0"),
          "sl_rna_met_space": ("mm_configs", "rna_met_s0"), "sl_rna_all_space": ("mm_configs", "rna_all_s0")}


def bare(a, X, var_names):
    return ad.AnnData(X=np.asarray(X, dtype=np.float32), obs=pd.DataFrame(index=a.obs_names),
                      var=pd.DataFrame(index=[str(v) for v in var_names]))


for run in ["10X001", "ACH0012", "ACH0013"]:
    rna = e0.load_slide(DATA, run)
    sc.pp.filter_genes(rna, min_cells=10)
    sc.pp.highly_variable_genes(rna, n_top_genes=2000, flavor="seurat_v3", layer="counts")
    sc.pp.normalize_total(rna)
    sc.pp.log1p(rna)
    pca_src = rna.copy()
    sc.pp.scale(pca_src, zero_center=False)
    sc.pp.pca(pca_src, n_comps=30, mask_var="highly_variable", random_state=e0.GRAPH_SEED)
    sc.pp.neighbors(pca_src, n_neighbors=15, n_pcs=30, random_state=e0.GRAPH_SEED)
    rna.obsm["X_pca"] = pca_src.obsm["X_pca"]
    rna.obsp["connectivities"] = pca_src.obsp["connectivities"]
    rna.obsp["distances"] = pca_src.obsp["distances"]
    names = rna.obs_names
    tf = ad.read_h5ad(DATA / "modalities" / f"collectri_ulm_{run}.h5ad")[names]
    tfm = bare(tf, tf.X, tf.var_names)
    tfm.layers["padj"] = np.asarray(tf.layers["padj"], dtype=np.float32)
    ccc = pd.read_csv(DATA / "modalities" / f"ccc_{run}.csv", index_col=0).loc[names]
    met = ad.read_h5ad(DATA / "modalities" / f"sccellfie_tasks_{run}.h5ad")[names]
    g = ad.read_h5ad(RES / "spagrn" / run / f"{run}_spagrn.h5ad")
    auc = g.obsm["auc_mtx"].loc[names]
    sp = bare(g[names], auc.to_numpy(), auc.columns)
    sp.uns["regulons"] = json.dumps({k: sorted(map(str, v)) for k, v in g.uns["regulon_dict"].items()})
    obs = rna.obs.copy()
    spots = pd.read_csv(DATA / "visium_tables" / f"{run}_spots.csv").set_index("barcode").loc[names]
    for c in spots.columns:
        if c not in obs.columns:
            obs[c] = spots[c].to_numpy()
    rna_slim = ad.AnnData(X=rna.X, obs=pd.DataFrame(index=names), var=rna.var, layers={"counts": rna.layers["counts"]},
                          obsm={"X_pca": rna.obsm["X_pca"]},
                          obsp={"connectivities": rna.obsp["connectivities"], "distances": rna.obsp["distances"]})
    mods = {"rna": rna_slim, "tf": tfm, "ccc": bare(rna, ccc.to_numpy(), ccc.columns),
            "met": bare(met, met.X.toarray() if hasattr(met.X, "toarray") else met.X, met.var_names), "spagrn": sp}
    m = md.MuData(mods)
    for c in obs.columns:
        m.obs[c] = obs[c].astype(str) if obs[c].dtype == object else obs[c]
    for key, (folder, col) in LABELS.items():
        lab = pd.read_csv(RES / folder / f"labels_{run}.csv.gz", index_col=0, dtype=str)[col]
        assert list(lab.index) == list(names), f"{run} {key}: spot order"
        m.obs[key] = pd.Categorical(lab.to_numpy())
    if run == "ACH0013":
        subl = pd.read_csv(RES / "mm_proto_rna_tf" / "vessel_band" / "vessel_band_sublabels.csv.gz", index_col=0, dtype=str)
        for k in ("k2", "k3"):
            m.obs[f"vessel_{k}"] = pd.Categorical(subl[f"{k}_s0"].reindex(names).fillna("outside").to_numpy())
    m.obsm["spatial"] = rna.obsm["spatial"]
    m.uns["description"] = __doc__
    m.write(OUT / f"kuppe_{run}_derived.h5mu")
    print(run, {k: v.shape for k, v in m.mod.items()}, f"obs columns {m.obs.shape[1]}", flush=True)
