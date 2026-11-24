"""Experiment 1: does the result depend on the Harmony covariates? (sections 1-4 of the plan)

Primary = the authors' covariates (patient + sample + batch); `none` and `batch` are sensitivity
variants. One PCA, corrected three ways; per variant a 15-NN graph, UMAP and Leiden at a fixed
K (9 = published cell-type niches, 12 = molecular niches), compared within K only. Every metric
is reported pooled and stratified by zone.

Input: the 27 processed Visium h5ads (raw counts in .raw), the authors' spot tables.
Spec name -> authors' column: visium_batch = batch, zone = major_labl, sample = hca_sample_id.

Expression recipe = exp0_sweep.preprocess on the joint object: genes shared by all 27 sections
(11,681), HVG seurat_v3 2,000 on counts with batch_key="sample", normalize_total, log1p,
scale(zero_center=False), PCA 30 (seed 0).

Metrics (sklearn ARI; harmonypy.compute_lisi, perplexity 30; sklearn silhouette, euclidean):
    ARI of clusters vs visium_batch, sample, zone (pooled only), celltype_niche, molecular_niche
    iLISI(visium_batch) scaled to [0, 1] (1 = fully mixed), cLISI(zone) scaled (1 = pure),
    silhouette(zone) on the embedding. Embedding metrics do not depend on K.
    Silhouette uses a stratified subsample (SIL_PER_ZONE spots per zone, seed 0).
Descriptive: silhouette-optimal K (silhouette of the Leiden labels on the authors' embedding,
K = 4..20), not used for the comparison.

Usage (from Derived_Modalities/):
    python 07_integration/exp1_harmony_covariates/exp1_covariates.py data/kuppe_2022 \
        results/kuppe_2022/exp1_harmony_covariates
"""
import sys
import time
from pathlib import Path

import anndata as ad
import harmonypy
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from sklearn.metrics import adjusted_rand_score as ari
from sklearn.metrics import silhouette_samples

from multiprocessing import get_context

sys.path.insert(0, str(Path(__file__).resolve().parent))
import exp1_workers as W  # noqa: E402

SEED = 0
KS = [9, 12]
VARIANTS = {"none": None, "batch": ["visium_batch"], "authors": ["patient", "sample", "visium_batch"]}
PRIMARY = "authors"
REN = {"batch": "visium_batch", "major_labl": "zone", "hca_sample_id": "sample"}
SIL_PER_ZONE = 3000
K_RANGE = range(4, 21)
ZONES = ["CTRL", "RZ", "BZ", "IZ", "FZ"]
N_PROC = 32  # cores of the allocation


def main():
    data, out = Path(sys.argv[1]), Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()


    def log(m):
        print(f"[{time.time() - t0:6.0f}s] {m}", flush=True)


    # ---- 1. joint object --------------------------------------------------------------------
    meta = pd.read_csv(data / "visium_raw" / "metadata-Visium.csv")
    parts = []
    for f in sorted((data / "processed1").glob("Visium_*.h5ad")):
        a = ad.read_h5ad(f)
        obs_names = a.obs_names
        raw = a.raw.to_adata()
        # the h5ad name does not identify the run (RZ_P6 has two), so match on the spot barcodes
        hits = [s for s in meta.hca_sample_id
                if (t := data / "visium_tables" / f"{s}_spots.csv").exists()
                and set(pd.read_csv(t, usecols=["barcode"]).barcode) == set(obs_names)]
        assert len(hits) == 1, f"{f.name}: spot tables matching its barcodes: {hits}"
        sample = hits[0]
        spots = pd.read_csv(data / "visium_tables" / f"{sample}_spots.csv").set_index("barcode").loc[obs_names]
        b = ad.AnnData(X=sp.csr_matrix(raw.X), var=pd.DataFrame(index=raw.var_names))
        b.obs = spots[["batch", "major_labl", "hca_sample_id", "patient", "patient_group",
                       "celltype_niche", "molecular_niche"]].rename(columns=REN).astype(str)
        b.obs_names = [f"{sample}_{x}" for x in obs_names]
        parts.append(b)
        log(f"{f.name} -> {sample}: {b.n_obs} spots, {b.n_vars} genes")
    assert len(parts) == 27
    a = ad.concat(parts, join="inner")
    a.layers["counts"] = a.X.copy()
    log(f"joint: {a.n_obs} spots x {a.n_vars} shared genes; zones {a.obs.zone.value_counts().to_dict()}")

    sc.pp.highly_variable_genes(a, n_top_genes=2000, flavor="seurat_v3", layer="counts", batch_key="sample")
    sc.pp.normalize_total(a)
    sc.pp.log1p(a)
    sc.pp.scale(a, zero_center=False)
    sc.pp.pca(a, n_comps=30, mask_var="highly_variable", random_state=SEED)
    log("PCA done (computed once for all variants)")

    # ---- 2. variants: Harmony on the same PCs -------------------------------------------------
    emb = {"none": a.obsm["X_pca"]}
    for v, cov in VARIANTS.items():
        if cov is None:
            continue
        ho = harmonypy.run_harmony(a.obsm["X_pca"], a.obs, cov, random_state=SEED, verbose=False)
        z = np.asarray(ho.Z_corr)
        emb[v] = z.T if z.shape[0] != a.n_obs else z
        log(f"Harmony {v} {cov} done")
    np.savez_compressed(out / "embeddings.npz", obs_names=np.asarray(a.obs_names), **emb)

    # ---- clustering: neighbours per variant, then every Leiden search and UMAP in parallel -------
    obs = a.obs
    rng = np.random.default_rng(SEED)
    sil_idx = np.concatenate([rng.choice(np.where(obs.zone == zn)[0], min(SIL_PER_ZONE, (obs.zone == zn).sum()),
                                         replace=False) for zn in ZONES])
    graphs = {}
    for v, z in emb.items():
        a.obsm["X_v"] = z
        sc.pp.neighbors(a, n_neighbors=15, use_rep="X_v", random_state=SEED, key_added=v)
        graphs[v] = {"connectivities": a.obsp[f"{v}_connectivities"], "distances": a.obsp[f"{v}_distances"],
                     "uns": dict(a.uns[v]), "emb": z}
        log(f"neighbours {v} done")
    jobs = sorted({(v, k) for v in emb for k in KS} | {(PRIMARY, k) for k in K_RANGE}, key=lambda j: -j[1])
    n_proc = min(N_PROC, len(jobs) + len(emb))
    log(f"{len(jobs)} Leiden searches + {len(emb)} UMAPs on {n_proc} processes")
    lab_by, res_rows, krows = {}, [], []
    with get_context("spawn").Pool(n_proc, initializer=W.init,
                                   initargs=(graphs, a.n_obs, sil_idx, {v: z[sil_idx] for v, z in emb.items()})) as pool:
        um = pool.map_async(W.umap, [(v, SEED) for v in emb])
        for v, k, r, ok, lab, sil in pool.imap_unordered(W.leiden_k, [(v, k, SEED) for v, k in jobs]):
            lab_by[(v, k)] = lab
            if k in KS:
                res_rows.append({"variant": v, "K": k, "resolution": r, "reached_K": ok,
                                 "n_clusters": len(set(lab)), "silhouette_clusters": sil})
            if v == PRIMARY:
                krows.append({"K": k, "reached_K": ok, "resolution": r, "silhouette": sil})
            log(f"{v} K={k}: resolution {r:.4f}, reached {ok}, silhouette(clusters) {sil:.4f}")
        umaps = dict(um.get())
        log("UMAPs done")
    labels = pd.DataFrame({f"{v}_K{k}": lab_by[(v, k)] for v in emb for k in KS}, index=a.obs_names)
    labels = pd.concat([labels, a.obs], axis=1)
    labels.to_csv(out / "exp1_labels.csv.gz")
    np.savez_compressed(out / "umaps.npz", **umaps)
    kdf = pd.DataFrame(krows).sort_values("K")
    kdf.to_csv(out / "exp1_silhouette_optimal_k.csv", index=False)
    log(f"silhouette-optimal K on {PRIMARY}: {int(kdf.loc[kdf.silhouette.idxmax(), 'K'])}")

    # ---- 3. metrics -------------------------------------------------------------------------
    emb_metrics = {}
    for v, z in emb.items():
        lisi = harmonypy.compute_lisi(z, obs[["visium_batch", "zone"]].reset_index(drop=True), ["visium_batch", "zone"])
        sil = np.full(a.n_obs, np.nan)
        sil[sil_idx] = silhouette_samples(z[sil_idx], obs.zone.to_numpy()[sil_idx])
        emb_metrics[v] = pd.DataFrame({"iLISI_batch": (lisi[:, 0] - 1) / (2 - 1),
                                       "cLISI_zone": (len(ZONES) - lisi[:, 1]) / (len(ZONES) - 1),
                                       "silhouette_zone": sil}, index=obs.index)
        log(f"embedding metrics {v} done")

    rows = []
    for v in emb:
        for k in KS:
            lab = labels[f"{v}_K{k}"]
            for zn in ["all"] + ZONES:
                m = np.ones(len(obs), bool) if zn == "all" else (obs.zone == zn).to_numpy()
                e = emb_metrics[v][m]
                o = obs[m]
                rows.append({
                    "variant": v, "primary": v == PRIMARY, "K": k, "zone": zn, "n_spots": int(m.sum()),
                    "n_10X_spots": int((o.visium_batch == "10X").sum()),
                    "ARI_batch": ari(o.visium_batch, lab[m]) if o.visium_batch.nunique() > 1 else np.nan,
                    "ARI_sample": ari(o["sample"], lab[m]),
                    "ARI_zone": ari(o.zone, lab[m]) if zn == "all" else np.nan,
                    "ARI_celltype_niche": ari(o.celltype_niche, lab[m]),
                    "ARI_molecular_niche": ari(o.molecular_niche, lab[m]),
                    "iLISI_batch": e.iLISI_batch.mean(),
                    "cLISI_zone": e.cLISI_zone.mean(),
                    "silhouette_zone": e.silhouette_zone.mean(),
                    "n_clusters_in_zone": lab[m].nunique(),
                })
    comp = pd.DataFrame(rows).merge(pd.DataFrame(res_rows), on=["variant", "K"])
    comp.to_csv(out / "exp1_covariate_comparison.csv", index=False)
    log("metrics written")

    # ---- 4. figures -------------------------------------------------------------------------
    COL = {"none": "#8c8c8c", "batch": "#4c78a8", "authors": "#e45756"}
    pooled = comp[comp.zone == "all"]
    per_zone = comp[comp.zone != "all"]
    fig, ax = plt.subplots(2, 2, figsize=(11, 9))
    for j, k in enumerate(KS):
        p = pooled[pooled.K == k]
        for _, r in p.iterrows():
            ax[0, j].scatter(1 - r.ARI_batch, r.ARI_zone, s=120, color=COL[r.variant], edgecolor="k",
                             marker="*" if r.primary else "o", label=r.variant, zorder=3)
            ax[0, j].annotate(f"{r.variant}\niLISI {r.iLISI_batch:.2f}", (1 - r.ARI_batch, r.ARI_zone),
                              xytext=(6, 6), textcoords="offset points", fontsize=8)
        ax[0, j].set(title=f"K = {k}: pooled", xlabel="batch removal: 1 - ARI(clusters, visium_batch)",
                     ylabel="biology: ARI(clusters, zone)")
        ax[0, j].grid(alpha=0.3)
    x = np.arange(len(ZONES))
    w = 0.27
    for i, v in enumerate(emb):
        pz = per_zone[(per_zone.variant == v) & (per_zone.K == KS[0])].set_index("zone").loc[ZONES]
        ax[1, 0].bar(x + (i - 1) * w, pz.silhouette_zone, w, color=COL[v], label=v)
        ax[1, 1].bar(x + (i - 1) * w, pz.iLISI_batch, w, color=COL[v], label=v)
    for axx, t in [(ax[1, 0], "biology by zone: silhouette(zone)"), (ax[1, 1], "batch mixing by zone: iLISI(visium_batch)")]:
        axx.set_xticks(x, [f"{z}\n{int(per_zone[(per_zone.zone == z)].n_10X_spots.iloc[0] / per_zone[(per_zone.zone == z)].n_spots.iloc[0] * 100)}% 10X" for z in ZONES])
        axx.set_title(t)
        axx.axhline(0, color="k", lw=0.5)
        axx.legend(fontsize=8)
    fig.suptitle("Exp. 1: batch removal vs biology preservation (star = primary, authors' covariates)")
    fig.tight_layout()
    fig.savefig(out / "fig1_tradeoff.png", dpi=150)
    fig.savefig(out / "fig1_tradeoff.pdf")

    colour_by = ["visium_batch", "zone", "celltype_niche"]
    fig, ax = plt.subplots(len(emb), len(colour_by), figsize=(15, 5 * len(emb)))
    order = rng.permutation(a.n_obs)
    for i, v in enumerate(emb):
        u = umaps[v][order]
        for j, c in enumerate(colour_by):
            cats = sorted(obs[c].unique())
            cmap = plt.get_cmap("tab10" if len(cats) <= 10 else "tab20")
            cidx = pd.Categorical(obs[c].to_numpy()[order], categories=cats).codes
            ax[i, j].scatter(u[:, 0], u[:, 1], c=[cmap(q % cmap.N) for q in cidx], s=0.3, rasterized=True)
            ax[i, j].set(xticks=[], yticks=[], title=f"{v}: {c}")
            if i == 0:
                ax[i, j].legend(handles=[plt.Line2D([], [], ls="", marker="o", color=cmap(q % cmap.N), label=cat)
                                         for q, cat in enumerate(cats)], fontsize=7, markerscale=1, loc="best")
    fig.tight_layout()
    fig.savefig(out / "fig2_umaps.png", dpi=130)

    fig, ax = plt.subplots(len(emb), len(KS), figsize=(13, 3.6 * len(emb)))
    zc = dict(zip(ZONES, plt.get_cmap("Set2").colors))
    for i, v in enumerate(emb):
        for j, k in enumerate(KS):
            ct = pd.crosstab(labels[f"{v}_K{k}"], obs.zone, normalize="index")[ZONES]
            ct = ct.loc[sorted(ct.index, key=int)]
            bottom = np.zeros(len(ct))
            for zn in ZONES:
                ax[i, j].bar(range(len(ct)), ct[zn], bottom=bottom, color=zc[zn], label=zn)
                bottom += ct[zn].to_numpy()
            ax[i, j].set(title=f"{v}, K = {k}", xticks=range(len(ct)), xticklabels=ct.index, ylabel="share of spots")
            if i == 0 and j == 0:
                ax[i, j].legend(fontsize=8, ncol=5, loc="upper center", bbox_to_anchor=(0.5, 1.25))
    fig.tight_layout()
    fig.savefig(out / "fig3_cluster_zone_composition.png", dpi=150)
    log("figures written; done")


if __name__ == "__main__":
    main()
