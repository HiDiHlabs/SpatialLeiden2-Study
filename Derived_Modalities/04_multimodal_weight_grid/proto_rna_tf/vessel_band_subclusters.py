"""Is there internal structure in the P9 rep2 vessel band (domain 6 of RNA + TF + space)?

Sub-clusters the 266 band spots with SpatialLeiden (RNA + TF + space, same weights as the
parent run: RNA 1, TF 0.5, spatial 0.75) at K = 2 and 3, 5 seeds. The RNA and TF PCA come
from the whole slide (experiment-0 preprocessing); kNN graphs (10 neighbours) are rebuilt on
the band spots only, the spatial graph is the Visium grid restricted to the band.
Then asks what a split follows:
    - position across the band (distance to the band edge, in grid steps: core vs rim)
    - position along the band (first principal axis of the band coordinates)
    - cell2location composition and top genes per sub-cluster
Seed stability (mean pairwise ARI) tells whether a split is reproducible at all.

Usage (env spatialleiden, from Derived_Modalities/):
    python 04_multimodal_weight_grid/proto_rna_tf/vessel_band_subclusters.py data/kuppe_2022 results/kuppe_2022
"""

import importlib.util
import sys
from itertools import combinations
from pathlib import Path

HERE = Path(__file__).resolve()
spec = importlib.util.spec_from_file_location(
    "exp0_sweep", HERE.parents[2] / "03_single_modality_clustering" / "exp0_layer_ratio_sweep" / "exp0_sweep.py")
e0 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(e0)

import anndata as ad  # noqa: E402
import mudata as md  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
import spatialleiden as sl  # noqa: E402
import squidpy as sq  # noqa: E402
from sklearn.metrics import adjusted_rand_score  # noqa: E402

DATA, RES = Path(sys.argv[1]), Path(sys.argv[2])
OUT = RES / "mm_proto_rna_tf" / "vessel_band"
OUT.mkdir(parents=True, exist_ok=True)
RUN, KEY, BAND = "ACH0013", "sp0.75_tf0.5_s0", "6"
SEEDS = range(5)

lab = pd.read_csv(RES / "mm_proto_rna_tf" / f"labels_{RUN}.csv.gz", index_col=0, dtype=str)[KEY]
rna = e0.load_slide(DATA, RUN)
e0.preprocess(rna)
tf = e0.load_tf_layer(DATA, RUN)
e0.preprocess_tf(tf)
assert list(rna.obs_names) == list(lab.index) == list(tf.obs_names)
band = (lab == BAND).to_numpy()
b_rna, b_tf = rna[band].copy(), tf[band].copy()
for a in (b_rna, b_tf):
    sc.pp.neighbors(a, n_neighbors=10, use_rep="X_pca", random_state=0)
    sq.gr.spatial_neighbors(a, coord_type="grid", n_neighs=6)
m = md.MuData({"rna": ad.AnnData(obs=pd.DataFrame(index=b_rna.obs_names), obsp={"connectivities": b_rna.obsp["connectivities"]}),
               "tf": ad.AnnData(obs=pd.DataFrame(index=b_rna.obs_names), obsp={"connectivities": b_tf.obsp["connectivities"]})})
m.obsp["spatial_connectivities"] = b_rna.obsp["spatial_connectivities"]


def cluster(res_e, res_s, seed):
    sl.spatialleiden_multimodal(m, resolution={"rna": res_e, "tf": res_e, "spatial": res_s},
                                layer_weights={"rna": 1.0, "tf": 0.5, "spatial": 0.75},
                                directed={"rna": False, "tf": False, "spatial": False},
                                n_iterations=-1, random_state=seed, key_added="sub")
    return m.obs["sub"].astype(str).to_numpy()


# geometry of the band
sp = pd.read_csv(DATA / "visium_tables" / f"{RUN}_spots.csv").set_index("barcode").loc[lab.index][band]
pos = {(r, c) for r, c in zip(sp.array_row, sp.array_col)}
steps = [(0, 2), (0, -2), (1, 1), (1, -1), (-1, 1), (-1, -1)]
depth = {}  # grid steps to the nearest spot outside the band (BFS from the rim)
frontier = [(r, c) for r, c in pos if any((r + dr, c + dc) not in pos for dr, dc in steps)]
for p in frontier:
    depth[p] = 1
d = 1
while frontier:
    nxt = []
    for r, c in frontier:
        for dr, dc in steps:
            q = (r + dr, c + dc)
            if q in pos and q not in depth:
                depth[q] = d + 1
                nxt.append(q)
    frontier, d = nxt, d + 1
sp["depth"] = [depth[(r, c)] for r, c in zip(sp.array_row, sp.array_col)]
xy = sp[["pxl_col_fullres", "pxl_row_fullres"]].to_numpy(float)
xy0 = xy - xy.mean(0)
axis = np.linalg.svd(xy0, full_matrices=False)[2][0]
sp["along"] = xy0 @ axis / np.ptp(xy0 @ axis)  # -0.5..0.5 along the band
print(f"band: {band.sum()} spots; depth (grid steps from edge) counts: {sp.depth.value_counts().sort_index().to_dict()}")

props = [c for c in sp.columns if c.startswith("prop_")]
counts = ad.AnnData(rna.layers["counts"][band], obs=pd.DataFrame(index=b_rna.obs_names), var=pd.DataFrame(index=rna.var_names))
sc.pp.normalize_total(counts, 1e4)
sc.pp.log1p(counts)
rows, labels = [], {}
for k in (2, 3):
    labs = []
    for s in SEEDS:
        r_e, _ = e0.search_resolution_steps(lambda r: len(set(cluster(r, 1.0, s))), k, start=0.3)
        lab_k = cluster(r_e, 1.0, s)
        if len(set(lab_k)) != k:
            r_s, _ = e0.search_resolution_steps(lambda r: len(set(cluster(r_e, r, s))), k, start=1.0)
            lab_k = cluster(r_e, r_s, s)
        labs.append(lab_k)
        labels[f"k{k}_s{s}"] = lab_k
    ok = [x for x in labs if len(set(x)) == k]
    stab = np.mean([adjusted_rand_score(a, b) for a, b in combinations(ok, 2)])
    l0 = labs[0]
    core = (sp.depth >= 2).astype(int).to_numpy()
    half = (sp.along > 0).astype(int).to_numpy()
    print(f"\nK={k}: {len(ok)}/5 seeds reached K; seed stability ARI {stab:.2f}; "
          f"ARI vs core/rim {adjusted_rand_score(core, l0):.2f}; ARI vs two halves along the band "
          f"{adjusted_rand_score(half, l0):.2f}")
    summ = sp.assign(sub=l0).groupby("sub").agg(n=("depth", "size"), mean_depth=("depth", "mean"),
                                               along_mean=("along", "mean"), along_sd=("along", "std"))
    summ = summ.join(sp.assign(sub=l0).groupby("sub")[props].mean())
    print(summ.round(2).to_string())
    counts.obs["sub"] = l0
    sc.tl.rank_genes_groups(counts, "sub", method="wilcoxon")
    for g in sorted(set(l0)):
        df = sc.get.rank_genes_groups_df(counts, g).query("pvals_adj < 0.01 and logfoldchanges > 1").head(12)
        print(f"  sub {g} top genes: {', '.join(df.names)}")
        rows.append({"K": k, "sub": g, "seed_stability_ARI": round(stab, 3), "n": int((l0 == g).sum()),
                     "top_genes": ";".join(df.names)})
pd.DataFrame(rows).to_csv(OUT / "vessel_band_subclusters.csv", index=False)
pd.DataFrame(labels, index=b_rna.obs_names).assign(depth=sp.depth.values, along=sp.along.values).to_csv(
    OUT / "vessel_band_sublabels.csv.gz")
