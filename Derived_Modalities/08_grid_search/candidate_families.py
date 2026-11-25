"""Collapse the grid-search candidates (candidates.csv: domain instances, one per run) into distinct regions
("families"): two candidate domains belong to the same family if their spot sets overlap with Jaccard >= 0.5
(single linkage = connected components). Per family: how many instances, configurations, K and spatial
weights produce it, the spots found in >= 50 % of its instances (consensus region), the best zone and fold
over control, top cell types, and overlap with the stressed-cardiomyocyte domain of 24 Nov
(exp3_joint_rna_tf r0.75_s0, D5 by size).

Usage (from Derived_Modalities/): python 08_grid_search/candidate_families.py data/kuppe_2022 results/kuppe_2022
Output: results/kuppe_2022/grid_search/candidate_families.csv, family_consensus_spots.csv.gz
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components

DATA, RES = Path(sys.argv[1]), Path(sys.argv[2])
G = RES / "grid_search"
spots = pd.read_csv(G / "spots.csv.gz", index_col=0)
cand = pd.read_csv(G / "candidates.csv").reset_index(drop=True)
n = len(spots)
rows, cols = [], []
for i, r in cand.iterrows():
    lab = np.load(G / "runs" / f"{r.run_id}.npz")["labels"]
    idx = np.flatnonzero(lab == r.domain)
    rows += [i] * len(idx)
    cols += idx.tolist()
M = sp.csr_matrix((np.ones(len(rows), np.float32), (rows, cols)), shape=(len(cand), n))
inter = (M @ M.T).toarray()
size = np.asarray(M.sum(1)).ravel()
jac = inter / (size[:, None] + size[None, :] - inter)
ncomp, fam = connected_components(sp.csr_matrix(jac >= 0.5), directed=False)
cand["family"] = fam

ref = pd.read_csv(RES / "exp3_joint_rna_tf" / "exp2_labels.csv.gz", index_col=0, dtype=str)["r0.75_s0"]
order = ref.value_counts().index
ref = (ref.map({d: f"D{i + 1}" for i, d in enumerate(order)}) == "D5").loc[spots.index].to_numpy()
tab = []
for s in spots["sample"].unique():
    t = pd.read_csv(DATA / "visium_tables" / f"{s}_spots.csv").set_index("barcode")
    t.index = s + "_" + t.index
    tab.append(t[[c for c in t.columns if c.startswith("prop_")]])
props = pd.concat(tab).loc[spots.index]
rows_out, cons_cols = [], {}
for f, g in cand.groupby("family"):
    cnt = np.asarray(M[g.index.to_numpy()].sum(0)).ravel()
    consensus = cnt >= 0.5 * len(g)
    if consensus.sum() == 0:
        consensus = cnt > 0
    comp = props[consensus].median().sort_values(ascending=False)
    zone = g.best_zone.mode().iloc[0]
    rows_out.append({
        "family": f, "n_instances": len(g), "configs": ";".join(sorted(g.config.unique())),
        "n_configs": g.config.nunique(), "target_k": ";".join(map(str, sorted(g.target_k.unique()))),
        "spatial_weights": ";".join(map(str, sorted(g.spatial_weight.unique()))), "seeds": g.seed.nunique(),
        "consensus_spots": int(consensus.sum()), "best_zone": zone, "median_fold_vs_ctrl": g.fold_vs_ctrl.median(),
        "median_p_perm": g.p_perm.median(), "median_coherence_gap": g.coherence_gap.median(),
        "median_j_rna": g.j_rna.median(), "max_j_rna": g.j_rna.max(),
        "top_celltypes": "+".join(c.replace("prop_", "") for c in comp.index[:3]),
        "jaccard_stressed_CM": (consensus & ref).sum() / (consensus | ref).sum(),
        "zones_of_spots": ", ".join(f"{z} {v:.0%}" for z, v in
                                   spots.loc[consensus, "major_labl"].value_counts(normalize=True).head(3).items())})
    cons_cols[f"family_{f}"] = consensus
fams = pd.DataFrame(rows_out).sort_values(["n_configs", "n_instances"], ascending=False)
fams.to_csv(G / "candidate_families.csv", index=False)
pd.DataFrame(cons_cols, index=spots.index).to_csv(G / "family_consensus_spots.csv.gz")
pd.set_option("display.width", 250)
print(f"{len(cand)} candidate domains -> {ncomp} families")
print(fams.drop(columns=["configs", "target_k", "spatial_weights"]).round(2).head(20).to_string(index=False))
