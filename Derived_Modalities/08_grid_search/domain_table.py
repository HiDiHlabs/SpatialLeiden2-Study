"""Domain table for the grid search (strategy section 4.2 of
handoff/HANDOFF_2026-11-25_joint_multimodal_findings_and_strategy.md): one row per domain of every run in
results/kuppe_2022/grid_search/runs/, with the filter criteria, then the candidate list.

Per domain (slide = unit for shares; zones from the authors' major_labl):
  n_spots, n_slides_2pct, n_patients_2pct   slides / patients where it is >= 2 % of the slide
  coherence_gap                             share of its spatial grid edges that stay inside it, minus the
                                            same share expected by chance (its share of the slide's spots)
  share_<zone>                              mean share over the slides of that zone
  best_zone, fold_vs_ctrl                   patient zone with the highest mean share; ratio to control
                                            (control share floored at 0.005)
  p_perm                                    permutation test of best-zone mean share > control mean share
                                            (labels of the slides of those two zones permuted, 2,000 times)
  seed_repro                                number of the other seeds (same config, K, weight) with a domain
                                            at Jaccard >= 0.5
  j_rna                                     best Jaccard with any RNA-only domain at the same target K and seed
                                            (any spatial weight)
  j_ct_niche, j_mol_niche                   best Jaccard with the authors' cell-type / molecular niches
Candidate (multimodal, disease-specific) = seed_repro >= 3 of 4 (or all available) AND coherence_gap > 0.2
  AND n_slides_2pct >= 3 AND n_patients_2pct >= 3 AND fold_vs_ctrl >= 3 AND j_rna < 0.5
  AND j_ct_niche < 0.5 AND j_mol_niche < 0.5. "Disease domain" = the same without the three j < 0.5 rules.

Usage (from Derived_Modalities/): python 08_grid_search/domain_table.py results/kuppe_2022/grid_search
Output: domain_table.csv, candidates.csv (multimodal + disease), disease_domains.csv
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path(sys.argv[1])
ZONES = ["CTRL", "RZ", "BZ", "IZ", "FZ"]
rng = np.random.default_rng(0)
spots = pd.read_csv(OUT / "spots.csv.gz", index_col=0)
runs = pd.read_csv(OUT / "runs.csv").drop_duplicates("run_id", keep="last")
slide = pd.factorize(spots["sample"])[0]
slides = pd.unique(spots["sample"])
zone_of_slide = spots.groupby("sample")["major_labl"].first().reindex(slides).to_numpy()
patient_of_slide = spots.groupby("sample")["patient"].first().reindex(slides).to_numpy()
n_per_slide = np.bincount(slide)

# spatial grid edges (same construction as the runs: grid neighbours within a slide)
rows, cols = spots["array_row"].to_numpy() if "array_row" in spots else None, None
if rows is None:
    import squidpy as sq
    import anndata as ad
    g = ad.AnnData(obs=pd.DataFrame({"sample": pd.Categorical(spots["sample"])}, index=spots.index))
    g.obsm["spatial"] = spots[["pxl_col_fullres", "pxl_row_fullres"]].to_numpy(float)
    sq.gr.spatial_neighbors(g, library_key="sample", coord_type="grid", n_neighs=6)
    EI, EJ = g.obsp["spatial_connectivities"].nonzero()


def load(rid):
    return np.load(OUT / "runs" / f"{rid}.npz")["labels"].astype(np.int64)


def best_jaccard(lab, other):
    """For each domain of lab: best Jaccard with any domain of other."""
    a, b = pd.factorize(lab)[0], pd.factorize(other)[0]
    ct = np.zeros((a.max() + 1, b.max() + 1))
    np.add.at(ct, (a, b), 1)
    na, nb = ct.sum(1), ct.sum(0)
    jac = ct / (na[:, None] + nb[None, :] - ct)
    keys = pd.factorize(lab)[1]
    return dict(zip(keys, jac.max(1)))


ct_lab = spots["celltype_niche"].to_numpy()
mol_lab = spots["molecular_niche"].to_numpy()
cache = {}
records = []
for _, r in runs.iterrows():
    lab = load(r.run_id)
    rna_ids = runs[(runs.config == "rna") & (runs.target_k == r.target_k) & (runs.seed == r.seed)].run_id
    j_rna = {}
    if r.config != "rna":
        for rid in rna_ids:
            for d, j in best_jaccard(lab, load(rid)).items():
                j_rna[d] = max(j_rna.get(d, 0), j)
    sib = runs[(runs.config == r.config) & (runs.target_k == r.target_k) & (runs.spatial_weight == r.spatial_weight)
               & (runs.seed != r.seed)].run_id
    repro = {}
    for rid in sib:
        for d, j in best_jaccard(lab, load(rid)).items():
            repro[d] = repro.get(d, 0) + int(j >= 0.5)
    j_ct, j_mol = best_jaccard(lab, ct_lab), best_jaccard(lab, mol_lab)
    same = lab[EI] == lab[EJ]
    for d in np.unique(lab):
        m = lab == d
        per_slide = np.bincount(slide[m], minlength=len(slides)) / n_per_slide
        zs = {z: per_slide[zone_of_slide == z].mean() for z in ZONES}
        best_zone = max(ZONES[1:], key=lambda z: zs[z])
        fold = zs[best_zone] / max(zs["CTRL"], 0.005)
        sel = np.isin(zone_of_slide, ["CTRL", best_zone])
        x, grp = per_slide[sel], zone_of_slide[sel] == best_zone
        obs = x[grp].mean() - x[~grp].mean()
        perm = [(lambda g_: x[g_].mean() - x[~g_].mean())(rng.permutation(grp)) for _ in range(2000)]
        p_perm = (1 + sum(v >= obs for v in perm)) / 2001
        e_in = m[EI]
        coh = same[e_in].mean() if e_in.any() else np.nan
        exp_ = np.average(per_slide[slide[m]]) if m.any() else np.nan  # chance: neighbour in same domain
        big = per_slide >= 0.02
        records.append({"run_id": r.run_id, "config": r.config, "target_k": r.target_k, "spatial_weight": r.spatial_weight,
                        "seed": r.seed, "n_domains_run": r.n_domains, "domain": int(d), "n_spots": int(m.sum()),
                        "n_slides_2pct": int(big.sum()), "n_patients_2pct": int(len(set(patient_of_slide[big]))),
                        "coherence_gap": coh - exp_, **{f"share_{z}": zs[z] for z in ZONES}, "best_zone": best_zone,
                        "fold_vs_ctrl": fold, "p_perm": p_perm, "seed_repro": repro.get(d, 0), "n_other_seeds": len(sib),
                        "j_rna": j_rna.get(d, 1.0 if r.config == "rna" else 0.0), "j_ct_niche": j_ct[d], "j_mol_niche": j_mol[d]})
dt = pd.DataFrame(records)
dt.to_csv(OUT / "domain_table.csv", index=False)
need = np.minimum(3, dt.n_other_seeds)
disease = (dt.seed_repro >= need) & (dt.coherence_gap > 0.2) & (dt.n_slides_2pct >= 3) & (dt.n_patients_2pct >= 3) \
    & (dt.fold_vs_ctrl >= 3)
cand = disease & (dt.config != "rna") & (dt.j_rna < 0.5) & (dt.j_ct_niche < 0.5) & (dt.j_mol_niche < 0.5)
dt[disease].sort_values("fold_vs_ctrl", ascending=False).to_csv(OUT / "disease_domains.csv", index=False)
dt[cand].sort_values("fold_vs_ctrl", ascending=False).to_csv(OUT / "candidates.csv", index=False)
print(f"{len(runs)} runs, {len(dt)} domains; disease-enriched {int(disease.sum())}; multimodal-only candidates {int(cand.sum())}")
show = ["config", "target_k", "spatial_weight", "seed", "domain", "n_spots", "best_zone", "fold_vs_ctrl", "p_perm",
        "coherence_gap", "seed_repro", "j_rna", "j_ct_niche", "j_mol_niche"]
print(dt[cand].sort_values("fold_vs_ctrl", ascending=False)[show].head(30).round(3).to_string(index=False))
