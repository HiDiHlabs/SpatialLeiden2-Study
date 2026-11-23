"""Does each modality / combination recover the P9 rep2 artery (vessel band)?

Reference: the vessel domain of RNA + TF + space (mm_proto_rna_tf sp0.75_tf0.5_s0, domain 6,
266 spots) and its artery-wall core (nested K = 2 sub-cluster, 55 spots; vessel_band_subclusters.py).
For every run (all K = 9 authors' cell-type niches, 5 seeds) take the domain that best overlaps
the band and report:
    jaccard_band    overlap of that domain with the band (1 = same spots)
    wall_in_best    share of the 55 wall spots that fall in one single domain (the artery found as one unit)
    wall_purity     share of that wall-holding domain's spots that are band spots (the domain is the vessel, not
                    a large region that happens to contain it)
Means over seeds that reached K.

Usage (from Derived_Modalities/): python 04_multimodal_weight_grid/proto_rna_tf/vessel_band_across_modalities.py results/kuppe_2022
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

R = Path(sys.argv[1])
RUN = "ACH0013"
ref = pd.read_csv(R / "mm_proto_rna_tf" / f"labels_{RUN}.csv.gz", index_col=0, dtype=str)["sp0.75_tf0.5_s0"]
band = set(ref.index[ref == "6"])
sub = pd.read_csv(R / "mm_proto_rna_tf" / "vessel_band" / "vessel_band_sublabels.csv.gz", index_col=0, dtype=str)
wall = set(sub.index[sub.k2_s0 == sub.k2_s0[sub.k2_s0.index.isin(band)].value_counts().index[-1]])
RUNS = [("RNA only", "exp0_k_celltype_niche", "r0_s{s}"), ("RNA + space", "exp0_k_celltype_niche", "r0.75_s{s}"),
        ("TF only", "exp0tf_k_celltype_niche", "r0_s{s}"), ("TF + space", "exp0tf_k_celltype_niche", "r0.75_s{s}"),
        ("CCC only", "exp0ccc_k_celltype_niche", "r0_s{s}"), ("CCC + space 0.25", "exp0ccc_k_celltype_niche", "r0.25_s{s}"),
        ("CCC + space 0.75", "exp0ccc_k_celltype_niche", "r0.75_s{s}"),
        ("Metabolic only", "exp0met_k_celltype_niche", "r0_s{s}"), ("Metabolic + space", "exp0met_k_celltype_niche", "r0.75_s{s}"),
        ("RNA + TF + space", "mm_configs", "rna_tf_s{s}"), ("RNA + CCC + space 0.25", "mm_configs_sp0.25", "rna_ccc_s{s}"),
        ("RNA + CCC + space 0.75", "mm_configs", "rna_ccc_s{s}"), ("RNA + metabolic + space", "mm_configs", "rna_met_s{s}"),
        ("RNA + TF + CCC + metabolic + space", "mm_configs", "rna_all_s{s}")]
rows = []
for name, folder, key in RUNS:
    t = pd.read_csv(R / folder / f"labels_{RUN}.csv.gz", index_col=0, dtype=str)
    for s in range(5):
        lab = t[key.format(s=s)]
        if lab.nunique() != 9:
            continue
        doms = {d: set(lab.index[lab == d]) for d in lab.unique()}
        jac = max(len(v & band) / len(v | band) for v in doms.values())
        wd = lab[list(wall)].value_counts()
        best = doms[wd.index[0]]
        rows.append({"run": name, "seed": s, "jaccard_band": jac, "wall_in_best": wd.iloc[0] / len(wall),
                     "wall_purity": len(best & band) / len(best), "wall_domain_size": len(best)})
d = pd.DataFrame(rows)
out = d.groupby("run", sort=False).agg(seeds=("seed", "size"), jaccard_band=("jaccard_band", "mean"),
                                       wall_in_best=("wall_in_best", "mean"), wall_purity=("wall_purity", "mean"),
                                       wall_domain_size=("wall_domain_size", "median")).round(2)
(R / "mm_proto_rna_tf" / "vessel_band").mkdir(exist_ok=True)
d.to_csv(R / "mm_proto_rna_tf" / "vessel_band" / "vessel_band_across_modalities.csv", index=False)
print(f"band {len(band)} spots, wall {len(wall)} spots")
print(out.to_string())
