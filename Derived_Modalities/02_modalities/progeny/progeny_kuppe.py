"""PROGENy pathway activity per Visium spot, computed the way Kuppe et al. 2022 did.

Authors' code (saezlab/visium_heart, st_snRNAseq/utils/funcomics.R, add_path_activities,
called with top = 500 and assay = "RNA"):
    model  <- progeny::getModel(organism = "Human", top = 500)
    scores <- t(model)[, common_genes] %*% GetAssayData(slide, assay = "RNA")[common_genes, ]
    scores <- scale(t(scores))      # z-score each pathway across the slide's spots
Here the same weighted sum is applied to the log-normalised expression shipped in
data/kuppe_2022/processed1/*.h5ad (X, LogNormalization = Seurat RNA data slot). The model is
the PROGENy human top-500 table from decoupler (dc.op.progeny), saved next to the output.

Usage (from Derived_Modalities/, decoupler env):
    python 02_modalities/progeny/progeny_kuppe.py data/kuppe_2022
Writes data/kuppe_2022/modalities/progeny/progeny_<hca_sample_id>.csv (spots x 14 pathways).
"""

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

DATA = Path(sys.argv[1])
OUT = DATA / "modalities" / "progeny"
OUT.mkdir(parents=True, exist_ok=True)
MODEL = OUT / "progeny_human_top500.csv"
if not MODEL.exists():
    import decoupler as dc

    dc.op.progeny(organism="human", top=500).to_csv(MODEL, index=False)

# processed1 file -> Visium run (checked: identical barcode sets with visium_tables/*_spots.csv)
SLIDES = {
    "Visium_control_P7": "ACH003",
    "Visium_RZ_BZ_P3": "10X0026",
    "Visium_IZ_P3": "10X0017",
    "Visium_FZ_P14": "ACH005",
}

net = pd.read_csv(MODEL)
W = net.pivot_table(index="target", columns="source", values="weight", fill_value=0.0)
for fname, sid in SLIDES.items():
    a = ad.read_h5ad(DATA / "processed1" / f"{fname}.h5ad")
    spots = pd.read_csv(DATA / "visium_tables" / f"{sid}_spots.csv", usecols=["barcode"])
    assert set(a.obs_names) == set(spots.barcode), sid
    genes = W.index.intersection(a.var_names)
    X = a[:, genes].X
    S = np.asarray(X @ W.loc[genes].to_numpy())
    S = (S - S.mean(0)) / S.std(0, ddof=1)  # R scale()
    df = pd.DataFrame(S, index=a.obs_names, columns=W.columns)
    df.index.name = "barcode"
    df.to_csv(OUT / f"progeny_{sid}.csv")
    print(sid, df.shape, f"{len(genes)} model genes on slide")
