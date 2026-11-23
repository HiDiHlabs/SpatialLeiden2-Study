"""Synergistic Information Score (SIS; Richter et al. 2026, Nat. Mach. Intell., s42256-026-01328-3)
for our derived modalities: does TF activity, CCC or scCellFie metabolic tasks add task-relevant
information beyond RNA, or is it redundant?

SIS(Y; Z3) = [P(Y; Z3) - max(P(Y; Z1), P(Y; Z2))] / max(P(Y; Z1), P(Y; Z2))
with P = linear-probe performance (paper eq. 13, MMI redundancy): logistic regression -> F1 macro
for classification, ridge -> R2 for regression. Positive = fusion makes information accessible
beyond the strongest single layer; ~0 = redundant; < 0 = fusion hurts.

Representations (per slide, built exactly as for SpatialLeiden, exp0_sweep.py):
    Z1 = RNA PCA 30; Z2 = TF / CCC / metabolic PCA (<= 30); Z3 = concatenation [Z1, Z2]
    (the paper's "Concat" aggregation baseline; SpatialLeiden fuses graphs, not embeddings,
    so it has no representation to probe).
Null control: Z3 with Z2 rows permuted (same dimension, no paired information) -> SIS_null,
the gain expected from extra probe capacity alone.
Targets Y:
    celltype_niche   authors' cell-type niche (F1 macro)
    molecular_niche  authors' molecular niche (F1 macro; defined from RNA, so biased to RNA)
    composition      cell2location proportions of 11 cell types (R2, multi-output)
    neighbour_comp   mean composition of the spots 2-3 grid rings away (R2), the paper's
                     spatial-context task in Visium form
Validation: 5-fold spatial-block CV within each slide (k-means on coordinates into 5 blocks),
so test spots are not neighbours of training spots (the CCC layer is spatially smoothed).

Usage (env spatialleiden, from Derived_Modalities/):
    python 05_synergy_sis/sis_layers.py data/kuppe_2022 results/kuppe_2022/sis
"""

import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
spec = importlib.util.spec_from_file_location(
    "exp0_sweep", HERE.parents[1] / "03_single_modality_clustering" / "exp0_layer_ratio_sweep" / "exp0_sweep.py")
e0 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(e0)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from sklearn.cluster import KMeans  # noqa: E402
from sklearn.linear_model import LogisticRegression, Ridge  # noqa: E402
from sklearn.metrics import f1_score, r2_score  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402

DATA, OUT = Path(sys.argv[1]), Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)
SLIDES = ["10X001", "ACH0012", "ACH0013"]
STEPS = [(0, 2), (0, -2), (1, 1), (1, -1), (-1, 1), (-1, -1)]


def embed(run):
    rna = e0.load_slide(DATA, run)
    e0.preprocess(rna)
    z = {"RNA": rna.obsm["X_pca"]}
    for name, ld, pp in [("TF", e0.load_tf_layer, e0.preprocess_tf), ("CCC", e0.load_ccc_layer, e0.preprocess_ccc),
                         ("Metabolic", e0.load_sccellfie_layer, e0.preprocess_ccc)]:
        a = ld(DATA, run)
        pp(a)
        assert list(a.obs_names) == list(rna.obs_names)
        z[name] = a.obsm["X_pca"]
    return rna, z


def neighbour_comp(sp, props):
    pos = {(r, c): n for n, (r, c) in enumerate(zip(sp.array_row, sp.array_col))}
    P = sp[props].to_numpy()
    out = np.full_like(P, np.nan)
    for (r, c), n in pos.items():
        ring = [pos[(r + dr * k, c + dc * k)] for dr, dc in STEPS for k in (2, 3) if (r + dr * k, c + dc * k) in pos]
        if ring:
            out[n] = P[ring].mean(0)
    return out


def probe(X, y, folds, kind):
    pred = np.zeros_like(y, dtype=float) if kind == "reg" else np.empty(len(y), dtype=object)
    for f in np.unique(folds):
        tr, te = folds != f, folds == f
        sc_ = StandardScaler().fit(X[tr])
        if kind == "clf":
            m = LogisticRegression(max_iter=3000, C=1.0).fit(sc_.transform(X[tr]), y[tr])
        else:
            m = Ridge(alpha=1.0).fit(sc_.transform(X[tr]), y[tr])
        pred[te] = m.predict(sc_.transform(X[te]))
    if kind == "clf":
        return f1_score(y, pred, average="macro")
    return r2_score(y, pred, multioutput="uniform_average")


rows = []
rng = np.random.default_rng(0)
for run in SLIDES:
    rna, z = embed(run)
    sp = pd.read_csv(DATA / "visium_tables" / f"{run}_spots.csv").set_index("barcode").loc[rna.obs_names]
    props = [c for c in sp.columns if c.startswith("prop_")]
    folds = KMeans(5, n_init=10, random_state=0).fit_predict(sp[["pxl_col_fullres", "pxl_row_fullres"]].to_numpy(float))
    nb = neighbour_comp(sp, props)
    ok = ~np.isnan(nb).any(1)
    tasks = {"celltype_niche": (sp.celltype_niche.to_numpy(), "clf", np.ones(len(sp), bool)),
             "molecular_niche": (sp.molecular_niche.to_numpy(), "clf", np.ones(len(sp), bool)),
             "composition": (sp[props].to_numpy(), "reg", np.ones(len(sp), bool)),
             "neighbour_comp": (nb, "reg", ok)}
    for task, (y, kind, keep) in tasks.items():
        if kind == "clf":  # drop classes too rare to appear in every training split
            vc = pd.Series(y[keep]).value_counts()
            keep = keep & np.isin(y, vc.index[vc >= 10])
        p = {n: probe(Z[keep], y[keep], folds[keep], kind) for n, Z in z.items()}
        for other in ["TF", "CCC", "Metabolic"]:
            z3 = np.hstack([z["RNA"], z[other]])
            z3n = np.hstack([z["RNA"], z[other][rng.permutation(len(z[other]))]])
            p3, p3n = probe(z3[keep], y[keep], folds[keep], kind), probe(z3n[keep], y[keep], folds[keep], kind)
            best = max(p["RNA"], p[other])
            rows.append({"slide": run, "task": task, "metric": "F1_macro" if kind == "clf" else "R2",
                         "pair": f"RNA + {other}", "P_RNA": p["RNA"], "P_other": p[other], "P_concat": p3,
                         "P_concat_null": p3n, "SIS": (p3 - best) / best, "SIS_null": (p3n - best) / best})
        print(run, task, {k: round(v, 3) for k, v in p.items()}, flush=True)
d = pd.DataFrame(rows)
d.to_csv(OUT / "sis_layers.csv", index=False)
pd.set_option("display.width", 200)
print(d.round(3).to_string(index=False))
