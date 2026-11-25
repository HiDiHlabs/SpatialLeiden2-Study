"""UMAPs of the three Harmony-integrated layers (RNA, TF, CCC; all 27 slides, 88,704 spots) coloured by
region (zone), patient, the authors' cell-type and molecular niches, the joint SpatialLeiden domains of several combinations, and two activity read-outs;
a-c UMAPs per embedding; d stacked bars: mean share of each domain / niche / patient per region
(over the slides of that region). Plus a test of which clustering separates which activity best (table only, eta2_activity.csv) and the
ARI between the clusterings (ari_between_clusterings.csv). Figure in the style of Kuppe ED Fig. 3a-b:
no frames, UMAP arrows, each legend under its column.

Embeddings (no recomputation of the integration):
  RNA  results/kuppe_2022/exp1_harmony_covariates/embeddings.npz, key "authors"
  TF   results/kuppe_2022/tf_harmony/tf_embeddings.npz, key "X_harmony" (its stored UMAP is reused)
  CCC  results/kuppe_2022/ccc_harmony_ours/ccc_embeddings.npz, key "X_harmony" (our LIANA+ layer)
UMAP: 15-NN on the embedding (the clustering graph), scanpy UMAP, seed 0; cached in
results/kuppe_2022/joint_umaps/umaps.npz.
Domains (seed 0): RNA + space 0.75, TF + space 0.75, CCC + space 0.25, RNA + TF + space 0.75,
RNA + CCC + space 0.25. Colours: D1.. by size; each clustering's domains take the colour of the
RNA + space domain they mostly come from; a second domain from the same source gets its own colour.
A thin dashed ring (RNA map only) marks the UMAP area where the stressed-CM domain makes up > 40 % of the spots
(smoothed). Pink and the label "stressed CM" are reserved for the verified stressed-cardiomyocyte domain
(NPPB / ANKRD1 high): RNA + TF + space D5 and RNA + CCC + space D6.
Activity read-outs per spot:
  TF activity     number of TFs with ULM score > 0 and padj < 0.05 (622 TFs scored on all 27 slides)
  communication   mean local CCC score over the 286 ligand-receptor pairs used for the CCC layer
Test: eta^2 = share of the spot-level variance of each read-out explained by the domains (one-way
ANOVA R^2), per clustering. Higher = the clustering separates that activity more.

Patient and region colours: the authors' palettes (Kuppe et al. 2022 ED Fig. 3a-b), imported from
common/plot_patient_region_umap.py.

Version 2 (third argument "v2"): the activity columns are replaced by the clusterings without the
spatial graph (RNA only, TF only, CCC only: spatial weight 0, K = 9) -> fig_joint_umaps_v2.

Version 3 ("v3"): no niche columns; RNA / TF / CCC alone (no space), each with space, and the combinations
-> fig_joint_umaps_v3.

Usage (from Derived_Modalities/, env spatialleiden_joint):
    python 07_integration/fig_joint_umaps.py data/kuppe_2022 results/kuppe_2022 [v2|v3]
Output: results/kuppe_2022/joint_umaps/fig_joint_umaps.png/.pdf, eta2_activity.csv
"""

import sys
from pathlib import Path

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from plot_patient_region_umap import PATIENT_COLOURS, REGION_COLOURS  # noqa: E402  authors' palettes (ED Fig. 3a-b)

DATA, RES = Path(sys.argv[1]), Path(sys.argv[2])
OUT = RES / "joint_umaps"
OUT.mkdir(parents=True, exist_ok=True)
MM = 1 / 25.4
FS, FH, FL = 5, 6, 8  # Nature style: body text 5 pt, headings 6 pt bold, panel letters 8 pt bold (Arial)
mpl.rcParams.update({"font.family": "Arial", "font.size": FS, "axes.titlesize": FH, "pdf.fonttype": 42})
ZONES = ["CTRL", "RZ", "BZ", "IZ", "FZ"]
ZONE_COL = {"CTRL": REGION_COLOURS["control"], **{z: REGION_COLOURS[z] for z in ["RZ", "BZ", "IZ", "FZ"]}}
RNA_COLS = ["#1f4e9c", "#8fb8e0", "#c49c6b", "#d62728", "#ff9f40", "#8a5fbf", "#7f7f7f", "#d4c62a", "#cfcfcf"]
NEW_COLS = ["#e43fb1", "#17becf", "#2ca02c", "#9c6b30", "#4b0082", "#00b894", "#ff6f91", "#6b8e23"]
EMB = {"RNA": ("exp1_harmony_covariates/embeddings.npz", "authors"), "TF": ("tf_harmony/tf_embeddings.npz", "X_harmony"),
       "CCC": ("ccc_harmony_ours/ccc_embeddings.npz", "X_harmony")}
CLUST = {"RNA + space": ("exp2_joint_harmony", "r0.75_s0"), "TF + space": ("exp2_joint_harmony_tf", "r0.75_s0"),
         "CCC + space": ("exp2_joint_harmony_ccc_ours", "r0.25_s0"), "RNA + TF + space": ("exp3_joint_rna_tf", "r0.75_s0"),
         "RNA + CCC + space": ("exp3_joint_rna_ccc_ours", "r0.25_s0"),
         "RNA only": ("exp2_joint_harmony", "r0_s0"), "TF only": ("exp2_joint_harmony_tf", "r0_s0"),
         "CCC only": ("exp2_joint_harmony_ccc_ours", "r0_s0")}
VER = sys.argv[3] if len(sys.argv) > 3 else "v1"
V2 = VER in ("v2", "v3")  # v2: no-space clusterings instead of the activity columns
V3 = VER == "v3"  # v3: no niche columns; every layer alone, with space, and the combinations
ALL4 = (RES / "exp3b_joint_rna_tf_ccc" / "exp2_labels.csv.gz").exists()  # RNA + TF + CCC + space (exp3b), once run
if ALL4:
    CLUST["RNA + TF + CCC + space"] = ("exp3b_joint_rna_tf_ccc", "r0.75_s0")
SHOW = (["RNA only", "TF only", "CCC only"] if V2 else []) + ["RNA + space"] + \
    (["TF + space", "CCC + space"] if V3 else []) + ["RNA + TF + space", "RNA + CCC + space"] + \
    (["RNA + TF + CCC + space"] if ALL4 else [])

names = None
X = {}
for k, (f, key) in EMB.items():
    e = np.load(RES / f, allow_pickle=True)
    idx = pd.Index(e["obs_names"].astype(str))
    names = idx if names is None else names
    X[k] = pd.DataFrame(e[key], index=idx).loc[names].to_numpy(np.float32)
cache = OUT / "umaps.npz"
if cache.exists():
    U = {k: v for k, v in np.load(cache, allow_pickle=True).items() if k != "obs_names"}
else:
    U = {}
    for k, x in X.items():
        a = ad.AnnData(obs=pd.DataFrame(index=names))
        a.obsm["X_e"] = x
        sc.pp.neighbors(a, n_neighbors=15, use_rep="X_e", random_state=0)
        sc.tl.umap(a, random_state=0)
        U[k] = a.obsm["X_umap"]
        print(f"UMAP {k} done", flush=True)
    np.savez_compressed(cache, obs_names=np.asarray(names), **U)

run = names.str.split("_", n=1).str[0]
spots = []
for r in pd.unique(run):
    t = pd.read_csv(DATA / "visium_tables" / f"{r}_spots.csv", usecols=["barcode", "major_labl", "patient", "celltype_niche", "molecular_niche"])
    t.index = r + "_" + t.barcode
    spots.append(t)
spots = pd.concat(spots).loc[names]

# activity read-outs
pairs = (RES / "ccc_harmony_ours" / "ccc_pairs.txt").read_text().split()
comm, tfact = pd.Series(np.nan, index=names), pd.Series(np.nan, index=names)
shared_tf = None
tfs = {}
for r in pd.unique(run):
    c = pd.read_csv(DATA / "modalities" / "ccc_all27_ours" / f"ccc_{r}.csv", index_col=0)
    comm.loc[r + "_" + c.index] = c[pairs].mean(axis=1).to_numpy()
    t = ad.read_h5ad(DATA / "modalities" / "collectri_all27" / f"collectri_ulm_{r}.h5ad")
    tfs[r] = t
    shared_tf = set(t.var_names) if shared_tf is None else shared_tf & set(t.var_names)
shared_tf = sorted(shared_tf)
for r, t in tfs.items():
    t = t[:, shared_tf]
    act = ((np.asarray(t.X) > 0) & (np.asarray(t.layers["padj"]) < 0.05)).sum(1)
    tfact.loc[r + "_" + t.obs_names] = act
print(f"{len(pairs)} CCC pairs, {len(shared_tf)} TFs on all slides", flush=True)


def load(folder, col):
    lab = pd.read_csv(RES / folder / "exp2_labels.csv.gz", index_col=0, dtype=str)[col].loc[names]
    order = lab.value_counts().index
    return lab.map({d: f"D{i + 1}" for i, d in enumerate(order)})


labs = {k: load(*v) for k, v in CLUST.items()}
rna = labs["RNA + space"]
rna_doms = [f"D{i}" for i in range(1, rna.nunique() + 1)]
rna_col = dict(zip(rna_doms, RNA_COLS))
# stressed-cardiomyocyte domain (NPPB / ANKRD1 high; checked per slide, see the 24 Nov handoff)
STRESSED = {"RNA + TF + space": "D5", "RNA + CCC + space": "D6"}
if ALL4:  # its stressed-CM domain = the one overlapping RNA + TF D5 best, if Jaccard > 0.5
    ref = labs["RNA + TF + space"] == "D5"
    lab4 = labs["RNA + TF + CCC + space"]
    jac = {d: ((lab4 == d) & ref).sum() / ((lab4 == d) | ref).sum() for d in lab4.unique()}
    best = max(jac, key=jac.get)
    print(f"RNA + TF + CCC + space: best match to the stressed-CM domain {best}, Jaccard {jac[best]:.2f}")
    if jac[best] > 0.5:
        STRESSED["RNA + TF + CCC + space"] = best
cols, new_by_src = {"RNA + space": rna_col}, {}
for k in ["RNA + TF + space", "RNA + CCC + space"] + (["RNA + TF + CCC + space"] if ALL4 else []) + \
        ["TF + space", "CCC + space", "RNA only", "TF only", "CCC only"]:  # combinations first: pink = new
    lab = labs[k]
    doms = [f"D{i}" for i in range(1, lab.nunique() + 1)]
    src = pd.crosstab(lab, rna, normalize="index").reindex(index=doms, columns=rna_doms, fill_value=0)
    col, used = {}, set()
    for d in doms:
        m = src.loc[d].idxmax()
        if m not in used:
            col[d], _ = rna_col[m], used.add(m)
        elif STRESSED.get(k) == d:  # the verified stressed-cardiomyocyte domain: pink
            col[d] = NEW_COLS[0]
        else:  # any other new domain: its own colour (never pink)
            new_by_src[(k, d)] = NEW_COLS[1 + len(new_by_src) % (len(NEW_COLS) - 1)]
            col[d] = new_by_src[(k, d)]
    cols[k] = col


def eta2(y, g):
    y = pd.Series(np.asarray(y, float)); g = pd.Series(np.asarray(g))
    ok = y.notna().to_numpy()
    y, g = y[ok], g[ok]
    return float(((y.groupby(g).transform("mean") - y.mean()) ** 2).sum() / ((y - y.mean()) ** 2).sum())


rows = [{"clustering": k, "eta2_TF_activity": eta2(tfact, lab), "eta2_communication": eta2(comm, lab),
         "eta2_zone_dummy": np.nan} for k, lab in labs.items()]
et = pd.DataFrame(rows)
et.to_csv(OUT / "eta2_activity.csv", index=False)
print(et.round(3).to_string(index=False))

# agreement between the clusterings (each has its own 9 domains)
from itertools import combinations  # noqa: E402

from sklearn.metrics import adjusted_rand_score  # noqa: E402

ari = pd.DataFrame(np.eye(len(labs)), index=list(labs), columns=list(labs))
for a_, b_ in combinations(labs, 2):
    ari.loc[a_, b_] = ari.loc[b_, a_] = adjusted_rand_score(labs[a_], labs[b_])
ari.round(3).to_csv(OUT / "ari_between_clusterings.csv")
print(ari.round(2).to_string())

# domain descriptions: two cell types with the highest median cell2location proportion
props_tab = []
for r in pd.unique(run):
    t = pd.read_csv(DATA / "visium_tables" / f"{r}_spots.csv").set_index("barcode")
    t.index = r + "_" + t.index
    props_tab.append(t[[c for c in t.columns if c.startswith("prop_")]])
props_tab = pd.concat(props_tab).loc[names]


ABBR = {"Myeloid": "Myel", "vSMCs": "vSMC", "prolif": "Prolif", "Lymphoid": "Lymph", "Neuronal": "Neur",
        "Endo": "EC", "Adipo": "Adip"}


def describe(lab):
    med = props_tab.groupby(lab.to_numpy()).median()
    return {d: "+".join(ABBR.get(c.replace("prop_", ""), c.replace("prop_", ""))
                        for c in med.loc[d].sort_values(ascending=False).index[:2]) for d in med.index}


# ---- figure: Kuppe ED Fig. 3a-b style (no frames, UMAP arrows, legends under each column) ----
def ring(ax, u, mask, level=0.4):
    """Light outline of the UMAP area where the stressed-CM domain dominates: smoothed local share of
    its spots (200 x 200 grid, Gaussian sigma 3 bins) at `level`, only where spots are dense enough."""
    from scipy.ndimage import gaussian_filter
    xe = np.linspace(u[:, 0].min(), u[:, 0].max(), 201)
    ye = np.linspace(u[:, 1].min(), u[:, 1].max(), 201)
    h_all, _, _ = np.histogram2d(u[:, 0], u[:, 1], bins=[xe, ye])
    h_in, _, _ = np.histogram2d(u[mask, 0], u[mask, 1], bins=[xe, ye])
    h_all, h_in = gaussian_filter(h_all, 3), gaussian_filter(h_in, 3)
    frac = np.where(h_all > 0.5, h_in / np.maximum(h_all, 1e-9), 0)
    xc, yc = (xe[:-1] + xe[1:]) / 2, (ye[:-1] + ye[1:]) / 2
    ax.contour(xc, yc, frac.T, levels=[level], colors="black", linewidths=0.35, linestyles="--", alpha=0.7)


def umap_arrows(ax, draw_arrows):
    for sp_ in ax.spines.values():
        sp_.set_visible(False)
    ax.set(xticks=[], yticks=[])
    if not draw_arrows:  # arrows only on the first UMAP of each row
        return
    arrow = dict(arrowstyle="-|>", lw=0.6, color="black", mutation_scale=5, shrinkA=0, shrinkB=0)
    ax.annotate("", xy=(0.3, 0.0), xytext=(0.0, 0.0), xycoords="axes fraction", arrowprops=arrow)
    ax.annotate("", xy=(0.0, 0.3), xytext=(0.0, 0.0), xycoords="axes fraction", arrowprops=arrow)
    # labels start at the corner and run along their arrow (outside the corner), so they never meet
    ax.text(0.0, -0.025, "UMAP 1", transform=ax.transAxes, ha="left", va="top", fontsize=FS)
    ax.text(-0.025, 0.0, "UMAP 2", transform=ax.transAxes, ha="left", va="bottom", rotation=90,
            rotation_mode="anchor", fontsize=FS)


pat = sorted(spots.patient.unique(), key=lambda p: int(p[1:]))
pat_col = {p: PATIENT_COLOURS[p] for p in pat}
# authors' niche colours: ArchR "stallion" in niche order (as ed_fig3abc.py for cell-type niches)
STALLION = list(PATIENT_COLOURS.values())
CT = sorted(spots.celltype_niche.unique(), key=lambda x: int(x.split("_")[-1]))
MOLN = sorted(spots.molecular_niche.unique(), key=lambda x: int(x.split("_")[-1]))
ct_col = dict(zip(CT, STALLION))
mol_col = dict(zip(MOLN, STALLION))
colnames = ["Region", "Patient"] + ([] if V3 else ["Cell-type\nniche", "Molecular\nniche"]) + SHOW + \
    ([] if V2 else ["TF activity", "Communication"])
fig = plt.figure(figsize=(183 * MM, 160 * MM))
gs = fig.add_gridspec(3, 1, height_ratios=[3, 0.75, 1.1], hspace=0.12)
top = gs[0].subgridspec(3, len(colnames), wspace=0.12, hspace=0.18)
bars = gs[1].subgridspec(1, len(colnames), wspace=0.12)
leg = gs[2].subgridspec(1, len(colnames), wspace=0.0)
rng = np.random.default_rng(0)
order = rng.permutation(len(names))
vlim = {"t": np.nanquantile(tfact, [0.01, 0.99]), "c": np.nanquantile(comm, [0.01, 0.99])}
for i, (emb, u) in enumerate([(k, U[k]) for k in EMB]):
    u = u[order]
    for j, cn in enumerate(colnames):
        a = fig.add_subplot(top[i, j])
        if cn == "Region":
            c = spots.major_labl.map(ZONE_COL).to_numpy()[order]
        elif cn == "Patient":
            c = spots.patient.map(pat_col).to_numpy()[order]
        elif cn.startswith("Cell-type"):
            c = spots.celltype_niche.map(ct_col).to_numpy()[order]
        elif cn.startswith("Molecular"):
            c = spots.molecular_niche.map(mol_col).to_numpy()[order]
        elif cn in labs:
            c = labs[cn].map(cols[cn]).to_numpy()[order]
        elif cn == "TF activity":
            c = mpl.colormaps["viridis"](mpl.colors.Normalize(*vlim["t"], clip=True)(tfact.to_numpy()[order]))
        else:
            c = mpl.colormaps["magma"](mpl.colors.Normalize(*vlim["c"], clip=True)(comm.to_numpy()[order]))
        a.scatter(u[:, 0], u[:, 1], c=c, s=0.12, linewidths=0, rasterized=True)
        umap_arrows(a, j == 0)
        if cn in STRESSED and emb == "RNA":  # ring on the RNA map only (the domain forms a region there)
            ring(a, u, (labs[cn].to_numpy()[order] == STRESSED[cn]))
        if i == 0:
            a.set_title(cn.replace(" + space", "\n+ space").replace("RNA + TF + CCC", "RNA+TF+CCC")
                        .replace("RNA + CCC\n", "RNA+CCC\n").replace("RNA + TF\n", "RNA+TF\n") if cn in labs else
                        {"TF activity": "TF activity\n(# active)", "Communication": "Communi-\ncation (CCC)"}.get(cn, cn),
                        fontweight="bold", fontsize=FH)
        if j == 0:
            a.text(-0.22, 0.5, f"{emb} embedding", transform=a.transAxes, rotation=90, ha="right",
                   va="center", fontweight="bold", fontsize=FH)
            a.text(-0.42, 1.08, "abc"[i], transform=a.transAxes, ha="left", va="bottom", fontweight="bold", fontsize=FL)


# d: stacked bars, mean share of each domain / niche / patient per region (slide = unit)
def stacked(ax, lab, cmap, order, first):
    df = pd.DataFrame({"run": run, "zone": spots.major_labl.to_numpy(), "lab": np.asarray(lab)})
    share = pd.crosstab([df.run, df.zone], df.lab, normalize="index").reindex(columns=order, fill_value=0)
    mean = share.groupby(level="zone").mean().reindex(ZONES)
    bottom = np.zeros(len(ZONES))
    for k in order:
        ax.bar(range(len(ZONES)), mean[k].to_numpy(), bottom=bottom, color=cmap[k], width=0.8, linewidth=0)
        bottom += mean[k].to_numpy()
    ax.set_xticks(range(len(ZONES)), ["ctrl" if z == "CTRL" else z for z in ZONES], fontsize=FS, rotation=90)
    ax.set_ylim(0, 1)
    ax.set_yticks([0, 0.5, 1] if first else [])
    ax.tick_params(labelsize=FS, length=1.5, width=0.4, pad=1)
    for sp_ in ("top", "right"):
        ax.spines[sp_].set_visible(False)
    for sp_ in ("left", "bottom"):
        ax.spines[sp_].set_linewidth(0.4)
    if first:
        ax.set_ylabel("mean share\nper region", fontsize=FS)


first_bar = True
for j, cn in enumerate(colnames):
    a = fig.add_subplot(bars[0, j])
    if cn == "Patient":
        stacked(a, spots.patient.to_numpy(), pat_col, pat, first_bar)
    elif cn.startswith("Cell-type"):
        stacked(a, spots.celltype_niche.to_numpy(), ct_col, CT, first_bar)
    elif cn.startswith("Molecular"):
        stacked(a, spots.molecular_niche.to_numpy(), mol_col, MOLN, first_bar)
    elif cn in labs:
        doms = [f"D{k}" for k in range(1, labs[cn].nunique() + 1)]
        stacked(a, labs[cn].to_numpy(), cols[cn], doms, first_bar)
    else:
        a.axis("off")
        if cn == "Region":
            a.text(-0.42, 0.95, "d", transform=a.transAxes, ha="left", va="bottom", fontweight="bold", fontsize=FL)
        continue
    first_bar = False


def dot_legend(ax, items, title, ncol=1, fs=FS, x0=0.0):
    h = [plt.Line2D([], [], ls="", marker="o", markersize=2.4, color=c, label=l) for l, c in items]
    lg = ax.legend(handles=h, title=title, loc="upper left", bbox_to_anchor=(x0, 1.0), frameon=False, fontsize=fs,
                   title_fontsize=FH, ncol=ncol, handletextpad=0.0, handlelength=0.8, labelspacing=0.3,
                   columnspacing=0.3, borderaxespad=0, borderpad=0)
    lg._legend_box.align = "left"
    for t, (_l, c) in zip(lg.get_texts(), items):
        if "stressed CM" in _l:
            t.set_fontweight("bold")
    ax.axis("off")


for j, cn in enumerate(colnames):
    a = fig.add_subplot(leg[0, j])
    if cn == "Region":
        dot_legend(a, [("control" if z == "CTRL" else z, ZONE_COL[z]) for z in ZONES], "Region")
    elif cn == "Patient":
        dot_legend(a, [(p, pat_col[p]) for p in pat], "Patient", ncol=2)
    elif cn.startswith("Cell-type"):
        dot_legend(a, [(n.split("_")[-1], ct_col[n]) for n in CT], "Cell-type\nniche")
    elif cn.startswith("Molecular"):
        dot_legend(a, [(n.split("_")[-1], mol_col[n]) for n in MOLN], "Molecular\nniche", ncol=2)
    elif cn in labs:
        desc = describe(labs[cn])
        doms = [f"D{k}" for k in range(1, labs[cn].nunique() + 1)]
        dot_legend(a, [(f"{d} stressed CM" if STRESSED.get(cn) == d else f"{d} {desc[d]}", cols[cn][d])
                       for d in doms], "Domains",
                   x0=0.1 if colnames[j - 1] in STRESSED else 0.0)  # room after a bold "stressed CM" label
    else:
        key, cm = ("t", "viridis") if cn == "TF activity" else ("c", "magma")
        a.axis("off")
        cax = a.inset_axes([0.15, 0.82, 0.7, 0.06])
        cb = fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(*vlim[key]), cmap=cm), cax=cax,
                          orientation="horizontal")
        cb.ax.tick_params(labelsize=FS, length=1.5, width=0.4)
        cb.outline.set_linewidth(0.4)
        cb.set_label("# TFs, padj < 0.05" if key == "t" else "mean local score", fontsize=FS)
tag = {"v1": "fig_joint_umaps", "v2": "fig_joint_umaps_v2", "v3": "fig_joint_umaps_v3"}[VER]
fig.savefig(OUT / f"{tag}.png", dpi=250, bbox_inches="tight")
fig.savefig(OUT / f"{tag}.pdf", dpi=300, bbox_inches="tight")
