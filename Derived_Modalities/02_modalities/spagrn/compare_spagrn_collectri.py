"""SpaGRN regulons vs CollecTRI TF activity on P1 control, P9 IZ and P9 IZ rep2.

The two answer different questions:
    CollecTRI + ULM (decoupler)  prior-knowledge regulons (curated TF -> target, signed), the same
                                 for every slide; per-spot score = t-value of a linear fit of
                                 expression on the target weights
    SpaGRN                       regulons inferred per slide (spatially variable genes, TF-target
                                 co-expression + cisTarget motif pruning, activating regulons
                                 "TF(+)"); per-spot score = AUCell
Comparison, per slide:
    1. how many regulons, how many TFs shared with CollecTRI (TFs with >= 5 targets)
    2. target overlap for shared TFs (Jaccard of SpaGRN targets vs CollecTRI targets), against a
       null of mismatched TF pairs
    3. spot-wise agreement for shared TFs: Spearman correlation of AUCell vs ULM score, against a
       null of mismatched TF pairs (same slide)
    4. regulons found on all three slides
    5. the P9 rep2 artery band (vessel domain of RNA + TF + space): top regulons / TFs by
       standardised difference band vs rest, for both methods
Figure: counts + correlation distributions + vessel ranking + maps of selected TFs (both methods;
scales 1st-99th percentile per map: SpaGRN AUCell viridis, CollecTRI ULM magma).

Usage (env with decoupler >= 2.1 and scanpy, e.g. masld_atlas; from Derived_Modalities/):
    python 02_modalities/spagrn/compare_spagrn_collectri.py data/kuppe_2022 results/kuppe_2022
Output: results/kuppe_2022/spagrn/compare_collectri/
"""

import sys
from pathlib import Path

import anndata as ad
import decoupler as dc
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import EllipseCollection
from scipy.spatial import cKDTree
from scipy.stats import rankdata

DATA, RES = Path(sys.argv[1]), Path(sys.argv[2])
OUT = RES / "spagrn" / "compare_collectri"
OUT.mkdir(parents=True, exist_ok=True)
MM = 1 / 25.4
mpl.rcParams.update({"font.family": "Arial", "font.size": 6, "axes.titlesize": 6.5, "pdf.fonttype": 42})
SLIDES = {"10X001": "P1 control", "ACH0012": "P9 IZ", "ACH0013": "P9 IZ rep2"}
SLIDE_COL = {"10X001": "#0072B2", "ACH0012": "#E69F00", "ACH0013": "#D55E00"}
CANDIDATES = ["HIF1A", "RELA", "NFKB1", "MEF2C", "MYOCD", "SRF", "TCF21", "STAT3", "GATA4", "TBX20", "SPI1",
              "CEBPB", "EGR1", "FOSL2", "KLF4", "TEAD1", "NFE2L2", "IRF8", "RUNX1", "MEF2A"]
rng = np.random.default_rng(0)

net = dc.op.collectri(organism="human")
ctri = {tf: set(g.target) for tf, g in net.groupby("source")}


def colcorr(A, B):
    """Spearman correlation per column pair (spots x features)."""
    ra = np.apply_along_axis(rankdata, 0, A)
    rb = np.apply_along_axis(rankdata, 0, B)
    ra = (ra - ra.mean(0)) / ra.std(0)
    rb = (rb - rb.mean(0)) / rb.std(0)
    return (ra * rb).mean(0)


S, rows, corr_rows = {}, [], []
for run in SLIDES:
    g = ad.read_h5ad(RES / "spagrn" / run / f"{run}_spagrn.h5ad")
    auc = g.obsm["auc_mtx"].copy()
    auc.columns = [c.replace("(+)", "").replace("(-)", "") for c in auc.columns]
    regs = {k.replace("(+)", "").replace("(-)", ""): set(v) for k, v in g.uns["regulon_dict"].items()}
    t = ad.read_h5ad(DATA / "modalities" / f"collectri_ulm_{run}.h5ad")
    ulm = pd.DataFrame(np.asarray(t.X), index=t.obs_names, columns=t.var_names)
    common_spots = ulm.index.intersection(auc.index)
    auc, ulm = auc.loc[common_spots], ulm.loc[common_spots]
    keep = auc.columns[auc.std() > 0]
    shared = sorted(set(keep) & set(ulm.columns))
    r = colcorr(auc[shared].to_numpy(), ulm[shared].to_numpy())
    perm = rng.permutation(len(shared))
    r_null = colcorr(auc[shared].to_numpy(), ulm[np.array(shared)[perm]].to_numpy())
    jac = [len(regs[tf] & ctri.get(tf, set())) / len(regs[tf] | ctri.get(tf, set())) for tf in shared]
    jac_null = [len(regs[a] & ctri.get(b, set())) / len(regs[a] | ctri.get(b, set()))
                for a, b in zip(shared, np.array(shared)[perm])]
    S[run] = dict(auc=auc, ulm=ulm, regs=regs, shared=shared, r=pd.Series(r, index=shared), r_null=r_null,
                  xy=t.obs.loc[common_spots, ["pxl_col_fullres", "pxl_row_fullres"]].to_numpy(float))
    rows.append({"slide": run, "spagrn_regulons": len(regs), "collectri_tfs": ulm.shape[1], "shared_tfs": len(shared),
                 "median_regulon_size": int(np.median([len(v) for v in regs.values()])),
                 "spearman_median": float(np.median(r)), "spearman_null_median": float(np.median(r_null)),
                 "share_r_gt_0.3": float(np.mean(r > 0.3)), "jaccard_median": float(np.median(jac)),
                 "jaccard_null_median": float(np.median(jac_null)),
                 "share_any_target_overlap": float(np.mean([j > 0 for j in jac]))})
    for tf, rv, jv in zip(shared, r, jac):
        corr_rows.append({"slide": run, "tf": tf, "spearman": rv, "jaccard_targets": jv,
                          "n_targets_spagrn": len(regs[tf]), "n_targets_collectri": len(ctri.get(tf, ()))})
summary = pd.DataFrame(rows)
summary.to_csv(OUT / "summary.csv", index=False)
per_tf = pd.DataFrame(corr_rows)
per_tf.to_csv(OUT / "per_tf_agreement.csv", index=False)
all3 = sorted(set.intersection(*(set(S[r]["regs"]) for r in SLIDES)))
pd.Series(all3, name="regulon_on_all_three_slides").to_csv(OUT / "regulons_on_all_slides.csv", index=False)
print(summary.round(3).to_string(index=False))
print(f"regulons on all three slides: {len(all3)}; of them in CollecTRI: {len(set(all3) & set(ctri))}")
wide = per_tf.pivot(index="tf", columns="slide", values="spearman")
print("best-agreeing TFs (mean over slides where present):")
print(wide.mean(1).sort_values(ascending=False).head(12).round(2).to_string())

# vessel band on P9 rep2
lab = pd.read_csv(RES / "mm_proto_rna_tf" / f"labels_ACH0013.csv.gz", index_col=0, dtype=str)["sp0.75_tf0.5_s0"]
s = S["ACH0013"]
band = (lab.loc[s["auc"].index] == "6").to_numpy()


def zdiff(M):
    return ((M[band].mean() - M[~band].mean()) / M.std().replace(0, np.nan)).sort_values(ascending=False)


vz_sp, vz_ct = zdiff(s["auc"]), zdiff(s["ulm"])
pd.DataFrame({"spagrn_regulon": vz_sp.index[:20], "spagrn_z": vz_sp.values[:20].round(2),
              "collectri_tf": vz_ct.index[:20], "collectri_z": vz_ct.values[:20].round(2)}).to_csv(
    OUT / "vessel_band_top_tfs.csv", index=False)
print("vessel band, top SpaGRN regulons:", ", ".join(f"{k} {v:.2f}" for k, v in vz_sp.head(10).items()))
print("vessel band, top CollecTRI TFs:  ", ", ".join(f"{k} {v:.2f}" for k, v in vz_ct.head(10).items()))
for tf in ["MYOCD", "SRF", "TCF21"]:
    print(f"  {tf}: SpaGRN regulon on rep2 {'yes' if tf in s['regs'] else 'no'}"
          + (f", band z {vz_sp[tf]:.2f}, rank {list(vz_sp.index).index(tf) + 1}/{len(vz_sp)}" if tf in vz_sp else ""))

# ---- figure ----
in_any = [tf for tf in CANDIDATES if sum(tf in S[r]["shared"] for r in SLIDES) >= 2]
show = in_any[:4]
print("TFs shown in the maps:", show)


def draw(ax, xy, colours):
    pitch = np.median(cKDTree(xy).query(xy, k=2)[0][:, 1])
    d = np.full(len(xy), 1.05 * pitch)
    coll = EllipseCollection(d, d, np.zeros(len(xy)), units="xy", offsets=xy, offset_transform=ax.transData,
                             facecolors=colours, linewidths=0)
    coll.set_rasterized(True)
    ax.add_collection(coll)
    ax.set_xlim(xy[:, 0].min() - pitch, xy[:, 0].max() + pitch)
    ax.set_ylim(xy[:, 1].max() + pitch, xy[:, 1].min() - pitch)
    ax.set_aspect("equal")
    ax.axis("off")


fig = plt.figure(figsize=(183 * MM, 165 * MM))
gs = fig.add_gridspec(2, 1, height_ratios=[0.6, 1.6], hspace=0.3)
top = gs[0].subgridspec(1, 4, wspace=0.7, width_ratios=[0.8, 1.2, 0.8, 0.8])
ax = fig.add_subplot(top[0])
x = np.arange(len(SLIDES))
ax.bar(x - 0.2, summary.spagrn_regulons, 0.4, color="#7B1FA2", label="SpaGRN regulons")
ax.bar(x + 0.2, summary.shared_tfs, 0.4, color="#9E9E9E", label="shared with\nCollecTRI")
ax.set_xticks(x, [SLIDES[r] for r in SLIDES], rotation=30, ha="right")
ax.set_ylabel("number of TFs")
ax.set_ylim(0, 480)
ax.legend(frameon=False, fontsize=5, loc="upper left", ncol=1)
ax.set_title("a  Regulons", fontweight="bold", loc="left")
ax2 = fig.add_subplot(top[1])
for i, run in enumerate(SLIDES):
    for vals, pos, colr in [(S[run]["r"].to_numpy(), i - 0.18, "#7B1FA2"), (S[run]["r_null"], i + 0.18, "#9E9E9E")]:
        vp = ax2.violinplot([vals], positions=[pos], widths=0.32, showmedians=True, showextrema=False)
        for b in vp["bodies"]:
            b.set_facecolor(colr); b.set_alpha(0.7)
        vp["cmedians"].set_color("black")
ax2.axhline(0, color="black", lw=0.4)
ax2.set_xticks(range(3), [SLIDES[r] for r in SLIDES], rotation=30, ha="right")
ax2.set_ylabel("Spearman, AUCell vs ULM")
ax2.legend(handles=[mpl.patches.Patch(color="#7B1FA2", label="same TF"), mpl.patches.Patch(color="#9E9E9E", label="mismatched TF (null)")],
           frameon=False, fontsize=5, loc="upper left")
ax2.set_title("b  Spot-wise agreement", fontweight="bold", loc="left")
k = 10
axes_c = []
for col_i, (vz, name, colr) in enumerate([(vz_sp, "SpaGRN", "#7B1FA2"), (vz_ct, "CollecTRI", "#9E9E9E")]):
    a3 = fig.add_subplot(top[2 + col_i])
    a3.barh(np.arange(k)[::-1], vz.values[:k], color=colr)
    a3.set_yticks(np.arange(k)[::-1], vz.index[:k], fontsize=5.5)
    a3.set_xlabel("z, artery vs rest")
    a3.set_title(("c  P9 rep2 artery\n" if col_i == 0 else "\n") + name, fontweight="bold", loc="left")
    axes_c.append(a3)
for a_ in [ax, ax2] + axes_c:
    for sp_ in ("top", "right"):
        a_.spines[sp_].set_visible(False)
    a_.tick_params(labelsize=5.5, length=1.5, width=0.4)
bot = gs[1].subgridspec(len(show), 6, wspace=0.05, hspace=0.12)
for i, tf in enumerate(show):
    for j, run in enumerate(SLIDES):
        for m, (method, M, cm) in enumerate([("SpaGRN", S[run]["auc"], "viridis"), ("CollecTRI", S[run]["ulm"], "magma")]):
            axm = fig.add_subplot(bot[i, 2 * j + m])
            if tf in M.columns and M[tf].std() > 0:
                v = M[tf].to_numpy()
                nrm = mpl.colors.Normalize(np.quantile(v, 0.01), np.quantile(v, 0.99), clip=True)
                draw(axm, S[run]["xy"], mpl.colormaps[cm](nrm(v)))
                if method == "SpaGRN" and tf in S[run]["r"]:
                    axm.text(0.02, 0.02, f"ρ={S[run]['r'][tf]:.2f}", transform=axm.transAxes, fontsize=5)
            else:
                draw(axm, S[run]["xy"], [(0.9, 0.9, 0.9, 1)] * len(S[run]["xy"]))
                axm.text(0.5, 0.5, "no regulon" if method == "SpaGRN" else "not scored", transform=axm.transAxes,
                         ha="center", va="center", fontsize=5, color="#555555")
            if i == 0:
                axm.set_title(f"{SLIDES[run]}\n{method}", fontsize=6, fontweight="bold", color=SLIDE_COL[run])
            if j == 0 and m == 0:
                axm.text(-0.08, 0.5, tf, transform=axm.transAxes, rotation=90, ha="right", va="center",
                         fontweight="bold", fontsize=7)
fig.text(0.02, 0.6, "d", fontweight="bold", fontsize=8)
fig.savefig(OUT / "fig_spagrn_vs_collectri.png", dpi=250, bbox_inches="tight")
fig.savefig(OUT / "fig_spagrn_vs_collectri.pdf", dpi=300, bbox_inches="tight")
