#!/usr/bin/env python
"""Visualise LIANA+ spatial (local, bivariate) ligand-receptor scores across all samples.

Input: spot x ligand-receptor scores from all samples, in any one of these forms
  * one concatenated table (.h5ad, .csv, .tsv or .parquet) with a column naming the sample
    (auto-detected from: sample, sample_id, library_id, orig.ident, hca_sample_id, batch,
    or set it with --sample-col). Pair columns / var_names look like "LIGAND^RECEPTOR".
  * a folder of per-sample files (ccc_<sample>.csv, <sample>.csv, <sample>.h5ad ...);
    the sample name is the file name without the "ccc_" prefix.
A pair that a sample's run filtered out should be NaN (pd.concat does this). If missing pairs
were filled with 0 they cannot be told apart from real zeros.

Optional --spots-dir: the Kuppe Visium spot tables (<sample>_spots.csv with barcode,
patient, major_labl, pxl_row_fullres, pxl_col_fullres). They give each sample its patient
and zone (CTRL / RZ / BZ / IZ / FZ) and the coordinates for the spatial maps. Without it,
use --meta (CSV with columns sample, group[, patient]) and spatial maps come from
obsm["spatial"] of an .h5ad, if present.

Outputs (in --out):
  01_heatmap_top_pairs.png/pdf   sample-mean score of the most variable pairs, z-scored per pair
  02_dotplot_by_group.png/pdf    per group: mean score (colour) and share of spots > --thr (size)
  03_pair_coverage.png/pdf       how many pairs each sample has; in how many samples each pair is
  04_boxplots_selected.png/pdf   sample means of the selected pairs by group
  05_spatial_<pair>.png          maps of one pair on every sample, one shared colour scale
  sample_pair_summary.csv        long table: sample, group, patient, pair, mean, frac_gt_thr
  pair_group_test.csv            Kruskal-Wallis of sample means across groups, per pair (BH-adjusted)

Example:
  python plot_liana_all_samples.py liana_all_samples.h5ad --spots-dir data/kuppe_2022/visium_tables \
      --out liana_figs --pairs "SPP1^CD44" "FN1^ITGAV_ITGB1"

Needs: numpy, pandas, matplotlib, scipy; anndata for .h5ad; pyarrow for .parquet.
"""

import argparse
import re
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import EllipseCollection
from scipy import stats
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial import cKDTree

SAMPLE_COLS = ["sample", "sample_id", "library_id", "orig.ident", "hca_sample_id", "batch", "slide"]
GROUP_ORDER = ["CTRL", "RZ", "BZ", "IZ", "FZ"]
GROUP_COL = {"CTRL": "#0072B2", "RZ": "#009E73", "BZ": "#56B4E9", "IZ": "#E69F00", "FZ": "#D55E00"}
DEFAULT_PAIRS = ["SPP1^CD44", "FN1^ITGAV_ITGB1", "THBS1^CD36", "APP^CD74", "C3^LRP1", "VEGFB^NRP1"]

mpl.rcParams.update({"font.size": 7, "axes.titlesize": 8, "pdf.fonttype": 42,
                     "axes.spines.top": False, "axes.spines.right": False})


# ----------------------------------------------------------------------------- loading
def read_table(path):
    """Return (DataFrame of spots x columns, spatial coords array or None)."""
    suf = path.suffix.lower()
    if suf == ".h5ad":
        import anndata as ad
        a = ad.read_h5ad(path)
        X = a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X)
        df = pd.DataFrame(X.astype(np.float32), index=a.obs_names, columns=a.var_names)
        df = pd.concat([a.obs, df], axis=1)
        xy = np.asarray(a.obsm["spatial"], float)[:, :2] if "spatial" in a.obsm else None
        return df, xy
    if suf == ".parquet":
        return pd.read_parquet(path), None
    sep = "\t" if suf in (".tsv", ".txt") else ","
    return pd.read_csv(path, sep=sep, index_col=0), None


def load_scores(inp, sample_col, glob="*"):
    """Return (scores: spots x pairs float, sample: Series per spot, xy: DataFrame x/y or None)."""
    inp = Path(inp)
    if inp.is_dir():
        files = sorted(f for f in inp.glob(glob) if f.suffix.lower() in (".csv", ".tsv", ".parquet", ".h5ad"))
        if not files:
            raise SystemExit(f"No .csv/.tsv/.parquet/.h5ad files in {inp}")
        parts, xys = [], []
        for f in files:
            df, xy = read_table(f)
            if not any("^" in str(c) for c in df.columns):
                print(f"  skipped {f.name} (no LIGAND^RECEPTOR columns)")
                continue
            name = re.sub(r"^ccc_", "", f.stem)
            df.index = df.index.astype(str)
            df.insert(0, "__sample__", name)
            parts.append(df)
            xys.append(None if xy is None else pd.DataFrame(xy, index=df.index, columns=["x", "y"]))
            print(f"  {name}: {len(df)} spots")
        if not parts:
            raise SystemExit(f"No file in {inp} matching '{glob}' has LIGAND^RECEPTOR columns")
        df = pd.concat(parts, axis=0)  # pairs missing in a sample become NaN
        xy = pd.concat(xys) if all(x is not None for x in xys) else None
        sample_col = "__sample__"
    else:
        df, xy = read_table(inp)
        df.index = df.index.astype(str)
        xy = None if xy is None else pd.DataFrame(xy, index=df.index, columns=["x", "y"])
        if sample_col is None:
            sample_col = next((c for c in SAMPLE_COLS if c in df.columns), None)
        if sample_col is None or sample_col not in df.columns:
            other = [c for c in df.columns if "^" not in str(c)]
            raise SystemExit(f"No sample column found. Pass --sample-col; non-pair columns are: {other}")
    pairs = [c for c in df.columns if "^" in str(c)]
    if not pairs:
        raise SystemExit("No ligand-receptor columns (expected names like 'LIGAND^RECEPTOR').")
    sample = df[sample_col].astype(str)
    scores = df[pairs].apply(pd.to_numeric, errors="coerce").astype(np.float32)
    # make the index unique per spot (barcodes repeat across samples) but keep the barcode
    barcode = pd.Series(df.index, index=df.index)
    key = sample + "|" + barcode
    scores.index, sample.index = key.values, key.values
    if xy is not None:
        xy.index = key.values
    return scores, sample, barcode.values, xy


def strip_sample(bc, s):
    for pat in (f"{s}_", f"{s}-", f"{s}:"):
        if bc.startswith(pat):
            return bc[len(pat):]
    for pat in (f"_{s}", f"-{s}"):
        if bc.endswith(pat):
            return bc[: -len(pat)]
    return bc


def load_meta(samples, scores_index, barcodes, spots_dir, meta_csv, xy):
    """Per-sample group/patient, and per-spot coords (from spot tables if not in the input)."""
    meta = pd.DataFrame(index=samples, columns=["group", "patient"], dtype=object)
    spot_sample = pd.Series([k.split("|", 1)[0] for k in scores_index], index=scores_index)
    bc = pd.Series(barcodes, index=scores_index)
    coords = []
    if spots_dir is not None:
        for s in samples:
            f = Path(spots_dir) / f"{s}_spots.csv"
            if not f.exists():
                print(f"  ! no spot table for {s} ({f.name}); group/coords unknown")
                continue
            sp = pd.read_csv(f).set_index("barcode")
            if "major_labl" in sp:
                meta.loc[s, "group"] = sp["major_labl"].iloc[0]
            if "patient" in sp:
                meta.loc[s, "patient"] = sp["patient"].iloc[0]
            if xy is None:
                keys = spot_sample.index[spot_sample.values == s]
                b = bc.loc[keys].map(lambda x: strip_sample(x, s))
                hit = b.isin(sp.index)
                if hit.mean() < 0.9:
                    print(f"  ! {s}: only {hit.mean():.0%} of barcodes match the spot table")
                c = sp.loc[b[hit].values, ["pxl_col_fullres", "pxl_row_fullres"]].to_numpy(float)
                coords.append(pd.DataFrame(c, index=keys[hit.values], columns=["x", "y"]))
    if meta_csv is not None:
        m = pd.read_csv(meta_csv).astype({"sample": str}).set_index("sample")
        for c in ("group", "patient"):
            if c in m:
                meta[c] = m[c].reindex(meta.index).combine_first(meta[c])
    meta["group"] = meta["group"].fillna("NA")
    meta["patient"] = meta["patient"].fillna("NA")
    if xy is None and coords:
        xy = pd.concat(coords)
    return meta, xy


# ----------------------------------------------------------------------------- helpers
def order_samples(meta):
    rank = {g: i for i, g in enumerate(GROUP_ORDER)}
    return sorted(meta.index, key=lambda s: (rank.get(meta.loc[s, "group"], 99), meta.loc[s, "group"],
                                             meta.loc[s, "patient"], s))


def group_colour(g):
    return GROUP_COL.get(g, "#999999")


def save(fig, out, name):
    fig.savefig(out / f"{name}.png", dpi=250, bbox_inches="tight")
    fig.savefig(out / f"{name}.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)


def pretty(p):
    lig, rec = p.split("^", 1)
    return f"{lig} → {rec.replace('_', '/')}"


# ----------------------------------------------------------------------------- figures
def heatmap(mean, meta, samples, pairs, out, zscore=True):
    M = mean.loc[samples, pairs].T
    if zscore:
        M = M.sub(M.mean(axis=1), axis=0).div(M.std(axis=1).replace(0, np.nan), axis=0)
    filled = M.fillna(0).to_numpy()
    if len(pairs) > 2:
        M = M.iloc[leaves_list(linkage(filled, "average", metric="correlation" if zscore else "euclidean"))]
    h = max(3, 0.13 * len(pairs) + 1.2)
    fig, (axc, ax) = plt.subplots(2, 1, figsize=(max(4, 0.22 * len(samples) + 2.2), h),
                                  gridspec_kw={"height_ratios": [0.25, len(pairs) * 0.13], "hspace": 0.02})
    groups = meta.loc[samples, "group"]
    axc.imshow([[mpl.colors.to_rgb(group_colour(g)) for g in groups]], aspect="auto")
    axc.set_xticks([]), axc.set_yticks([])
    for sp in axc.spines.values():
        sp.set_visible(False)
    cmap = mpl.colormaps["RdBu_r" if zscore else "viridis"].copy()
    cmap.set_bad("#dddddd")
    lim = np.nanquantile(np.abs(M.to_numpy()), 0.98) if zscore else None
    im = ax.imshow(np.ma.masked_invalid(M.to_numpy()), aspect="auto", cmap=cmap,
                   vmin=-lim if zscore else None, vmax=lim)
    ax.set_yticks(range(len(M)), [pretty(p) for p in M.index], fontsize=5.5)
    ax.set_xticks(range(len(samples)),
                  [f"{s} ({meta.loc[s, 'patient']})" for s in samples], rotation=90, fontsize=5.5)
    for t, g in zip(ax.get_xticklabels(), groups):
        t.set_color(group_colour(g))
    cb = fig.colorbar(im, ax=[axc, ax], shrink=0.4, pad=0.01)
    cb.set_label("z-score of sample mean (per pair)" if zscore else "sample mean score")
    handles = [mpl.patches.Patch(color=group_colour(g), label=g) for g in dict.fromkeys(groups)]
    axc.legend(handles=handles, ncol=len(handles), loc="lower left", bbox_to_anchor=(0, 1.05),
               frameon=False, fontsize=6)
    axc.set_title(f"Top {len(pairs)} most variable ligand-receptor pairs across samples "
                  "(grey = pair not in that sample's results)", loc="left", pad=14)
    save(fig, out, "01_heatmap_top_pairs")


def dotplot(scores, sample, meta, pairs, thr, out):
    grp = sample.map(meta["group"])
    groups = [g for g in GROUP_ORDER if g in set(grp)] + sorted(set(grp) - set(GROUP_ORDER))
    S = scores[pairs]
    mean = S.groupby(grp.values).mean().reindex(groups)
    frac = (S > thr).astype(float).where(S.notna()).groupby(grp.values).mean().reindex(groups)
    order = mean.T.idxmax(axis=1).map({g: i for i, g in enumerate(groups)}).sort_values(kind="stable").index
    fig, ax = plt.subplots(figsize=(1.0 + 0.45 * len(groups), max(3, 0.13 * len(pairs) + 1)))
    yy, xx = np.meshgrid(range(len(order)), range(len(groups)), indexing="ij")
    m = mean[order].T.to_numpy()
    f = frac[order].T.to_numpy()
    sc = ax.scatter(xx.ravel(), yy.ravel(), s=np.nan_to_num(f.ravel()) * 60 + 0.5, c=m.ravel(),
                    cmap="viridis", vmin=0, vmax=np.nanquantile(m, 0.98), linewidths=0)
    ax.set_xticks(range(len(groups)), groups)
    for t, g in zip(ax.get_xticklabels(), groups):
        t.set_color(group_colour(g))
    ax.set_yticks(range(len(order)), [pretty(p) for p in order], fontsize=5.5)
    ax.set_ylim(len(order) - 0.5, -0.5)
    ax.set_xlim(-0.6, len(groups) - 0.4)
    fig.colorbar(sc, ax=ax, shrink=0.3, pad=0.02, label="mean score (all spots in group)")
    for v in (0.25, 0.5, 1.0):
        ax.scatter([], [], s=v * 60 + 0.5, c="grey", label=f"{v:.0%}")
    ax.legend(title=f"spots > {thr}", frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1), fontsize=6)
    ax.set_title("Ligand-receptor scores by zone", loc="left")
    save(fig, out, "02_dotplot_by_group")


def coverage(scores, sample, meta, samples, out):
    present = scores.notna().groupby(sample.values).any().reindex(samples)
    n_pairs = present.sum(axis=1)
    n_samples = present.sum(axis=0)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(8, 2.6), gridspec_kw={"width_ratios": [2.2, 1], "wspace": 0.3})
    a1.bar(range(len(samples)), n_pairs, color=[group_colour(meta.loc[s, "group"]) for s in samples])
    a1.set_xticks(range(len(samples)), samples, rotation=90, fontsize=5.5)
    a1.set_ylabel("ligand-receptor pairs scored")
    a1.set_title("Pairs per sample", loc="left")
    a2.hist(n_samples, bins=np.arange(0.5, len(samples) + 1.5), color="#555555")
    a2.set_xlabel("number of samples the pair is scored in")
    a2.set_ylabel("pairs")
    a2.set_title(f"{len(n_samples)} pairs in total, {int((n_samples == len(samples)).sum())} in all samples",
                 loc="left")
    save(fig, out, "03_pair_coverage")
    return n_samples


def boxplots(mean, meta, pairs, out):
    groups = [g for g in GROUP_ORDER if g in set(meta["group"])] + sorted(set(meta["group"]) - set(GROUP_ORDER))
    ncol = min(4, len(pairs))
    nrow = int(np.ceil(len(pairs) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.0 * ncol, 2.0 * nrow), squeeze=False)
    rng = np.random.default_rng(0)
    for ax, p in zip(axes.flat, pairs):
        data = [mean.loc[meta.index[meta["group"] == g], p].dropna().to_numpy() for g in groups]
        ax.boxplot(data, widths=0.55, showfliers=False, medianprops={"color": "black"},
                   boxprops={"linewidth": 0.6}, whiskerprops={"linewidth": 0.6}, capprops={"linewidth": 0.6})
        for i, (g, d) in enumerate(zip(groups, data), 1):
            ax.scatter(i + rng.uniform(-0.15, 0.15, len(d)), d, s=9, color=group_colour(g), zorder=3,
                       linewidths=0)
        ax.set_xticks(range(1, len(groups) + 1), groups, fontsize=6)
        ax.set_title(pretty(p), fontsize=7)
        ax.set_ylabel("sample mean score", fontsize=6)
    for ax in axes.flat[len(pairs):]:
        ax.axis("off")
    fig.tight_layout()
    save(fig, out, "04_boxplots_selected")


def spatial(scores, sample, xy, meta, samples, pair, out):
    vals = scores[pair]
    vmax = float(np.nanquantile(vals, 0.99)) or 1.0
    norm = mpl.colors.Normalize(0, vmax, clip=True)
    ncol = min(7, len(samples))
    nrow = int(np.ceil(len(samples) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(1.6 * ncol, 1.75 * nrow), squeeze=False,
                             gridspec_kw={"wspace": 0.05, "hspace": 0.25})
    for ax, s in zip(axes.flat, samples):
        ax.axis("off")
        keys = xy.index.intersection(sample.index[sample.values == s])
        g = meta.loc[s, "group"]
        ax.set_title(f"{s}\n{meta.loc[s, 'patient']} · {g}", fontsize=6, color=group_colour(g))
        if len(keys) < 2:
            ax.text(0.5, 0.5, "no coordinates", ha="center", va="center", transform=ax.transAxes, fontsize=5)
            continue
        c = xy.loc[keys].to_numpy()
        v = vals.loc[keys].to_numpy()
        pitch = np.median(cKDTree(c).query(c, k=2)[0][:, 1])
        if np.isnan(v).all():
            col = [(0.88, 0.88, 0.88, 1)] * len(c)
            ax.text(0.5, 0.5, "pair not scored", ha="center", va="center", transform=ax.transAxes, fontsize=5)
        else:
            col = mpl.colormaps["viridis"](norm(np.nan_to_num(v)))
        d = np.full(len(c), 1.05 * pitch)
        coll = EllipseCollection(d, d, np.zeros(len(c)), units="xy", offsets=c, offset_transform=ax.transData,
                                 facecolors=col, linewidths=0)
        coll.set_rasterized(True)
        ax.add_collection(coll)
        ax.set_xlim(c[:, 0].min() - pitch, c[:, 0].max() + pitch)
        ax.set_ylim(c[:, 1].max() + pitch, c[:, 1].min() - pitch)
        ax.set_aspect("equal")
    for ax in axes.flat[len(samples):]:
        ax.axis("off")
    cax = fig.add_axes([0.92, 0.35, 0.01, 0.3])
    fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap="viridis"), cax=cax, label="local score (0 to p99)")
    fig.suptitle(pretty(pair), x=0.1, y=1.0 + 0.12 / nrow, ha="left", fontsize=9)
    fig.savefig(out / f"05_spatial_{re.sub(r'[^A-Za-z0-9]+', '_', pair)}.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help="concatenated table (.h5ad/.csv/.tsv/.parquet) or folder of per-sample files")
    ap.add_argument("--out", default="liana_figures", help="output folder")
    ap.add_argument("--glob", default="*", help="folder input only: which files to read, e.g. 'ccc_*.csv'")
    ap.add_argument("--sample-col", default=None, help="column naming the sample (auto-detected if omitted)")
    ap.add_argument("--spots-dir", default=None, help="folder with <sample>_spots.csv (Kuppe visium_tables)")
    ap.add_argument("--meta", default=None, help="CSV with columns sample, group[, patient]")
    ap.add_argument("--pairs", nargs="*", default=None, help="pairs for box plots and spatial maps")
    ap.add_argument("--top", type=int, default=40, help="number of most variable pairs in heatmap/dotplot")
    ap.add_argument("--min-samples", type=float, default=0.5,
                    help="a pair must be scored in at least this share of samples to enter the top list")
    ap.add_argument("--thr", type=float, default=0.1, help="score threshold for 'spot is active' in the dotplot")
    ap.add_argument("--no-spatial", action="store_true", help="skip the spatial maps")
    a = ap.parse_args()

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    print("Loading scores ...")
    scores, sample, barcodes, xy = load_scores(a.input, a.sample_col, a.glob)
    samples_all = sorted(sample.unique())
    print(f"{len(scores)} spots, {scores.shape[1]} pairs, {len(samples_all)} samples")
    meta, xy = load_meta(samples_all, scores.index, barcodes, a.spots_dir, a.meta, xy)
    samples = order_samples(meta)
    meta = meta.loc[samples]

    mean = scores.groupby(sample.values).mean().reindex(samples)
    frac = (scores > a.thr).astype(float).where(scores.notna()).groupby(sample.values).mean().reindex(samples)

    n_samples = coverage(scores, sample, meta, samples, out)
    eligible = n_samples.index[n_samples >= a.min_samples * len(samples)]
    top = mean[eligible].var().sort_values(ascending=False).index[: a.top].tolist()
    print(f"{len(eligible)} pairs scored in >= {a.min_samples:.0%} of samples; plotting top {len(top)} by variance")
    heatmap(mean, meta, samples, top, out)
    dotplot(scores, sample, meta, top, a.thr, out)

    # per-pair test across groups on sample means (samples are the replicates, not spots)
    groups = [g for g in meta["group"].unique() if g != "NA"]
    rows = []
    if len(groups) > 1:
        for p in mean.columns:
            data = [mean.loc[meta.index[meta["group"] == g], p].dropna() for g in groups]
            data = [d for d in data if len(d) > 0]
            if len(data) > 1 and sum(len(d) for d in data) > len(data):
                try:
                    h, pv = stats.kruskal(*data)
                except ValueError:  # all values identical
                    h, pv = np.nan, np.nan
                rows.append({"pair": p, "n_samples": int(sum(len(d) for d in data)), "H": h, "p": pv,
                             **{f"median_{g}": float(np.median(mean.loc[meta.index[meta["group"] == g], p].dropna()))
                                if mean.loc[meta.index[meta["group"] == g], p].notna().any() else np.nan
                                for g in groups}})
    test = pd.DataFrame(rows)
    if len(test):
        ok = test["p"].notna()
        test.loc[ok, "p_adj_BH"] = stats.false_discovery_control(test.loc[ok, "p"])
        test.sort_values("p").to_csv(out / "pair_group_test.csv", index=False)

    pairs = a.pairs or [p for p in DEFAULT_PAIRS if p in mean.columns]
    if not pairs and len(test):
        pairs = test.sort_values("p")["pair"].head(8).tolist()
    missing = [p for p in pairs if p not in mean.columns]
    if missing:
        print(f"  ! pairs not found in any sample, skipped: {missing}")
    pairs = [p for p in pairs if p in mean.columns]
    if pairs:
        boxplots(mean, meta, pairs, out)
        if not a.no_spatial:
            if xy is None:
                print("  ! no coordinates (give --spots-dir or an .h5ad with obsm['spatial']); no spatial maps")
            else:
                for p in pairs:
                    spatial(scores, sample, xy, meta, samples, p, out)

    long = (mean.stack(future_stack=True).rename("mean").to_frame()
            .join(frac.stack(future_stack=True).rename(f"frac_gt_{a.thr}")).reset_index())
    long.columns = ["sample", "pair", "mean", f"frac_gt_{a.thr}"]
    long = long.dropna(subset=["mean"]).merge(meta.rename_axis("sample").reset_index(), on="sample")
    long.to_csv(out / "sample_pair_summary.csv", index=False)
    print(f"Done. Figures and tables in {out.resolve()}")


if __name__ == "__main__":
    main()
