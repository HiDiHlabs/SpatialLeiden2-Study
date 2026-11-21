"""Helper functions for comparing spatial-domain methods on the coronal ABC atlas.

Import in a notebook (run from SpatialLeiden2-Study/):

    import matplotlib.pyplot as plt
    from abc_coronal_comparison import *

    ccf, parcellation = load_ccf()                       # CCF annotation per cell

    # one entry per benchmarked tool (column name = where the tool stored its domains)
    methods = {
        "SpatialLeiden": load_labels("results/.../scvi_r1.0.parquet", domain_col="spatialleiden"),
        "GraphST":       load_labels("results/.../graphst.h5ad", domain_col="mclust"),
        "STAligner":     load_labels("results/.../staligner.h5ad", domain_col="mclust"),
        "SEDR":          load_labels("results/.../sedr.h5ad", domain_col="SEDR"),
    }
    labels = methods["GraphST"]

    # 1. heatmap: domains x CCF parcellation index
    fig, ax = plt.subplots(figsize=(90 * mm, 70 * mm))
    plot_ccf_heatmap(ax, labels, ccf, parcellation)

    # 2. example sections: CCF (left) vs domains (right)
    plot_mirror_sections(labels, ccf, parcellation)

    # 3. continuity: domain proportion per section
    fig, ax = plt.subplots(figsize=(90 * mm, 30 * mm))
    plot_continuity(ax, labels)

Every plot function takes an optional `palette`. Use `match_palette` to give a
method the colours of a reference method, so the same region gets the same colour.

Expected input per method: one row per cell, indexed by `cell_label`, with a
domain column, `brain_section_label`, `x` and `y`.
"""

from pathlib import Path

import colorcet as cc
import matplotlib as mpl
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

ABC_DIR = Path("data") / "ABC_atlas"
CCF_FILES = {"Zhuang-ABCA-1": "CCF.tsv", "Zhuang-ABCA-2": "CCF_Zhuang-ABCA-2.tsv"}

# plotting constants (same as utils.py)
mm = 0.1 / 2.54
fontsize = 5
fontsize_l = 7
mpl.rcParams.update({
    "axes.titlesize": fontsize_l,
    "axes.labelsize": fontsize_l,
    "xtick.labelsize": fontsize,
    "ytick.labelsize": fontsize,
    "legend.fontsize": fontsize,
    "legend.title_fontsize": fontsize,
    "pdf.fonttype": 42,
})


# -----------------------------------------------------------------------------
# Loading
# -----------------------------------------------------------------------------
def load_ccf(datasets=("Zhuang-ABCA-1",)):
    """CCF parcellation per cell (+ division, structure, ...) and the parcellation table."""
    ccf = pd.concat([
        pd.read_table(
            ABC_DIR / CCF_FILES[ds], index_col="cell_label",
            usecols=["cell_label", "parcellation_index"],
            converters={"cell_label": str},  # IDs look numeric, keep them as strings
        )
        for ds in datasets
    ])
    parcellation = pd.read_table(ABC_DIR / "parcellation.tsv", index_col="parcellation_index")
    ccf = ccf.join(parcellation, on="parcellation_index")
    return ccf, parcellation


def load_labels(src, domain_col="spatialleiden", datasets=("Zhuang-ABCA-1",), spatial_key="spatial"):
    """Read one tool's result and return columns
    domain, brain_section_label, section (int), x, y — coronal cells only.

    src can be:
      - an AnnData object or .h5ad file (e.g. GraphST / STAligner / SEDR output):
        domains from adata.obs[domain_col], coordinates from adata.obsm[spatial_key]
      - a .parquet / .csv / .tsv file or DataFrame indexed by cell_label with
        columns domain_col, brain_section_label, x, y
    """
    if isinstance(src, (str, Path)) and Path(src).suffix == ".h5ad":
        import anndata as ad
        src = ad.read_h5ad(src, backed="r")
    if hasattr(src, "obs") and hasattr(src, "obsm"):  # AnnData
        df = src.obs[[domain_col, "brain_section_label"]].copy()
        df[["x", "y"]] = np.asarray(src.obsm[spatial_key])[:, :2]
    elif isinstance(src, pd.DataFrame):
        df = src.copy()
    else:
        path = Path(src)
        if path.suffix == ".parquet":
            df = pd.read_parquet(path)
        else:
            df = pd.read_csv(path, sep="\t" if path.suffix == ".tsv" else ",", index_col=0, dtype={0: str})
    df.index = df.index.astype(str)

    df = df.rename(columns={domain_col: "domain"})[["domain", "brain_section_label", "x", "y"]]
    sections = df["brain_section_label"].astype(str)
    df = df[sections.str.rsplit(".", n=1).str[0].isin(datasets)].copy()
    df["section"] = sections.str.extract(r"(\d+)$", expand=False).astype(int)
    df["domain"] = pd.Categorical(df["domain"], categories=sorted(df["domain"].unique()))
    return df


# -----------------------------------------------------------------------------
# Metrics
# -----------------------------------------------------------------------------
def domain_enrichment(labels, anno, anno_key="parcellation_index"):
    """Enrichment (FC) of each annotation category in each domain:
    fraction of the domain's cells in the category / fraction of all cells in it."""
    n = labels[["domain"]].join(anno[[anno_key]], how="inner").groupby(["domain", anno_key], observed=True).size()
    domain_frac = n / n.groupby("domain", observed=True).transform("sum")
    total_frac = n.groupby(anno_key).transform("sum") / n.sum()
    return (domain_frac / total_frac).rename("FC").reset_index()


def ccf_agreement(labels, ccf, level="division"):
    """ARI / NMI between domains and a CCF level (cells outside the brain excluded)."""
    df = labels[["domain"]].join(ccf[[level, "parcellation_index"]], how="inner")
    df = df[df["parcellation_index"] != 0]
    a, b = df[level].astype(str), df["domain"].astype(str)
    return {"ARI": adjusted_rand_score(a, b), "NMI": normalized_mutual_info_score(a, b)}


# -----------------------------------------------------------------------------
# Colours
# -----------------------------------------------------------------------------
def domain_palette(labels):
    cats = labels["domain"].cat.categories
    return dict(zip(cats, sns.color_palette(cc.glasbey, n_colors=len(cats))))


def division_palette(parcellation):
    return dict(zip(parcellation["division"].unique(), cc.glasbey_dark))


def match_palette(labels, ref_labels, ref_palette=None):
    """Colour `labels`' domains like the best-overlapping domain of `ref_labels`."""
    ref_palette = ref_palette or domain_palette(ref_labels)
    common = labels.index.intersection(ref_labels.index)
    overlap = pd.crosstab(ref_labels.loc[common, "domain"], labels.loc[common, "domain"])
    rows, cols = linear_sum_assignment(-overlap.to_numpy())
    palette = {overlap.columns[c]: ref_palette[overlap.index[r]] for r, c in zip(rows, cols)}

    extra = iter(sns.color_palette(cc.glasbey, n_colors=len(ref_palette) + len(labels["domain"].cat.categories))[len(ref_palette):])
    return {d: palette[d] if d in palette else next(extra) for d in labels["domain"].cat.categories}


def add_legend(ax, palette, title=None, ncol=2, **kwargs):
    handles = [
        mlines.Line2D([], [], color=c, marker="o", linestyle="None", markersize=3, label=str(d))
        for d, c in palette.items()
    ]
    kwargs = {"loc": "center left", "bbox_to_anchor": (1.0, 0.5), "frameon": False,
              "handletextpad": 0.3, "columnspacing": 1, **kwargs}
    return ax.legend(handles=handles, title=title, ncol=ncol, **kwargs)


# -----------------------------------------------------------------------------
# 1. Heatmap: domains x CCF parcellation index
# -----------------------------------------------------------------------------
def plot_ccf_heatmap(ax, labels, ccf, parcellation, vmax=10, cbar=True, title=None):
    """Enrichment heatmap. Columns are ordered by CCF division, then parcellation
    index (same order for every method); rows are sorted to give a diagonal.
    A division colour strip is drawn above the heatmap."""
    enr = domain_enrichment(labels, ccf, "parcellation_index")
    matrix = enr.pivot(index="domain", columns="parcellation_index", values="FC").fillna(0)

    div_rank = {d: i for i, d in enumerate(parcellation["division"].unique())}
    cols = sorted(matrix.columns, key=lambda p: (div_rank[parcellation.loc[p, "division"]], p))
    matrix = matrix[cols]
    rows = np.lexsort((-matrix.max(axis=1).to_numpy(), matrix.to_numpy().argmax(axis=1)))
    matrix = matrix.iloc[rows]

    sns.heatmap(matrix, ax=ax, vmin=0, vmax=vmax, cmap="rocket_r", xticklabels=False, yticklabels=True,
                cbar=cbar, cbar_kws={"label": "Enrichment score"}, rasterized=True)
    ax.set(xlabel="CCF (parcellation index)", ylabel="Domain", title=title)
    ax.tick_params(axis="y", length=0, rotation=0)
    for spine in ax.spines.values():
        spine.set_visible(True)

    # division colour strip above the heatmap
    div_pal = division_palette(parcellation)
    strip = ax.inset_axes([0, 1.01, 1, 0.04])
    rgb = [mpl.colors.to_rgb(div_pal[parcellation.loc[p, "division"]]) for p in cols]
    strip.imshow([rgb], aspect="auto", interpolation="nearest")
    strip.set_axis_off()
    if title:
        ax.set_title(title, pad=8)
    return matrix


# -----------------------------------------------------------------------------
# 2. Example sections: CCF annotation (left) vs domains (right)
# -----------------------------------------------------------------------------
def plot_mirror_section(ax, labels, ccf, parcellation, section, palette=None, s=0.05):
    """One section drawn twice: CCF division on the left, mirrored domains on the right.
    `section` is a brain_section_label (e.g. "Zhuang-ABCA-1.064") or a section number."""
    palette = palette or domain_palette(labels)
    div_pal = division_palette(parcellation)
    col = "brain_section_label" if isinstance(section, str) else "section"
    df = labels[labels[col] == section].join(ccf[["division"]], how="inner")

    x = df["x"] - df["x"].max()
    x = x + x.min() * 0.05  # small gap between the halves
    y = -df["y"]
    ax.scatter(x, y, s=s, c=df["division"].map(div_pal).tolist(), linewidth=0, rasterized=True)
    ax.scatter(-x, y, s=s, c=df["domain"].astype(object).map(palette).tolist(), linewidth=0, rasterized=True)
    ax.axvline(0, color="black", linestyle="dashed", lw=0.6)
    ax.set(aspect=1)
    ax.axis("off")


def plot_mirror_sections(labels, ccf, parcellation, sections=range(20, 140, 10), palette=None,
                         nrows=3, title=None):
    """Grid of mirror sections for one method, with domain legend."""
    palette = palette or domain_palette(labels)
    sections = [s for s in sections if (labels["section"] == s).any()]
    ncols = int(np.ceil(len(sections) / nrows))

    fig, axs = plt.subplots(nrows, ncols + 1, figsize=(180 * mm, 80 * mm),
                            gridspec_kw={"width_ratios": [1] * ncols + [0.8], "wspace": 0.05})
    for ax, sec in zip(axs[:, :-1].flat, sections):
        plot_mirror_section(ax, labels, ccf, parcellation, sec, palette)
        ax.set_title(sec, fontsize=fontsize, pad=1)
    for ax in axs[:, :-1].flat[len(sections):]:
        ax.axis("off")
    axs[0, 0].text(0.0, 1.0, "CCF", transform=axs[0, 0].transAxes, fontsize=fontsize, va="top")
    axs[0, 0].text(1.0, 1.0, "Domains", transform=axs[0, 0].transAxes, fontsize=fontsize, va="top", ha="right")

    for ax in axs[:, -1]:
        ax.axis("off")
    add_legend(axs[1, -1], palette, title="Domain", ncol=3, loc="center")
    if title:
        fig.suptitle(title)
    return fig


# -----------------------------------------------------------------------------
# 3. Continuity: domain proportion per section
# -----------------------------------------------------------------------------
def plot_continuity(ax, labels, palette=None, title=None):
    """Stacked bars of domain proportion per section (anterior -> posterior)."""
    palette = palette or domain_palette(labels)
    frac = pd.crosstab(labels["section"], labels["domain"], normalize="index")
    bottom = np.zeros(len(frac))
    for d in frac.columns:
        ax.bar(frac.index, frac[d], 1, bottom=bottom, color=palette[d], linewidth=0)
        bottom += frac[d].to_numpy()
    ax.set(xlim=(frac.index.min() - 0.5, frac.index.max() + 0.5), ylim=(0, 1),
           xlabel="Section", ylabel="Domain proportion", title=title)
    ax.set_yticks([0, 0.5, 1])
    sns.despine(ax=ax)
