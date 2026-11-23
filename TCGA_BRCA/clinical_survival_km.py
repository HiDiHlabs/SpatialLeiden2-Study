#!/usr/bin/env python3
"""Kaplan-Meier analysis of progression-free interval (PFI) by tumour stage
for the selected TCGA-BRCA patients.

Each stage grouping is analysed separately:
  ajcc_stage                    collapsed stage (II, III)
  ajcc_pathologic_tumor_stage   substage (IIA, IIB, IIIA, ...)

Per grouping: KM curves of PFI (event = PFI, time = PFI.time in days, shown
in months = days / 30.44) with 95% CI bands, number at risk below the plot
and the log-rank p-value (multivariate log-rank test across all groups);
median PFI per group (KM median); and a clinical summary table per group.

INPUT (under /data/cephfs-2/unmirrored/projects/hgsoc-gen/tcga_brca)
  clinical/TCGA_BRCA_CDR_endpoints_selected.csv   from filter_clinical_selected.py

OUTPUT (in clinical_out/, per grouping <g>)
  km_pfi_<g>.pdf / .png          KM plot
  median_pfi_<g>.csv             group, n, events, median PFI (months)
  clinical_summary_<g>.csv       per group: n, age median [IQR],
                                 histological_type n (%), PFI events n (%),
                                 median follow-up (reverse KM, months)

--------------------------------------------------------------------------
USAGE
--------------------------------------------------------------------------
  cd /data/cephfs-1/work/groups/krieger/users/frka13_c/proj-spatialleiden/SpatialLeiden2-Study
  source /data/cephfs-1/work/groups/krieger/users/frka13_c/miniforge/etc/profile.d/conda.sh && conda activate survival
  python3 TCGA_BRCA/clinical_survival_km.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.plotting import add_at_risk_counts
from lifelines.statistics import multivariate_logrank_test


# --------------------------------------------------------------------------
# STYLE -- edit here to change the look of every figure
# Nature Methods theme: https://research-figure-guide.nature.com/figures/
# --------------------------------------------------------------------------
def mm(x):
    return x / 25.4                                  # mm -> inches


FONT_FAMILY = ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"]  # first installed wins
FONT_SIZE = {"label": 7, "tick": 6, "legend": 6, "annot": 6}  # pt; Nature: 5-7

FIGSIZE = {"km": (mm(89), mm(70))}                   # at-risk table is added below the axes
DPI = 450
FORMATS = ("pdf", "png")

LINE_WIDTH = 0.75                                    # pt; Nature: 0.25-1
AXIS_WIDTH = 0.5
CI_ALPHA = 0.15                                      # transparency of the 95% CI bands
LEGEND_LOC = "upper right"
PVAL_POS = (0.0, 1.03)                               # axes coords of the log-rank p text (above the axes)
AT_RISK_YPOS = -0.3                                  # axes coords of the at-risk table (lifelines default -0.6)

WONG = {
    "black": "#000000", "orange": "#E69F00", "sky_blue": "#56B4E9",
    "green": "#009E73", "yellow": "#F0E442", "blue": "#0072B2",
    "vermillion": "#D55E00", "purple": "#CC79A7",
}
PALETTE = [WONG[k] for k in ("blue", "orange", "green", "vermillion", "purple", "sky_blue", "yellow", "black")]

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": FONT_FAMILY,
    "font.size": FONT_SIZE["label"],
    "axes.titlesize": FONT_SIZE["label"],
    "axes.labelsize": FONT_SIZE["label"],
    "xtick.labelsize": FONT_SIZE["tick"],
    "ytick.labelsize": FONT_SIZE["tick"],
    "legend.fontsize": FONT_SIZE["legend"],
    "legend.frameon": False,
    "lines.linewidth": LINE_WIDTH,
    "axes.linewidth": AXIS_WIDTH,
    "xtick.major.width": AXIS_WIDTH,
    "ytick.major.width": AXIS_WIDTH,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "axes.grid": False,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "savefig.dpi": DPI,
    "pdf.fonttype": 42,                              # editable text in Illustrator/Inkscape
    "ps.fonttype": 42,
})
# --------------------------------------------------------------------------

BASE_DIR = Path("/data/cephfs-2/unmirrored/projects/hgsoc-gen/tcga_brca")
CSV_IN = BASE_DIR / "clinical" / "TCGA_BRCA_CDR_endpoints_selected.csv"
OUT_DIR = BASE_DIR / "clinical_out"

EVENT_COL = "PFI"
TIME_COL = "PFI.time"
DAYS_PER_MONTH = 30.44

# stage column -> label format for its values
STAGE_COLS = {
    "ajcc_stage": "Stage {}",
    "ajcc_pathologic_tumor_stage": "{}",
}
AGE_COL = "age_at_initial_pathologic_diagnosis"
HIST_COL = "histological_type"


def save(fig, stem: Path) -> None:
    for ext in FORMATS:
        fig.savefig(f"{stem}.{ext}", dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"-> {stem}.{{{','.join(FORMATS)}}}")


def fmt_median(x: float) -> str:
    return "not reached" if np.isinf(x) else f"{x:.1f}"


def km_plot(df: pd.DataFrame, col: str, label_fmt: str, groups: list[str]) -> None:
    fig, ax = plt.subplots(figsize=FIGSIZE["km"])
    fitters = []
    for group, color in zip(groups, PALETTE):
        sub = df[df[col] == group]
        kmf = KaplanMeierFitter()
        kmf.fit(sub["months"], sub[EVENT_COL], label=f"{label_fmt.format(group)} (n={len(sub)})")
        kmf.plot_survival_function(ax=ax, ci_show=True, ci_alpha=CI_ALPHA, color=color)
        fitters.append(kmf)

    test = multivariate_logrank_test(df["months"], df[col], df[EVENT_COL])
    ax.text(*PVAL_POS, f"Log-rank p = {test.p_value:.3g}", transform=ax.transAxes,
            fontsize=FONT_SIZE["annot"])
    print(f"{col}: log-rank p = {test.p_value:.3g} (df = {test.degrees_of_freedom})")

    ax.set_xlabel("Time (months)")
    ax.set_ylabel("PFI probability")
    ax.set_ylim(0, 1.05)
    ax.set_xlim(left=0)
    ax.legend(loc=LEGEND_LOC)
    add_at_risk_counts(*fitters, ax=ax, rows_to_show=["At risk"], ypos=AT_RISK_YPOS)
    save(fig, OUT_DIR / f"km_pfi_{col}")


def median_pfi(df: pd.DataFrame, col: str, groups: list[str]) -> None:
    rows = []
    for group in groups:
        sub = df[df[col] == group]
        kmf = KaplanMeierFitter().fit(sub["months"], sub[EVENT_COL])
        rows.append({
            "group": group,
            "n": len(sub),
            "events": int(sub[EVENT_COL].sum()),
            "median_pfi_months": fmt_median(kmf.median_survival_time_),
        })
    out = pd.DataFrame(rows)
    out.to_csv(OUT_DIR / f"median_pfi_{col}.csv", index=False)
    print(out.to_string(index=False))
    print(f"-> {OUT_DIR / f'median_pfi_{col}.csv'}")


def clinical_summary(df: pd.DataFrame, col: str, groups: list[str]) -> None:
    hist_types = sorted(df[HIST_COL].dropna().unique())
    table: dict[str, dict[str, str]] = {}
    for group in groups:
        sub = df[df[col] == group]
        n = len(sub)
        age = sub[AGE_COL]
        events = int(sub[EVENT_COL].sum())
        followup = KaplanMeierFitter().fit(sub["months"], 1 - sub[EVENT_COL]).median_survival_time_
        column = {
            "n": str(n),
            "Age (years), median [IQR]":
                f"{age.median():.0f} [{age.quantile(0.25):.0f}-{age.quantile(0.75):.0f}]",
        }
        for h in hist_types:
            k = int((sub[HIST_COL] == h).sum())
            column[f"{HIST_COL}: {h}, n (%)"] = f"{k} ({100 * k / n:.0f}%)"
        column["PFI events, n (%)"] = f"{events} ({100 * events / n:.0f}%)"
        column["Median follow-up (months, reverse KM)"] = fmt_median(followup)
        table[group] = column
    out = pd.DataFrame(table)
    out.index.name = "variable"
    out.to_csv(OUT_DIR / f"clinical_summary_{col}.csv")
    print(out.to_string())
    print(f"-> {OUT_DIR / f'clinical_summary_{col}.csv'}")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(CSV_IN, index_col=0)
    df["months"] = df[TIME_COL] / DAYS_PER_MONTH
    print(f"{len(df)} patients from {CSV_IN}")

    for col, label_fmt in STAGE_COLS.items():
        groups = sorted(df[col].unique())
        print(f"\n=== {col}: {', '.join(groups)} ===")
        km_plot(df, col, label_fmt, groups)
        median_pfi(df, col, groups)
        clinical_summary(df, col, groups)


if __name__ == "__main__":
    main()
