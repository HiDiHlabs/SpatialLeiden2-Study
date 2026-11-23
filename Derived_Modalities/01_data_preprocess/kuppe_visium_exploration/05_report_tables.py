"""Write LaTeX tables for the report from the exploration outputs.

Usage: python 05_report_tables.py <results dir> <report dir>
"""

import sys
from pathlib import Path

import pandas as pd

RES, REP = Path(sys.argv[1]), Path(sys.argv[2])
ZONE_NAME = {
    "CTRL": "control",
    "RZ": "remote",
    "BZ": "border",
    "IZ": "ischaemic",
    "FZ": "fibrotic",
}


def esc(s):
    return str(s).replace("_", r"\_").replace("%", r"\%").replace("&", r"\&")


# per-section table
sec = pd.read_csv(RES / "visium_sections.csv")
lines = [
    r"\begin{longtable}{llllrrrrrc}",
    r"\toprule",
    r"Section & Run & Patient & Zone & Spots & UMIs & Genes & \% mt & Cells & In h5ad \\",
    r"\midrule",
    r"\endhead",
]
for _, r in sec.iterrows():
    lines.append(
        " & ".join(
            [
                esc(r.patient_region_id),
                esc(r.hca_sample_id),
                r.patient,
                ZONE_NAME[r.zone],
                f"{r.n_spots:,}",
                f"{r.median_umi:,.0f}",
                f"{r.median_genes:,.0f}",
                f"{r.median_mt:.0f}",
                f"{r.median_cells:.1f}",
                r.h5ad,
            ]
        )
        + r" \\"
    )
lines += [r"\bottomrule", r"\end{longtable}"]
(REP / "table_sections.tex").write_text("\n".join(lines) + "\n")

# per-zone table
z = pd.read_csv(RES / "depth_by_zone.csv")
lines = [
    r"\begin{tabular}{lrrrrrrr}",
    r"\toprule",
    r"Zone & Sections & Spots & UMIs & Genes & Cells & UMIs per cell & \% mt \\",
    r"\midrule",
]
for _, r in z.iterrows():
    lines.append(
        " & ".join(
            [
                ZONE_NAME[r.zone],
                str(r.n_sections),
                f"{r.n_spots:,}",
                f"{r.median_umi:,.0f}",
                f"{r.median_genes:,.0f}",
                f"{r.median_cells:.2f}",
                f"{r.median_umi_per_cell:,.0f}",
                f"{r.median_mt:.0f}",
            ]
        )
        + r" \\"
    )
lines += [r"\bottomrule", r"\end{tabular}"]
(REP / "table_zones.tex").write_text("\n".join(lines) + "\n")

# download status from the inventory
inv = pd.read_csv(RES / "data_inventory.csv")
inv["status"] = inv.loaded.fillna("").map(
    lambda s: "incomplete" if "incomplete" in s else "complete"
)
g = (
    inv.groupby(["set", "status"])
    .agg(n=("path", "size"), gb=("size_gb", "sum"))
    .reset_index()
)
lines = [
    r"\begin{tabular}{llrr}",
    r"\toprule",
    r"Set & Status & Files & GB \\",
    r"\midrule",
]
for _, r in g.iterrows():
    lines.append(f"{esc(r.set)} & {r.status} & {r.n} & {r.gb:.1f}" + r" \\")
lines += [r"\bottomrule", r"\end{tabular}"]
(REP / "table_download.tex").write_text("\n".join(lines) + "\n")
print("tables written")
