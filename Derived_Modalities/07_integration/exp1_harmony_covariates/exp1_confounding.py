"""Experiment 1, step 0: confounding audit of the candidate Harmony covariates.

Crosstabs of Visium batch, zone, pathologist group, patient and sample, counted in sections
(one row per Visium run in metadata-Visium.csv; n = 28 incl. both RZ_P6 runs) and in spots
(the authors' spot tables of the 27 runs in the processed h5ads). Column names follow the authors:
    visium_batch = batch (10X / ACH), zone = major_labl, sample = hca_sample_id.
Usage (from Derived_Modalities/): python 07_integration/exp1_harmony_covariates/exp1_confounding.py \
    data/kuppe_2022 results/kuppe_2022/exp1_harmony_covariates
"""
import sys
from pathlib import Path

import pandas as pd
from sklearn.metrics import adjusted_mutual_info_score as ami
from scipy.stats import chi2_contingency

data, out = Path(sys.argv[1]), Path(sys.argv[2])
REN = {"batch": "visium_batch", "major_labl": "zone", "hca_sample_id": "sample"}
meta = pd.read_csv(data / "visium_raw" / "metadata-Visium.csv").rename(columns=REN)
spots = pd.concat(
    [pd.read_csv(data / "visium_tables" / f"{s}_spots.csv",
                 usecols=["batch", "major_labl", "patient_group", "patient", "hca_sample_id"])
     for s in meta["sample"] if s != "ACH0021"]).rename(columns=REN)
PAIRS = [("visium_batch", "zone"), ("visium_batch", "patient_group"), ("visium_batch", "patient"),
         ("patient", "zone"), ("sample", "zone"), ("zone", "patient_group")]
rows, summ = [], []
for unit, df in [("sections", meta), ("spots", spots)]:
    for a, b in PAIRS:
        ct = pd.crosstab(df[a], df[b])
        if unit == "sections" or a != "sample":
            print(f"--- {unit}: {a} x {b}\n{ct}\n")
        long = ct.stack().rename("n").reset_index().rename(columns={a: "row_value", b: "col_value"})
        rows.append(long.assign(unit=unit, row_var=a, col_var=b))
        chi2, p, *_ = chi2_contingency(ct)
        v = (chi2 / (ct.values.sum() * (min(ct.shape) - 1))) ** 0.5
        summ.append({"unit": unit, "row_var": a, "col_var": b, "cramers_v": v, "chi2_p": p,
                     "AMI": ami(df[a], df[b])})
pd.concat(rows)[["unit", "row_var", "col_var", "row_value", "col_value", "n"]].to_csv(
    out / "exp1_confounding.csv", index=False)
s = pd.DataFrame(summ)
s.to_csv(out / "exp1_confounding_summary.csv", index=False)
print(s.round(3).to_string(index=False))
