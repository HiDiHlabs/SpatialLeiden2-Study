#!/usr/bin/env python3
"""Filter the TCGA clinical sheets down to the patients in selected_sample_ids.txt.

The patient code is the first three dash-separated fields of each slide ID
(TCGA-A1-A0SK-01Z-00-DX1.... -> TCGA-A1-A0SK).

INPUT (under /data/cephfs-2/unmirrored/projects/hgsoc-gen/tcga_brca)
  selected_sample_ids.txt
  clinical/TCGA_BRCA_CDR_endpoints.csv        patient column: case_id
  clinical/TCGA-CDR-SupplementalTableS1.xlsx  patient column: bcr_patient_barcode

OUTPUT (in clinical/)
  TCGA_BRCA_CDR_endpoints_selected.csv
  TCGA-CDR-SupplementalTableS1_selected.xlsx  tabs with a bcr_patient_barcode
                                              column are filtered, all other
                                              tabs are copied unchanged
  missing_patients.txt                        selected patients not found in a
                                              sheet (sheet<TAB>patient), also
                                              printed to the terminal

--------------------------------------------------------------------------
USAGE
--------------------------------------------------------------------------
  cd /data/cephfs-1/home/users/frka13_c/work/proj-spatialleiden/SpatialLeiden2-Study
  conda activate survival
  python3 TCGA_BRCA/filter_clinical_selected.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

BASE_DIR = Path("/data/cephfs-2/unmirrored/projects/hgsoc-gen/tcga_brca")
CLINICAL_DIR = BASE_DIR / "clinical"
SAMPLE_IDS = BASE_DIR / "selected_sample_ids.txt"

CSV_IN = CLINICAL_DIR / "TCGA_BRCA_CDR_endpoints.csv"
CSV_OUT = CLINICAL_DIR / "TCGA_BRCA_CDR_endpoints_selected.csv"
CSV_PATIENT_COL = "case_id"

XLSX_IN = CLINICAL_DIR / "TCGA-CDR-SupplementalTableS1.xlsx"
XLSX_OUT = CLINICAL_DIR / "TCGA-CDR-SupplementalTableS1_selected.xlsx"
XLSX_PATIENT_COL = "bcr_patient_barcode"

MISSING_OUT = CLINICAL_DIR / "missing_patients.txt"


def read_patients(path: Path) -> list[str]:
    ids = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    return list(dict.fromkeys("-".join(i.split("-")[:3]) for i in ids))


def missing_from(patients: list[str], present: pd.Series) -> list[str]:
    found = set(present.astype(str))
    return [p for p in patients if p not in found]


def main() -> None:
    patients = read_patients(SAMPLE_IDS)
    print(f"{len(patients)} selected patients from {SAMPLE_IDS}")
    missing: list[tuple[str, str]] = []

    csv = pd.read_csv(CSV_IN, index_col=0)
    csv_sel = csv[csv[CSV_PATIENT_COL].isin(patients)]
    csv_sel.to_csv(CSV_OUT)
    missing += [(CSV_IN.name, p) for p in missing_from(patients, csv[CSV_PATIENT_COL])]
    print(f"{CSV_IN.name}: {len(csv)} -> {len(csv_sel)} rows -> {CSV_OUT}")

    sheets = pd.read_excel(XLSX_IN, sheet_name=None, header=None)
    with pd.ExcelWriter(XLSX_OUT) as writer:
        for name, sheet in sheets.items():
            header = sheet.iloc[0].tolist() if len(sheet) else []
            if XLSX_PATIENT_COL in header:
                col = header.index(XLSX_PATIENT_COL)
                body = sheet.iloc[1:]
                body_sel = body[body[col].isin(patients)]
                out = pd.concat([sheet.iloc[:1], body_sel])
                missing += [(f"{XLSX_IN.name}:{name}", p) for p in missing_from(patients, body[col])]
                print(f"{XLSX_IN.name}:{name}: {len(body)} -> {len(body_sel)} rows")
            else:
                out = sheet
                print(f"{XLSX_IN.name}:{name}: no {XLSX_PATIENT_COL} column, copied unchanged")
            out.to_excel(writer, sheet_name=name, header=False, index=False)
    print(f"-> {XLSX_OUT}")

    with MISSING_OUT.open("w") as fh:
        fh.write("sheet\tpatient\n")
        for sheet_name, patient in missing:
            fh.write(f"{sheet_name}\t{patient}\n")
    if missing:
        print(f"{len(missing)} selected patient(s) missing from a sheet:")
        for sheet_name, patient in missing:
            print(f"  {sheet_name}\t{patient}")
    else:
        print("All selected patients found in every filtered sheet.")
    print(f"-> {MISSING_OUT}")


if __name__ == "__main__":
    main()
