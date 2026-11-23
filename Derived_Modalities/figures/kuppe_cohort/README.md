# Kuppe et al. 2022 cohort

![Kuppe 2022 cohort](kuppe_cohort.png)

**Kuppe et al. 2022 cohort.** **a**, Samples contributed by each of 23 patients (rows) across five myocardial zones. Filled markers: matched snRNA-seq, snATAC-seq and Visium; open markers: snRNA-seq only (external validation). ×2, two samples of one zone. Bottom: samples per zone, and those with Visium. **b**, Time from symptom onset to sampling. Dashed line, 15 days. n.a., control donors. Data from Kuppe et al., Nature 608, 766–777 (2022), Supplementary Table 1.

## Files

| File | Content |
|---|---|
| `kuppe_cohort.pdf` | Vector figure, 183 × 96 mm, Arial 5.5–8 pt (Nature format) |
| `kuppe_cohort.png` | Same figure at 450 dpi, for slides |
| `cohort_table.csv` | One row per sample: patient, zone, patient group, disease, age, sex, days after onset, data types, Visium runs and depth |

## Numbers, checked against the paper

- **Patients:** 23 = 4 control donors + 19 patients (7 acute MI, 9 transplant recipients, 3 validation-only).
- **Samples:** 31. Per zone: control 4, remote 6, border 3, ischaemic 12, fibrotic 6. This matches the paper's Fig. 1a.
- **Visium:** 28 runs. After merging the two runs of patient P6's remote zone, there are 27 sections: 4 / 5 / 3 / 9 / 6 per zone. This matches the per-zone counts in the paper's figure legends.

## How to rebuild

From `Derived_Modalities/`:

```
python 01_data_preprocess/kuppe_cohort/cohort_figure.py data/kuppe_2022 results/kuppe_2022/cohort
```

The script needs the Kuppe metadata sheets (`downloads/kuppe_2022.py --sets metadata visium_raw`) and Supplementary Table 1. The table is not on Zenodo; the download link is in the script's docstring. Visium depth columns are filled only if the Visium exploration (`01_data_preprocess/kuppe_visium_exploration/run_all.sh`) has been run.
