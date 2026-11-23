# 05_evaluation

Compare the single-modality runs (step 03) and the weight-grid runs (step 04).

| Subfolder | Question | Output |
|---|---|---|
| `domain_maps/` | What do the domains look like? | Maps per modality and for selected grid runs |
| `redundancy_ari/` | Does a modality add anything beyond RNA? | ARI of every run against every other run |
| `gain_readouts/` | Does a modality improve agreement with readouts not used for clustering? | Metric per run vs. Kuppe zones, patient group, histology |
| `functional_only_example/` | Is there a domain only functional modalities find? | One worked example: map, defining features, interpretation |

A readout must not have been used to build any modality.
