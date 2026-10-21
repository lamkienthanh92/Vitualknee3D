# Analysis code

`run_all_analyses.py` reproduces every number, table and figure of the manuscript except Fig. 3 (screenshots of the web app). It runs as a single Google Colab cell or as a plain Python 3 script (NumPy, SciPy, pandas, Matplotlib; all preinstalled in Colab).

## Inputs (one folder)

| File | Content | Included here |
|---|---|---|
| `lit_params_v5.json` | Published parameters used by the rotation model (values transcribed from the cited studies) | yes (`data/`) |
| `template_knee_S193761.npz` | Femur and tibia vertices, landmarks and robot passive-flexion path of Dryad specimen S193761 (CC0) | yes (`data/`) |
| `Dataset_DCC_DCTN_1sheet.xlsx` | Vietnamese specimen measurements (22 knees) | no — available from the corresponding author |
| `vn_triplicates.csv` | Three repeated readings per rotation and Lachman measurement | no — available from the corresponding author |

In Colab the script mounts Google Drive and reads `MyDrive/Sim_DCC_DCTN/`. Outside Colab, set the folder with `SIM_BASE=/path/to/folder/ python run_all_analyses.py`.

## Outputs (`outputs/`; a copy of the manuscript run is in `results/`)

- `Table1_vietnamese_vs_reference.csv`, `sectioning_effects.csv`, `reliability.csv` — Part 1 (Table 1, Fig. 1)
- `AdditionalFile1_leave_one_out.csv`, `procedure_outcomes.csv` — Part 2 (Fig. 2c, Additional file 1)
- `AdditionalFile1_geometry_checks.csv`, `ring_sweep.csv`, `insertion_curves.csv`, `tunnel_convergence.csv` — Part 3 (Fig. 2a, b, d)
- `Fig1.tif`, `Fig2.tif` (600 dpi) and PNG copies — Part 4

The random seed is fixed (2026), so the population-model proportions are reproducible.
