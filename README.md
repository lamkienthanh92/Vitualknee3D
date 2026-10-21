# Virtual Knee 3D

A browser-based 3D knee for **teaching** and **anterolateral reconstruction planning** (ALL reconstruction / lateral extra-articular tenodesis combined with ACL reconstruction), with a Vietnamese reference anatomy and a patient-scan module that reads knee CT/MRI DICOM directly in the browser.

> **Not a medical device.** This software is for education and research. Outputs are simulations and planning support only; every value must be checked by a qualified clinician. It has not been clinically validated.

## Quick start

- **Online:** enable GitHub Pages for this repository (Settings → Pages → deploy from the `main` branch, root folder) and open the published URL.
- **Offline:** download the repository and double-click `index.html`. An internet connection is needed the first time for the three.js, JSZip and font libraries loaded from public CDNs.
- **Try the scan module:** open the *Patient scan* tab and load `samples/synthetic_knee_CT.zip`, then press *Extract bones from CT* and *Plan on this patient's knee*.

## Interface

Icon-first workspace: mode switcher at the top (Teaching, Graft planning, Vietnamese anatomy, Patient scan); floating, draggable windows on the right (📋 surgical plan, 📈 graft length chart, 🎛️ controls, 🧍 patient factors and flags, 🧱 display); view tools on the left (☰ side panel, 🏷️ labels, ⤢ fit, ▶️ play flexion, 🩻 see-through bone); a flexion slider along the bottom. The full side panel is one click away. In Patient scan, ⊞ shows the best sagittal, coronal, axial and ACL-oblique MRI series side by side; landmarks can be placed in any panel. 🪄 builds the 3D knee automatically: from CT by bone extraction, from MRI by matching the template bone shapes to the images (cortex-contrast score, multi-start pattern search over position, rotation and size, then femur and tibia refined separately; same approach as the sibling shoulder app). The matched outlines are drawn on every slice; a low image-match score is flagged. The MRI result is the template shape at the patient's size and pose, not a segmentation.

## What it does

| Tab | Purpose |
|---|---|
| **Teaching** | Real knee geometry (bones, cartilage, ACL, PCL, MCL, LCL from cadaver MRI) moving through 0–90° with robot-measured kinematics. Toggle ACL / ALL / Kaplan-fibre injuries, Segond and tibial-spine avulsion fractures (types I–III), and ACLR ± ALLR / LET; internal rotation under 5 N·m comes from a population "virtual knee" model. |
| **Graft planning** | Move femoral and tibial fixation points, choose fixation angle, see graft length change across flexion, an isometry map, femoral tunnel convergence between ACL and ALL tunnels, patient-factor indication flags with sources, and print a one-page plan (PDF). Defaults to the Vietnamese reference anatomy. |
| **Vietnamese anatomy** | ALL origin ring at the measured distance from the lateral epicondyle, insertion at the Vietnamese ratio between Gerdy's tubercle and the fibular head, ALL/ITB dimensions and footprints, ALL–lateral meniscus attachment and ITB–ALL fibre continuity. |
| **Patient scan** | Load a DICOM folder (.zip) from a hospital CD/PACS. Series are detected and ranked automatically, identifiers ignored, images viewed with multiplanar reformats. Landmarks can be placed manually; tibial slope is measured from four points. **CT:** bones are extracted automatically (threshold, cavity filling, bone separation, surface nets), labelled, landmarked and fitted to the model at true size, so planning runs on the patient's own knee. **MRI or difficult cases:** load STL files segmented in 3D Slicer or similar. |

A second, simpler tool (`tools/rotation-risk-simulator.html`, Vietnamese interface) shows the population model: the probability that a knee remains loose, is restored, or is over-constrained after ACLR alone, ACLR + ALLR, or ACLR + LET.

## Inputs

| Input | Support |
|---|---|
| DICOM CT or MRI, uncompressed (Explicit/Implicit VR Little Endian), as a .zip of the folder or individual files | Yes |
| Compressed DICOM (JPEG, JPEG 2000, RLE) | Not yet; the app reports these and skips them |
| STL (binary or ASCII), LPS or RAS coordinates | Yes |
| Screenshots, JPEG, PDF reports, photos of films | No (they carry no physical scale) |

Planning reliability is graded from the data: **high** for 3D or CT data at ≤1 mm, **moderate** up to 2 mm, **low** for thick 2D slices.

**Privacy:** all files are processed locally in the browser; nothing is uploaded. Patient name, ID, birth date, institution and physician fields are detected and never displayed or stored.

## How it works (short)

See [`docs/METHODS.md`](docs/METHODS.md) for details and validation.

- **Geometry:** MRI segmentation of cadaver specimen S193761 (Dryad, CC0), left knee mirrored to a right knee; colour surface scans registered by trimmed ICP for display.
- **Kinematics:** passive flexion path of the same knee measured on a 6-DOF robot (Grood–Suntay convention; translations applied in the tibial frame), 0–70°, held constant beyond.
- **Rotational model:** lumped-parameter "virtual knee" (springs in parallel) with Monte Carlo variation, calibrated to robotic sectioning studies.
- **Vietnamese reference anatomy:** aggregate measurements from a cadaver study of 22 Vietnamese knees (see `DATA_SOURCES.md`).
- **Patient bones (CT):** threshold → marrow-cavity filling → two-voxel erosion to separate touching bones → connected components → label regrowth → smoothed surface nets; femur, tibia, fibula, patella and small fragments are labelled by size and position; epicondyles, Gerdy's tubercle and the fibular-head apex are located automatically and fitted rigidly (no scaling) to the model frame.

## Limitations

- One cadaver knee provides the template geometry and kinematics; the ITB and Kaplan fibres are drawn schematically.
- Kinematics are rigid-body and data-driven up to 70°; soft-tissue deformation is not modelled. Graft "strain" is length change, not force.
- Automatic CT bone extraction is threshold-based: osteoporotic bone, very narrow joint spaces, comminution and osteophytes may need manual checking or STL import.
- MRI bones are not segmented automatically yet.
- Fractures in the Teaching tab are simulated patterns, not read from a scan.

## Roadmap

Compressed DICOM decoding · in-app correction tools for segmentation · MRI segmentation model · validation of automatic extraction against expert segmentation (Dice, surface distance) · additional modules (patellar instability, PCL / posterolateral corner, meniscal root repair).

## Data and licences

- **Code:** MIT licence (`LICENSE`).
- **Template knee:** Andreassen T. et al., *Subject-specific knee models, data, and results for specimen S193761*, Dryad, https://doi.org/10.5061/dryad.zcrjdfnpv (CC0). Please cite the original works listed in `DATA_SOURCES.md`.
- **Vietnamese reference anatomy:** see `DATA_SOURCES.md`.
- **Libraries (loaded from CDN, not bundled):** three.js r128 (MIT), JSZip 3.10.1 (MIT), IBM Plex fonts (SIL OFL).
- **Analysis code:** `analysis/run_all_analyses.py` reproduces the manuscript's tables and figures (see `analysis/README.md`).
- **Sample:** `samples/synthetic_knee_CT.zip` is synthetic (generated by `scripts/make_synthetic_ct.py` from the CC0 meshes) and contains no patient data.

## Tóm tắt (Tiếng Việt)

Gối ảo 3D chạy trên trình duyệt, dùng để **giảng dạy** và **lập kế hoạch tái tạo phức hợp trước ngoài** (ALLR/LET kèm tái tạo DCCT). Có bộ giải phẫu tham chiếu người Việt (mặc định khi lập kế hoạch) và thẻ đọc ảnh DICOM CT/MRI ngay trên máy người dùng; với CT, app tự tách xương, tự tìm mốc và lập kế hoạch trên chính gối bệnh nhân. Đây là công cụ giảng dạy và nghiên cứu, **không phải thiết bị y tế**; mọi kết quả cần bác sĩ kiểm tra.

## Archiving on Zenodo (to obtain a DOI)

1. Push this folder to a public GitHub repository.
2. Sign in to zenodo.org with GitHub, open *Account → GitHub* and switch the repository **on**.
3. On GitHub, create a release (for example `v1.3.0`). Zenodo archives it automatically and issues a DOI; metadata come from `.zenodo.json`.
4. Copy the DOI into `CITATION.cff` and into the manuscript (Data availability and Code availability).

Fill in the author names, affiliation and ORCID in `.zenodo.json` and `CITATION.cff` before the first release. Patient images and individual Vietnamese specimen data must never be added to this repository.
