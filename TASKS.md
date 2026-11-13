# English task guide and data availability

## SWB-01 · Soma-free 3D neurite tracing

**Inputs and task.** A fluorescence TIFF volume contains distal neuronal processes without a soma. Select seeds, account for anisotropic voxels, trace branches and export an SWC skeleton.

**Grading.** SWC and graph validity (12), geometric path agreement (35), global topology (30), path length and subtree completeness (15), and local topology (8). A visually close trace can still score poorly if it creates false connections between adjacent fibers.

**Included materials.** Includes the original TIFF volume and reference SWC. The automated verifier measures physical-space coverage, branch matches, endpoint and junction agreement, and length. The original rubric adds scientific topology criteria for manual review. [English instruction](dataset/swb-01-neurite/instruction.md).

## SWB-02 · AGO2 perturbation and RNA kinetics

**Inputs and task.** Compare whole, nascent and existing RNA in a compact guide-level pseudobulk dataset derived from PerturbSci-Kinetics. Choose the analysis unit and controls, test whether the AGO2 response is specific, and assess whether the evidence supports direct global transcriptional repression.

**Grading.** Layer dependence, normalization and gene-level attribution (40), internal counterfactuals and specificity (35), and QC, robustness, dimensionality reduction and consistency (25).

**Included materials.** Includes the original guide-level H5AD, metadata, provenance and reference tables. The automated verifier compares aligned numeric and categorical estimates with the frozen reference. Use the original scientific rubric when assessing alternative estimators or claims about causality. Upstream source: [GEO GSE218566](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE218566). [English instruction](dataset/swb-02-ago2/instruction.md).

## SWB-03 · In vivo calcium imaging

**Inputs and task.** Use fluorescence video, measured frame timestamps, a reference frame and ROIs to align images, detect unreliable frames and separate motion, background, neuropil and intensity drift from stimulus responses.

**Grading.** Output schemas (10), motion and abnormal frames (18), time series for six cells (32), response metrics (25), responder calls (5), and integrated visualization (10).

**Included materials.** Includes the original recording, frame timestamps, regions and reference outputs. The automated verifier checks table layout, motion, fluorescence traces, summary metrics, responder calls and PNG integrity. Inspect figure content separately when applying the full original rubric. [English instruction](dataset/swb-03-calcium/instruction.md).

## SWB-04 · Dual-gate FLIM cardiac electrophysiology

**Inputs and task.** Analyze simulated photon counts representing a dual-gated SPAD acquisition. Estimate time-varying fluorescence lifetimes, calibrate membrane potentials, screen low-confidence pixels, align action potentials and quantify repolarization across fields of view.

**Grading.** Quantitative lifetime, voltage and action-potential results (55), adaptive denoising and pixel screening (25), reproducible outputs (10), and biological and physical interpretation (10).

**Included materials.** Includes all files available in the original export. The photon-count table contains sparse coordinate/frame samples. The Harbor adaptation scores coverage, output structure and recognition of unestimable results. The original illustrative reference files do not establish full-field ground truth. [English instruction](dataset/swb-04-flim/instruction.md).

## SWB-05 · Organic structure elucidation with NMR, NOESY and HRMS

**Inputs and task.** Combine NMR peak tables, accurate mass and a phase-sensitive NOESY image to infer three molecular structures. The third compound requires a relative stereochemical assignment.

**Grading.** Formula (15), ring counts (15), functional groups (15), constitutional connectivity (30), and stereochemistry for the third compound (25).

**Included materials.** Includes the original spectral tables, NOESY image and reference structures. The verifier reuses the original RDKit-based score for formula, ring count, functional groups, constitution and relative stereochemistry. Upstream source: [associated article and supporting information](https://doi.org/10.1002/ajoc.202400242). [English instruction](dataset/swb-05-structure/instruction.md).

## SWB-06 · Molecular generation claims and docking verification

**Inputs and task.** Re-dock positive controls, compare recalculated results with the reported benchmark, examine significance and ablation claims, and assess missing data and selection bias in a mock peer review.

**Grading.** Docking verification (8), performance and statistical analysis (26), ablations and interactions (28), success criteria and data limitations (24), and critical synthesis (14).

**Included materials.** Includes all original structure, ligand and study-table inputs. The verifier uses structured findings and report checks from the existing Harbor conversion. The full scientific rubric also assesses analysis, docking methods and evidence quality. The bundled protein structures are [PDB 1H00](https://www.rcsb.org/structure/1H00) and [PDB 6E9W](https://www.rcsb.org/structure/6E9W), subject to the [PDB usage policy](https://www.rcsb.org/pages/usage-policy). [English instruction](dataset/swb-06-docking/instruction.md).

## SWB-07 · Raman spectroscopy of Fe-doped monolayer MoS₂

**Inputs and task.** Compare doped and undoped spectra, follow polarization-resolved branches across magnetic fields, estimate an effective phonon magnetic moment, and interpret temperature-dependent low-frequency scattering.

**Grading.** Sample comparison (20), magnetic and polarization response (25), magnetic moment (15), low-frequency temperature response (15), interpretation and evidence limits (20), and output consistency (5).

**Included materials.** Includes experimental inputs, the English prompt and rubric, original materials and reference fits. The three CSV files in [environment/data](dataset/swb-07-raman/environment/data/) are byte-identical to the archived files below. Only the filenames were changed. [English instruction](dataset/swb-07-raman/instruction.md).

| Original file | Task input |
| --- | --- |
| Fig.2 - Raman spectra MoS2 vs FeMoS2.csv | `sample_comparison.csv` |
| Fig.3 - Field Dependent Raman.csv | `field_series.csv` |
| Fig.4 - Temperature Dependent Raman.csv | `temperature_series.csv` |

**Data attribution.** Wencan Jin and Rui He (2025), *Dilute Magnetic Impurity-Induced Effective Phonon Magnetic Moment in Fe-doped Monolayer MoS₂*, v1, [Zenodo dataset](https://doi.org/10.5281/zenodo.16867939), [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). See also the [associated paper](https://doi.org/10.1088/2053-1583/ae04fc).

## SWB-08 · Flare-aware transit validation with TESS

**Inputs and task.** Use AU Mic Sector 1 and 27 light curves, event and flare windows, target-pixel summaries, apertures and neighbor data. Test the periodic dimming signal using event timing, source localization, sensitivity to stellar activity and injection–recovery completeness.

**Included materials.** Includes MAST-derived inputs, the task-specific verifier and a reference pipeline. See the [provenance manifest](dataset/swb-08-tess/environment/data/provenance_manifest.json) and [MAST data-use policy](https://archive.stsci.edu/publishing/data-use). [English instruction](dataset/swb-08-tess/instruction.md).

## SWB-09 · Quantum reaction rates with ring-polymer surface hopping

**Inputs and task.** Reproduce the adiabatic, nonadiabatic and bead-approximation rates at 143 K and 250 K from a specified one-dimensional two-state model. Compute open-boundary scattering, thermal averages, ring-polymer sampling and nonadiabatic trajectories, and quantify statistical uncertainty.

**Included materials.** Includes the model, scientific specification, verifier and frozen reference calculations. The [English scientific specification](dataset/swb-09-tully/environment/data/scientific_spec.md) paraphrases Shushkov, Li and Tully (2012). The oracle rebuilds deliverables from frozen calculations. [English instruction](dataset/swb-09-tully/instruction.md).

## SWB-10 · Interface-localized electronic modulation in GdTe₂

**Inputs and task.** Analyze STM/STS topography, conductance maps and point spectra. Separate spatial modulation components, locate the curved domain boundary, compare sampled biases and validate the anomaly across forward and retrace scans.


**Included materials.** Includes the original STM/STS maps, spectra, task-specific grader and reference workflow. [English instruction](dataset/swb-10-gdte2/instruction.md).

