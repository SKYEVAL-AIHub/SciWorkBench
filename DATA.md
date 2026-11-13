# Data sources and reuse

The MIT license covers SciWorkBench-authored code and documentation. Bundled
datasets, images, scientific reference outputs and other third-party material
are outside its scope. Upstream terms apply to third-party files. For bundled
scientific assets without a stated reuse license, contact the relevant rights
holder before redistributing or reusing them beyond applicable permissions.

| Task | Bundled materials and source |
| --- | --- |
| SWB-01 | Fluorescence TIFF and reference SWC from the original benchmark task. |
| SWB-02 | Guide-level pseudobulk derived from PerturbSci-Kinetics, [GEO GSE218566](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE218566). The [provenance table](dataset/swb-02-ago2/environment/data/source_provenance.tsv) records the original download and checksum. |
| SWB-03 | In vivo recording, timestamps, reference frame and ROI annotations from the original benchmark task. |
| SWB-04 | Simulated photon counts and illustrative reference files from the original benchmark export. The inputs contain sparse samples. |
| SWB-05 | Spectral tables and NOESY image supplied with the benchmark task. See the [associated article and supporting information](https://doi.org/10.1002/ajoc.202400242). |
| SWB-06 | PDB entries [1H00](https://www.rcsb.org/structure/1H00) and [6E9W](https://www.rcsb.org/structure/6E9W) are covered by the [wwPDB CC0 policy](https://www.rcsb.org/pages/usage-policy). Ligand files and study tables come from the benchmark task. |
| SWB-07 | Wencan Jin and Rui He, *Dilute Magnetic Impurity-Induced Effective Phonon Magnetic Moment in Fe-doped Monolayer MoS₂*, [Zenodo v1](https://doi.org/10.5281/zenodo.16867939), [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). The three CSVs were renamed without modifying their contents. See [the filename mapping](TASKS.md#swb-07--raman-spectroscopy-of-fe-doped-monolayer-mos). |
| SWB-08 | TESS observations and catalog products for AU Mic from MAST. The [provenance manifest](dataset/swb-08-tess/environment/data/provenance_manifest.json) records the source products and checksums. Follow the [MAST data-use and acknowledgment guidance](https://archive.stsci.edu/publishing/data-use). |
| SWB-09 | Model specification based on Shushkov, Li and Tully (2012), with benchmark-generated frozen reference calculations. See the [scientific specification](dataset/swb-09-tully/environment/data/scientific_spec.md). |
| SWB-10 | Original STM/STS maps, topography and spectra supplied with the benchmark task. |
