# Organic structure elucidation from NMR, NOESY and HRMS

Work in `/app` with the supplied `data/` and write your submission to `output/`.

Infer the structures of compounds A, B and C from `data/compound_A.csv`, `data/compound_B.csv`, `data/compound_C.csv` and `data/compound_C_noesy.png`. Compound C is a tricyclic product of radical cascade cyclization of compound A. The CSVs list proton and carbon NMR peaks (500/126 MHz, CDCl₃) and measured HRMS m/z. Use the `note` field for ion type and equivalent carbons. Use accurate mass and NMR counts to establish composition. Use the phase-sensitive NOESY spectrum and proton assignments to determine the relative configuration of C.

Write only `output/structures.csv` with the header `compound_id,smiles,inchikey` and rows for `compound_A`, `compound_B` and `compound_C`. Derive each standard InChIKey programmatically from the same row's SMILES. Report the constitution of A and B. Specify all six stereocenters of C using SMILES stereochemical marks. Either enantiomer of the correct relative stereoisomer is accepted.
