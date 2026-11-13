# Independent docking verification and molecular-generation peer review

Work in `/app` with the supplied `data/` and write your submission to `output/`.

Review the claims of a molecular generation study evaluated on 100 CrossDocked2020 targets using the supplied benchmark tables and molecular structures in `data/`.

Prepare the CDK2 and ROCK1 receptors from their PDB structures. Justify box centers and extents using the corresponding co-crystal ligands. Dock the two positive controls and fourteen test ligands. Quantify co-crystal redocking RMSD and calculate QED, SA, LogP, TPSA, HBA, HBD, Fsp3 and rotatable bonds for all sixteen ligands. Apply each SR condition from `sr_criteria.txt` and compare your results against the study tables.

1. Determine which methods lead across metrics. Examine one-sided inference, effect sizes, intervals and consistency between aggregate and paired pocket-level results.
2. Quantify ablation effects and interactions. Correct for multiple testing and assess whether synergy generalizes across metrics, timesteps and backbones.
3. Audit the SR criteria. Explain what can be calibrated with the available Vina Score and Vina Dock fields and what is limited by missing data.
4. Evaluate whether the selected case-study molecules support claims of broad applicability. Use the independent results to quantify agreement and discuss selection bias.
5. Identify the least supported central claim and propose one concrete supplementary experiment.

Write `docking_results.csv`, `property_results.csv`, `sr_screening.csv`, `review_report.md`, `analysis_script.py`, and at least six legible figures to `output/`. Include ligand IDs, target IDs and values with units and justified precision. Also submit `review_findings.json` with the following keys:

- Booleans: `method_i_is_full_model`, `method_i_best_vina`, `method_i_not_best_sa_or_lipinski`, `c1_pairwise_calibration_feasible`, `active_fsp3_systematic_failure`, `case_study_selected`, `backbone_independence_complete`.
- An eight-metric direction map `metric_directions` with `lower` or `higher`, plus lists `significant_methods`, `hle_negative_methods` and `fdr_significant_metrics`.
- A method ID `sign_reversal_method`, numeric fields `method_h_aggregate_delta` and `method_h_mpd`, and a `degradation_vina` map keyed by `w/o_S`, `w/o_PQ`, `w/o_V`, `w/o_P` and `w/o_Q`.
- Numeric fields `pq_synergy_vina_score`, `sa_p_gain`, `sa_q_gain`, `sa_both_gain`, `sr_criteria_count`, `full_sr_rate` (fraction), `sa_target`, `sa_actual`, `sa_sr_threshold`, `optimal_timestep`, `cdk2_redocking_rmsd` and `rock1_redocking_rmsd`.
- `most_vulnerable_claim` with either `full_sr` or `multi_metric_synergy`.
