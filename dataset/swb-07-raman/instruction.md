# Harbor workspace

Work in `/app`. Read the supplied inputs from `/app/data` and write the requested deliverables to `/app/output`. Use the bundled inputs for the analysis.

# Raman spectroscopy of Fe-doped monolayer MoS₂

English translation of the original task instructions. All quantitative results must be calculated independently from the supplied data. Do not substitute literature values for peak positions, magnetic moments or fit parameters.

## Inputs

- `data/sample_comparison.csv` contains comparison spectra for undoped and Fe-doped monolayer MoS₂.
- `data/field_series.csv` contains circularly polarized Raman spectra for Fe-doped MoS₂ in out-of-plane fields from −7 to +7 T.
- `data/temperature_series.csv` contains Raman spectra for the doped sample from approximately 12 to 200 K.

R and L denote right and left circular polarization. The two letters indicate incident and detected scattered polarization, in that order. RL and LR are cross-circular channels. LL is a co-circular channel.

The experiment uses a 633 nm continuous-wave laser in a backscattering geometry. The field series was measured at 12 K. Intensities are in arbitrary units with no common absolute intensity calibration. Positive and negative Raman shifts are retained.

## Objectives

1. Compare the spectra of undoped and Fe-doped monolayer MoS₂ and quantify the main changes due to Fe doping.
2. Analyze the magnetic-field and polarization dependence of distinguishable Raman responses in the doped sample. Extract physically meaningful field-induced frequency shifts and effective phonon magnetic moments. Report key results, uncertainties and evidence for your choices.
3. Analyze the temperature dependence of the low-frequency Raman continuum. Derive a quantitative measure related to the low-frequency static response and assess the reliability of its trend and physical interpretation.
4. Integrate the sample comparison, field dependence and temperature dependence into a consistent physical conclusion. Identify the claims supported directly by these data and the questions requiring further experiments.

## Deliverables

Write the following files under `output/`.

- `raman_analysis.png`, a 300 dpi figure showing the data, fits and quantitative results behind the main conclusions.
- `fit_results.csv`, a table of quantitative results with experimental conditions, uncertainties and units.
- `raman_results.md`, a report describing methods, main quantitative findings, the interpretation across datasets and its evidential limits.


### Harbor output convention

Use the exact column names `component,helicity,metric,value,uncertainty,unit` in
`fit_results.csv` (additional columns are welcome). For the two doped-sample
features, use the component IDs `original_hybridized_mode` and
`dopant_induced_mode`. Give one `center` row per component with unit `cm^-1`,
and one `mu_eff` row per component and channel (`RL` and `LR`) with unit
`mu_B`. Report positive uncertainty estimates for these quantities. Include
your remaining results and the scientific interpretation as requested above.
