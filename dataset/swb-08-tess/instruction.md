# Flare-aware validation of a periodic dimming signal in AU Mic

Use only the NASA/MAST TESS Sector 1 and 27 products derived in `/workspace/data`. Write all deliverables to `/workspace/output`. Classify the signal as `robust` when the period, event timing, source localization and recovery support it, `ambiguous` when evidence conflicts or is insufficient, or `artifact` when stellar activity, systematics or a neighboring source explains it better.

1. Use BTJD = BJD_TDB − 2457000, preserve the original QUALITY bits, explain masking, and integrate the model over the 20 s and 120 s exposures separately.
2. Use all 32 short-cadence flare windows to define activity criteria. Compare at least two detrending methods that preserve events and distinguish spots, flares and dimmings.
3. Search for periods with BLS or an equivalent open method. Quantify the window function. Adjudicate P/2, P and 2P using coverage, supporting events and contradictions.
4. Fit the midtime, depth and duration with uncertainties for five events using at least 300 bootstrap replicates or equivalent resampling. Fit a linear ephemeris and O−C offsets, and perform leave-one-event-out fits.
5. Estimate five difference-image centroids independently in pixels and arcseconds (21 arcsec per pixel), and assess dilution or contamination with the ICRS neighbor catalog and Tmag flux ratios.
6. Run all 120 cases in `config/injection_grid.csv` and all eight sensitivity scenarios. Derive the final classification from the full evidence chain.

Output BTJD and BJD_TDB midtimes, angles in degrees and original flux in electrons per second. Record uncertainties, masks, thresholds and outlier handling. Use seed `20260819`. The analysis must run on a CPU within two hours using local inputs and open tools.

Write exactly the ten files specified below. In `results.json`, cross-reference the R12–R17 core evidence. In `analysis_report.md`, lead with your conclusion, cite NASA/MAST TESS and Gilbert et al. (2022), DOI [10.3847/1538-3881/ac23ca](https://doi.org/10.3847/1538-3881/ac23ca), and discuss limitations, the BLS box model and correlated noise. Provide six readable diagnostic panels for time series, aliases, O−C, centroids, recovery and sensitivity. Make `analysis_pipeline.py` reproduce the other nine deliverables. In `replay_manifest.json`, record the command, environment, seed, runtime and input/output SHA-256 hashes.

```json
{
  "schema_version": "1.0.0",
  "task_id": "AST-TESS-001",
  "required_files": [
    "analysis_pipeline.py",
    "analysis_report.md",
    "centroid_audit.csv",
    "diagnostic_summary.png",
    "event_metrics.csv",
    "injection_recovery.csv",
    "period_search.csv",
    "replay_manifest.json",
    "results.json",
    "sensitivity.csv"
  ],
  "results_required_keys": [
    "task_id",
    "random_seed",
    "time_system",
    "selected_period_days",
    "classification",
    "event_count",
    "accepted_event_count",
    "centroid_event_count",
    "injection_count",
    "injection_recovery_fraction",
    "catalog_dilution_fraction_42arcsec",
    "linear_ephemeris",
    "leave_one_out_periods",
    "core_chain",
    "runtime_seconds"
  ],
  "classification_enum": ["robust", "ambiguous", "artifact"],
  "csv_required_columns": {
    "event_metrics.csv": ["event_id", "epoch", "sector", "midtime_btjd", "midtime_bjd_tdb", "midtime_error_days", "depth", "depth_error", "duration_days", "exposure_seconds", "oc_seconds", "accepted"],
    "period_search.csv": ["alias_label", "period_days", "bls_snr", "bls_t0_btjd", "bls_duration_days", "bls_depth", "window_function_power", "covered_events", "supported_events", "contradictions", "event_consistency", "adjudication_score", "selected"],
    "centroid_audit.csv": ["event_id", "difference_x_pix", "difference_y_pix", "reference_x_pix", "reference_y_pix", "offset_pix", "offset_arcsec", "difference_snr", "catalog_dilution_fraction_42arcsec", "on_target"],
    "injection_recovery.csv": ["injection_id", "depth", "duration_days", "period_days", "phase_offset_days", "recovered_period_days", "recovered_t0_btjd", "recovered_depth", "detection_snr", "phase_error_days", "recovered"],
    "sensitivity.csv": ["scenario_id", "quality_bits", "flare_sigma", "detrend_days", "exposure_integrated", "selected_period_days", "selected_alias", "supported_events", "classification_stable"]
  },
  "row_count_contract": {
    "event_metrics.csv": 5,
    "period_search.csv": 3,
    "centroid_audit.csv": 5,
    "injection_recovery.csv": 120,
    "sensitivity.csv": 8
  },
  "diagnostic_png_minimum": {"width_px": 1200, "height_px": 700}
}
```
