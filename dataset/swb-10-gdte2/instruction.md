# Interface-localized electronic modulation in GdTe₂

Use low-temperature STM/STS data from two neighboring domains with distinct periodic modulations. Derive real and reciprocal coordinates, detect the interface and test for an energy-selective local modulation. Derive the boundary, active biases and anomalous point spectrum from the data.

Inputs are read-only in `/app/data`: `metadata.json`; `region1_topography.npz` (40 nm × 40 nm); `map_E001.npz` through `map_E004.npz` (same field of view at four biases, see metadata); and `spectra_set_A.csv` (P01/P02/P03, bias, in-phase conductance proxy, quadrature and current). Arrays contain `forward`, `retrace`, `x_nm` and `y_nm`. The retrace x-axis has already been flipped to match the forward direction. Check residual registration. x increases by column, y by stored row. Report wavevectors in cycles/nm, with rad/nm only as an additional conversion. Round wavevector components to 3 decimals, lengths to 2, biases to 1 and dimensionless scores to 2.

1. Check units, pixel scale and forward/retrace registration. Describe preprocessing.
2. Use the topographic FFT to report at least one common modulation vector and two side-selective vectors per side.
3. Use a local Fourier or lock-in amplitude, or an equivalent complex-amplitude method, to segment the domains and fit x_b(y). Report x at y = 5.00, 20.00 and 35.00 nm and a localization uncertainty.
4. Use a common interface and anomaly statistic for all maps. Identify enhanced sampled biases without interpolating boundaries between them. Report the spatial FWHM or an explicitly defined equivalent width.
5. Straighten the interface into (u,y), compare active and control Fourier power, and report the main feature's q_u, q_y, magnitude and angle.
6. Identify the spectrum with a localized positive-bias feature relative to the average of the other two. Use reproducible baseline and half-maximum rules to report its peak, edges and width.
7. Apply the same criteria to forward and retrace scans and compare active biases and fitted boundary locations.

Write `analysis_report.md`, `summary.json`, `analysis.py` and the figures `reciprocal_space.png`, `domain_segmentation.png`, `bias_dependence.png`, `boundary_feature_fft.png`, `spectroscopy_comparison.png` and `validation_retrace.png` to `/app/output`. Include numerical and visual evidence, formulas, uncertainty, assumptions and limitations. The analysis script must recreate the results offline using the supplied inputs. Equivalent scientific methods are accepted when defined reproducibly. Use bias strings with one decimal as keys in the two per-bias maps.

```json
{
  "schema_version": "1.0",
  "reciprocal_convention": "cycles/nm",
  "common_modulation_vector": {"qx": null, "qy": null, "magnitude": null},
  "left_domain_vectors": [{"qx": null, "qy": null, "magnitude": null}],
  "right_domain_vectors": [{"qx": null, "qy": null, "magnitude": null}],
  "boundary": {
    "sample_y_nm": [5.0, 20.0, 35.0],
    "sample_x_nm": [null, null, null],
    "uncertainty_nm": null
  },
  "boundary_anomaly": {
    "detected": null,
    "active_sampled_biases_mV": [],
    "per_bias_boundary_excess_z": {},
    "active_spatial_fwhm_nm": {}
  },
  "fourier_feature": {
    "q_u_cycles_per_nm": null,
    "q_y_cycles_per_nm": null,
    "magnitude_cycles_per_nm": null,
    "angle_deg_from_u": null
  },
  "spectroscopy": {
    "anomalous_point_id": null,
    "feature_peak_mV": null,
    "left_edge_mV": null,
    "right_edge_mV": null,
    "width_mV": null
  },
  "retrace_validation": {
    "forward_active_biases_mV": [],
    "retrace_active_biases_mV": [],
    "boundary_forward_retrace_rmse_nm": null
  }
}
```
