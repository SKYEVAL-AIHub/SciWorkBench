# Quantitative in vivo calcium imaging

Work in `/app` with the supplied `data/` and write your submission to `output/`.

A short chemical stimulus is applied from 5.5 s to 6.0 s while recording six peripheral sensory neurons. Align each frame to the motion-free reference, identify unreliable frames and separate motion, background, neuropil and slow fluorescence drift from the stimulus response. Use the measured times in `data/timestamps.csv`. Follow the time windows, region statistics, quality thresholds and responder criteria in `data/recording_notes.json`. The input includes `recording.avi`, `reference_frame.tif`, `roi_map.png` and `regions/`; the region-label map is `regions/region_labels.csv`.

Submit exactly these files in `output/`:

- `motion_qc.csv`: 240 rows with `frame,time_s,shift_x_px,shift_y_px,global_mean_ratio,saturated_pixel_fraction,valid`. Report the x/y shifts applied to align each frame. Use `yes` or `no` for `valid`.
- `calcium_traces.csv`: 240 rows with `frame,time_s,cell1,cell2,cell3,cell4,cell5,cell6`. Report the corrected relative fluorescence in percent.
- `summary.csv`: six rows with `cell,F0,baseline_sd_percent,peak_dFF0_percent,peak_time_s,latency_half_peak_s,AUC_percent_s,responder`. Use `yes` or `no` for `responder`.
- `analysis.png`: a readable figure with the reference image and six labeled ROIs, both motion components with rejected frames marked, and six labeled response curves with seconds, percent and the stimulus interval.
