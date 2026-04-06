# Block 6 File Guide

This file explains the exported Block 6 artifacts.

## Top-level files

- `block6_summary.md`: Short Markdown interpretation summary for Block 6.
- `block6_influence_top.csv`: Top influential leave-one-out / leave-last-k-out cases for Block 6.
- `block6_family_sensitivity.csv`: Top-level latency-family sensitivity comparison table for Block 6.
- `block6_bayesian_joint_summary.csv`: Top-level Bayesian joint-uncertainty summary table for Block 6.
- `block6_supported_horizon_overlay.csv`: Top-level overlay comparing model-implied latency medians against Block 4 primary-supported and latest-stable horizons.
- `block6_plot_family_sensitivity.png`: Top-level latency-family sensitivity plot rebuilt from the exported family-sensitivity CSV.
- `block6_plot_bayesian_cure_vs_sigma.png`: Top-level Bayesian cure-vs-sigma scatter plot rebuilt from the exported Bayesian draw CSV.
- `block6_plot_bayesian_cure_vs_uncured_median.png`: Top-level Bayesian cure-vs-uncured-median scatter plot rebuilt from the exported Bayesian draw CSV.
- `README.md`: Top-level guide to the Block 6 exported files.

## Supporting tables

- `block6_baseline_validation.csv`: Supporting-table validation of reconstructed lognormal baseline metrics against saved Block 1 fits.
- `block6_influence_log_normal.csv`: Long-form leave-one-out / leave-last-k-out diagnostics for Block 6.
- `block6_bayesian_joint_draws.csv`: Supporting-table posterior-draw export for Block 6 Bayesian joint uncertainty.

## Technical files

- `block6_code_reference.csv`: Reference table for the Block 6 entry and core analysis scripts.
- `block6_export_manifest.csv`: Manifest of exported Block 6 files.
