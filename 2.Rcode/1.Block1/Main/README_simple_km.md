# KM Output README

## Overview

This folder contains Kaplan-Meier benchmark outputs generated from the merged preprocessed CHR-P dataset.

Current default export root:

- `/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis of CHR-P Using a Mixture Cure Model/1.Block1/1.KM`

Analysis input:

- Source CSV: `/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis of CHR-P Using a Mixture Cure Model/0.Data/2.Preprocessed data/Preprocessed_Merged_PNUH_SNUH_Data.csv`
- Cohorts analyzed: `PNU`, `SNU`, and `merged`
- Time scale: `time_year = days_followup / 365.25`
- Main event definition: `status_num == 1` is treated as `transition`
- Censoring for the main KM analysis: `status_num %in% c(0, 2)` where `0 = right censoring` and `2 = remission`

The fitted KM objects are saved so the model does not need to be re-fit every time downstream summaries or plots are regenerated.

## Main Files

### Tables

- `simple_km_yearly_estimates.csv`
  - Annual KM summaries for each dataset at `1, 2, ..., 10` years.
- `simple_km_dataset_summary.csv`
  - Simple cohort-level counts and maximum follow-up summaries.
- `simple_km_plot_manifest.csv`
  - Registry of generated plot files and matching RDS objects.

### Fitted Object

- `simple_km_fit_objects.rds`
  - A named list with entries `PNU`, `SNU`, and `merged`.
  - Each entry is a `survival::survfit` object fit from `Surv(time_year, event_main) ~ 1`.

### Plot Objects

- Plot `.png` files are image exports.
- Plot `.rds` files are saved plot objects for reuse in R.
  - Standard curve plots are `ggplot` objects.
  - Curve-plus-bar-count plots are `patchwork` objects.

## Table Column Definitions

### `simple_km_yearly_estimates.csv`

- `dataset`
  - Cohort name: `PNU`, `SNU`, or `merged`.
- `horizon_year`
  - Evaluation horizon in years.
- `estimated_survival_probability`
  - KM estimate of remaining transition-free survival by that horizon.
- `estimated_survival_lower_95`, `estimated_survival_upper_95`
  - Greenwood-based 95% confidence interval for the survival estimate where available.
- `estimated_risk_probability`
  - `1 - estimated_survival_probability`.
  - Interpreted as the cumulative transition risk on the main transition-only scale.
- `estimated_risk_lower_95`, `estimated_risk_upper_95`
  - Complement-transformed 95% confidence interval for the cumulative risk estimate.

### `simple_km_dataset_summary.csv`

- `dataset`
  - Cohort name.
- `n_subjects`
  - Number of rows/subjects in that cohort.
- `n_events_transition`
  - Number of subjects with `status_num == 1`.
- `n_events_remission`
  - Number of subjects with `status_num == 2`.
- `n_right_censored`
  - Number of subjects with `status_num == 0`.
- `max_followup_days`
  - Maximum observed follow-up in days.
- `max_followup_years`
  - Maximum observed follow-up in years.
  - For the current PNU data, this value is the observed endpoint used for the `observed_range` PNU plots.

### `simple_km_plot_manifest.csv`

- `scenario_key`
  - Scenario or file family identifier.
- `plot_type`
  - Plot subtype within that scenario.
- `png_file`
  - Absolute path to the exported PNG image.
- `rds_file`
  - Absolute path to the saved R plot object.

## Plot Families

### 1. Standard curve plots

These show the KM curve, annual point markers, and available 95% confidence interval guides.

Files:

- `simple_km_survival_plot_all_cohorts.*`
- `simple_km_risk_plot_all_cohorts.*`
- `simple_km_survival_plot_pnu_only.*`
- `simple_km_risk_plot_pnu_only.*`
- `simple_km_survival_plot_snu_only.*`
- `simple_km_risk_plot_snu_only.*`
- `simple_km_survival_plot_merged_only.*`
- `simple_km_risk_plot_merged_only.*`

Meaning:

- `survival`
  - Estimated transition-free survival probability.
- `risk`
  - Estimated cumulative transition risk, equal to `1 - survival`.
- `all_cohorts`
  - PNU, SNU, and merged on the same figure.
- `pnu_only`, `snu_only`, `merged_only`
  - Single-dataset versions.
- dotted step lines
  - Lower and upper 95% confidence interval guides for the KM curve.
- yearly interval bars
  - 95% confidence intervals at annual horizons where estimable.

Alias files:

- `simple_km_survival_plot.png`
  - Same content as the `all_cohorts` survival plot.
- `simple_km_risk_plot.png`
  - Same content as the `all_cohorts` risk plot.

### 2. PNU observed-range standard plots

These are PNU-only plots where the x-axis ends at the actual observed maximum PNU follow-up instead of extending to 10 years.

Files:

- `simple_km_survival_plot_pnu_observed_range.*`
- `simple_km_risk_plot_pnu_observed_range.*`

Meaning:

- Intended for interpretation that avoids visually emphasizing extrapolated time beyond the observed PNU follow-up window.
- The x-axis ends at `max_followup_years` from `simple_km_dataset_summary.csv`.

### 3. Curve-plus-bar-count plots

These show:

- The KM curve on the top panel.
- A bar panel for `Number at risk` below.
- A second bar panel for `Cumulative transitions` below that.

Files:

- `simple_km_survival_plot_pnu_with_counts.*`
- `simple_km_risk_plot_pnu_with_counts.*`
- `simple_km_survival_plot_snu_with_counts.*`
- `simple_km_risk_plot_snu_with_counts.*`
- `simple_km_survival_plot_merged_with_counts.*`
- `simple_km_risk_plot_merged_with_counts.*`

Interpretation:

- Bars are shown every 6 months.
- `Number at risk`
  - Number of subjects still under observation and event-free just before each displayed time point, operationalized here as subjects with `time_year >= t`.
- `Cumulative transitions`
  - Number of observed transition events accumulated up to and including time `t`, operationalized as `status_num == 1` and `time_year <= t`.
- Bar colors:
  - Blue bars = `Number at risk`
  - Orange bars = `Cumulative transitions`

### 4. PNU observed-range curve-plus-bar-count plots

These are the PNU bar-count plots restricted to the actual observed PNU follow-up range.

Files:

- `simple_km_survival_plot_pnu_with_counts_observed_range.*`
- `simple_km_risk_plot_pnu_with_counts_observed_range.*`

Meaning:

- Same information as the standard `pnu_with_counts` plots.
- The x-axis ends at the observed PNU maximum follow-up rather than at 10 years.

## PNU Vertical Line

Whenever PNU appears in a plot, a vertical dashed line is drawn at the observed PNU maximum follow-up.

Purpose:

- To separate the observed follow-up region from the visually extrapolated tail of the x-axis.

## How To Reuse The Saved Objects

### Fitted KM objects

Load:

```r
fit_objects <- readRDS("simple_km_fit_objects.rds")
fit_objects$PNU
fit_objects$SNU
fit_objects$merged
```

### Plot objects

Load:

```r
p <- readRDS("simple_km_survival_plot_pnu_with_counts.rds")
print(p)
```

## Practical Reading Guide

- Use `all_cohorts` plots for side-by-side comparison across datasets.
- Use `pnu_only` plots when the PNU curve is visually compressed by the longer SNU/merged follow-up.
- Use `pnu_observed_range` plots when you want the x-axis to stop at the actual PNU observation window.
- Use `with_counts` plots when you want to see how much information remains over time.
- Use `with_counts_observed_range` for the clearest PNU-specific interpretation.
