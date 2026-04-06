# Block4 Metric Provenance Note

This note is intended to accompany the Block4 interpretation outputs for reviewer tracing.

## Code package scope

- Main Block4 script:
  `/Users/ido/Documents/GitHub/bayesian-cure-model-chrp-transition/2.Rcode/4.Block4/1.Follow-up adequacy and tail stability assessment.r`
- There is currently no separate helper-function file for Block4.
- The helper functions used for the requested metrics are defined inline inside the same main Block4 script.

## Requested metrics

### 1. `latest_stable_horizon_year`

- Computed in `build_block4_tail_window_summary()` at lines 2661-2688.
- Exact rule at line 2675:
  `max(horizon_year[completeness_plot_status == 'stable'], na.rm = TRUE)`
- Infinite results are converted to `NA` at line 2685.
- This metric is exported in:
  - `sub_supporting_tables/block4_tail_window_summary.csv`
  - `block4_late_tail_narrative_summary.csv`

### 2. `first_flagged_horizon_year`

- Computed in `build_block4_tail_window_summary()` at lines 2661-2688.
- Exact rule at lines 2677-2680:
  `min(horizon_year[completeness_plot_status %in% c('low_fraction', 'sparse', 'very_sparse', 'masked_no_eligible', 'masked_no_risk_set')], na.rm = TRUE)`
- Infinite results are converted to `NA` at line 2687.
- This metric is exported in:
  - `sub_supporting_tables/block4_tail_window_summary.csv`
  - `block4_late_tail_narrative_summary.csv`

### 3. `normalized_auc_*`

- The underlying per-time quantity is `difference_width`, defined in `build_block4_betensky_curve_data()` at line 2655 as:
  `pmax(as.numeric(conf_high) - as.numeric(conf_low), 0)`
- Trapezoidal integration is implemented in `compute_trapezoid_area()` at lines 296-313.
- Normalization is implemented in `compute_normalized_difference_area()` at lines 316-329.
- The normalization divisor is the observed time span:
  `max(time_year) - min(time_year)`
- The summary metrics are assigned in `summarise_block4_betensky_metrics()` at lines 2692-2747.
- Exact output fields:
  - `normalized_auc_full` at line 2741
  - `normalized_auc_late_tail` at line 2744
  - `normalized_auc_0_2` at line 2745
  - `normalized_auc_2_5` at line 2746
  - `normalized_auc_5_10` at line 2747
- Time-window rules used for these summaries:
  - Full curve: all available `time_year` values for the selected group
  - Late tail: `time_year >= late_tail_start_year`, where `late_tail_start_year = coalesce(primary_supported_max_year, 1)` at line 2688
  - `0-2y`: `time_year <= 2` at lines 2714-2717
  - `2-5y`: `time_year > 2 & time_year <= 5` at lines 2718-2721
  - `5-10y`: `time_year > 5` at lines 2722-2725
- These metrics are exported in:
  - `block4_betensky_summary.csv`
  - `block4_late_tail_narrative_summary.csv`

### 4. `site_dominance_*`

- The merged-site contribution profile is built in `compute_site_contribution_profile()` at lines 1791 onward.
- This logic is only applied for the merged dataset and only for `overall` or `sex` subgroup summaries.
- Configuration values are defined near the top of the script:
  - `site_dominance_warning_fraction <- 0.80` at line 69
  - `site_dominance_warning_min_horizon_year <- 3L` at line 70
  - `site_dominance_warning_exempt_primary_supported <- TRUE` at line 71
- Core status classification is assigned at lines 1839-1844:
  - `no_eligible_subjects` if `total_eligible_n <= 0`
  - `single_site_only` if `eligible_site_count <= 1`
  - `single_site_dominant` if `dominant_site_fraction >= 0.80`
  - otherwise `multi_site_balanced`
- Descriptive note text for the raw status is assigned at lines 1846-1850.
- Warning escalation logic is separated from the raw status and implemented in `should_raise_site_dominance_warning()` at lines 1752-1766.
- Exact warning behavior:
  - `single_site_only` always raises a warning
  - `single_site_dominant` raises a warning only when:
    - `horizon_year >= 3`, and
    - the support tier is not exempted primary-supported support
  - `not_applicable`, `not_computed`, `multi_site_balanced`, and `no_eligible_subjects` do not raise warnings
- Reviewer-facing explanatory text for the warning layer is generated in `make_site_dominance_warning_note()` at lines 1769-1788.
- These variables are exported in:
  - `stage2_followup_horizon_summary.csv`
  - `sub_supporting_tables/stage2_merged_site_contribution.csv`
- Relevant exported fields include:
  - `dominant_site`
  - `dominant_site_n`
  - `dominant_site_fraction`
  - `site_dominance_status`
  - `site_dominance_note`
  - `site_dominance_warning_flag`
  - `site_dominance_warning_note`

## Narrative-layer reuse

- `block4_late_tail_narrative_summary.csv` re-derives the detailed horizon fields from `horizon_summary` at lines 2804-2818.
- It then coalesces those detailed fields back into the exported narrative summary at lines 2836-2838:
  - `latest_stable_horizon_year = coalesce(latest_stable_horizon_year_detail, latest_stable_horizon_year)`
  - `first_instability_horizon_year = coalesce(first_instability_horizon_year_detail, first_flagged_horizon_year)`

## Reviewer package recommendation

For line-trace review, send these together:

1. `1.Follow-up adequacy and tail stability assessment.r`
2. `Block4_metric_provenance.md`
3. `stage2_followup_horizon_summary.csv`
4. `sub_supporting_tables/block4_tail_window_summary.csv`
5. `block4_betensky_summary.csv`
6. `block4_late_tail_narrative_summary.csv`
7. `sub_supporting_tables/stage2_merged_site_contribution.csv`
