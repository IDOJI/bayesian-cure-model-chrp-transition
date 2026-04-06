# Block 5

`1.Cure-appropriateness screening.r` rebuilds Block 5 from the completed legacy Stage 6 outputs while enforcing the new 5-block specification. It is now treated as an archived legacy-reuse reference workflow rather than the published default.

`2.Raw-data recompute screening.r` is a self-contained recompute-oriented alternative built from the merged preprocessed CSV. It now aligns the published Block 5 headline screen with the main modeling family by using a prespecified `lognormal` non-cure versus cure comparison as the primary screen. The older multi-family RECeUS-AIC-style family comparison is retained as exploratory sensitivity material for candidate-map interpretation. It keeps Maller-Zhou contextual only and does not use a Shen proxy for Xie. If legacy Stage 6 Xie outputs are available, it can import those as the contradiction check; otherwise it leaves Xie as not implemented instead of using a proxy.

Main points:

- Published primary screening follows the prespecified `lognormal` main-model family.
- Exploratory family comparison keeps the multi-family `RECeUS-AIC` idea for sensitivity/candidate-map interpretation.
- `Xie` is used only as a contradiction-oriented follow-up check.
- `Maller-Zhou` is reported as an older contextual reference only.
- The default export root is `/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis of CHR-P Using a Mixture Cure Model/5.Block5`.
- The published export root keeps only the interpretation-facing files at top level:
  `block5_screening_summary.csv`, `block5_screening_overview.png`,
  `block5_receus_candidate_map.png`, and `block5_visual_summary.pdf`.
- The published export root also syncs the interpretation-supporting raw tables
  requested for close reading of the candidate map and AIC comparison:
  `block5_raw_cure_candidate_summary.csv`,
  `block5_raw_primary_model_fit_summary.csv`,
  `block5_raw_sensitivity_model_fit_summary.csv`,
  `block5_raw_receus_primary_summary.csv`,
  `block5_raw_maller_zhou_contextual_summary.csv`,
  `block5_raw_xie_status_summary.csv`, and
  `block5_raw_tail_summary.csv`.
- The raw recompute workflow now also writes plot-ready CSVs for every figure:
  `block5_raw_model_fit_plot_data.csv`,
  `block5_raw_receus_primary_map_plot_data.csv`,
  `block5_raw_screening_overview_plot_data.csv`, and
  `block5_raw_receus_candidate_map_plot_data.csv`.
- `block5_raw_screening_overview_plot_data.csv` includes the exact in-tile
  numeric labels used in the overview heatmap so the published PNG can be
  recreated without consulting intermediate objects.
- `block5_raw_receus_candidate_map_plot_data.csv` now includes the RECeUS threshold columns so the candidate-map PNG can be recreated from that single CSV without consulting code constants.
- `2.Raw-data recompute screening.r` writes its full audit bundle into
  `/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis of CHR-P Using a Mixture Cure Model/5.Block5/sub_raw_recompute`
  and syncs the published top-level files.
- `1.Cure-appropriateness screening.r` now writes by default into
  `/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis of CHR-P Using a Mixture Cure Model/5.Block5/archive_legacy_reuse`.
- Legacy supporting tables, bootstrap reuse files, and metadata are kept under the archive folder rather than the published root.
- Heavy bootstrap reruns are avoided. The script rebuilds `Xie` bootstrap summaries from the legacy Stage 6 bootstrap table when those records already exist.

Inputs used by default:

- `/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis of CHR-P Using a Mixture Cure Model/old/stage1_Backbone lock`
- `/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis of CHR-P Using a Mixture Cure Model/old/stage6_Cure-appropriateness screening`

Environment overrides:

- `CHR_COHORTS_DROPBOX_ROOT`
- `BLOCK5_STAGE1_DATA_PATH`
- `BLOCK5_LEGACY_STAGE6_EXPORT_PATH`
- `BLOCK5_EXPORT_PATH`
- `BLOCK5_RAW_MERGED_FILE`
- `BLOCK5_RAW_LEGACY_STAGE6_EXPORT_PATH`
- `BLOCK5_RAW_EXPORT_PATH`
- `BLOCK5_RAW_PUBLISHED_DIR`
