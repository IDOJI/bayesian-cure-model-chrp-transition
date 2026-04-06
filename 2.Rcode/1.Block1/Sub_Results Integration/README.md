# Results Integration README

## Active workflow

This folder now keeps only the active Block 1 results-integration scripts at the top level:

- `run_block1_results_integration_exports.r`
- `export_mixture_cure_fraction_tables.r`
- `export_simple_model_probability_comparison_tables.r`
- `integrate_simple_model_probability_plots.r`
- `plot_km_vs_mle_cure_integrated_curves.r`
- `plot_km_vs_mle_nocure_cure_integrated_curves.r`

The single recommended entry point is:

- `run_block1_results_integration_exports.r`

That wrapper does three things in order:

1. Archives the current contents of `/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis of CHR-P Using a Mixture Cure Model/1.Block1/5.Results Integration` into `legacy/<timestamp>/`
2. Regenerates the active cross-model outputs
3. Writes a fresh top-level `README.md` in the export folder

## Active export folders

The rebuilt output bundle contains these active folders:

- `mixture_cure_fraction_tables`
- `simple_model_probability_comparison_tables`
- `integrated_model_probability_plots`
- `km_vs_mle_cure_integrated_curves`
- `km_vs_mle_nocure_cure_integrated_curves`

These correspond to the required deliverables:

- cure fraction summaries and comparison tables
- probability comparison tables
- related comparison plots

## Legacy code

Scripts that were primarily intermediate or superseded are moved into `legacy/` under this code folder.

They are kept for traceability, but they are no longer part of the default export workflow.
