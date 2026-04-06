# Block 6

`1.Cure decomposition stability and latency fragility.r` is the Block 6 entry script.

Block 6 standardizes four interpretation-facing diagnostics that qualify cure-fraction and uncured-latency claims:

- `leave-one-out / leave-last-k-out` influence analysis
- `latency family sensitivity`
- `Bayesian joint uncertainty` between cure fraction and latency scale
- `supported-horizon overlay` against Block 4 primary-supported and latest-stable follow-up windows

Default Dropbox export root:

- `/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/Survival Analysis of CHR-P Using a Mixture Cure Model/6.Block6`

Export layout:

- Top-level files: interpretation-facing summary tables, plots, `block6_supported_horizon_overlay.csv`, and `README.md`
- `2.Supporting Tables/`: long-form influence diagnostics, baseline validation, and Bayesian draw exports
- `3.Technical/`: export manifest and code references

Export notes:

- top-level plots are rebuilt from exported CSV files
- `block6_file_guide.md` describes the generated files one by one

Implementation note:

- The entry script reuses the core diagnostics engine at `2.Rcode/1.Block1/Main/5.Latency instability priority diagnostics.r` and syncs the outputs into the Block 6 export structure.
