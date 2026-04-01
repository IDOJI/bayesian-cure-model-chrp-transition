# Stage 11 Simplified Transition-Only Cure Project

This sub-project is a lightweight reporting layer for a simplified transition-only survival modeling project.

It reuses already-generated outputs from:

- `stage3_KM benchmark classifier`
- `stage5_Individualized no-cure comparator`
- `stage7_Frequentist cure block`
- `stage8A_Bayesian transition-only cure`

The Stage 11 script does not refit those upstream models. Instead, it:

- restricts the comparison to four models only,
- treats remission as censoring throughout,
- standardizes reporting to four dataset versions,
- exports simplified tables and plots,
- and carries forward the Stage 8A Bayesian diagnostic PDF as the trace-plot reference file.

## Dataset versions

- `PNU`
- `SNU`
- `merged`
- `merged (site-adjusted)`

## Model set

- `Kaplan-Meier`
- `Non-cure log-normal`
- `Frequentist mixture cure (log-normal)`
- `Bayesian mixture cure (log-normal, anchor-informed)`

## Main run command

```bash
'/Library/Frameworks/R.framework/Resources/bin/Rscript' \
  '2.Rcode/stage11_simplified_transition_only_cure_project/run_stage11_simplified_transition_only_cure_project.R'
```

## Optional export override

Use `STAGE11_EXPORT_PATH` when testing outside Dropbox.

```bash
STAGE11_EXPORT_PATH=/tmp/stage11_demo \
'/Library/Frameworks/R.framework/Resources/bin/Rscript' \
  '2.Rcode/stage11_simplified_transition_only_cure_project/run_stage11_simplified_transition_only_cure_project.R'
```
