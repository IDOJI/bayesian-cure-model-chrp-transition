---
editor_options: 
  markdown: 
    wrap: 72
---

# Simplified Modeling Specification for a Transition-Only Cure Modeling Project

## 1. Scope and purpose

This document defines a simplified survival-modeling framework for the
current CHR-P project.

The goal is to keep the existing data dictionary and core backbone
intact while replacing the previous, more elaborate specification with a
lighter and more interpretable framework focused on **transition-only
survival modeling**. In this simplified project, remission is **not**
modeled as a separate outcome state in the main analysis. It is treated
as censoring.

Only the following four analyses are included:

1.  Kaplan–Meier survival estimation\
2.  Non-cure parametric survival model using **log-normal** only\
3.  Frequentist mixture cure model using **log-normal** latency only\
4.  Bayesian mixture cure model using the existing **Stage 8A-type
    transition-only structure** only

The Bayesian competing-risk remission-aware extension corresponding to
the former Stage 8B is excluded from this simplified project.

------------------------------------------------------------------------

## 2. Fixed backbone retained from the current project

The simplified project inherits the existing backbone without changing
variable definitions or dataset identity rules.

### 2.1 Analysis datasets

The following dataset views must be reported:

-   `SNU`
-   `PNU`
-   `merged`
-   `merged (site-adjusted)`

For `merged (site-adjusted)`, site is included as an adjustment
covariate in model-based analyses. Site effects are interpreted as
contextual proxy terms and **not** as clean treatment effects.

### 2.2 Unit of observation and identifiers

-   one row = one subject\
-   subject key = `site + id`

### 2.3 Time variables

-   main survival time variable: `days_followup`
-   reporting time scale: $$
    \text{time\_year} = \frac{\text{days\_followup}}{365.25}
    $$

`date_entry` is Day 0 and `days_followup` is elapsed time from
`date_entry`.

### 2.4 Core variables retained unchanged

-   `site`
-   `id`
-   `sex_num`
-   `age_exact_entry`
-   `days_followup`
-   `status_num`

### 2.5 Age scaling

For regression-based models, use the existing project scaling rule:

$$
\text{age\_s} = \frac{\text{age\_exact\_entry} - \operatorname{mean}(\text{age\_exact\_entry})}{2 \cdot \operatorname{sd}(\text{age\_exact\_entry})}
$$

### 2.6 Sex coding

-   `sex_num = 0` → Female
-   `sex_num = 1` → Male

------------------------------------------------------------------------

## 3. Outcome definition for the simplified project

This project uses a **single primary risk scale**.

### 3.1 Primary risk scale

**Label:** `transition_only_main`

-   event: `status_num == 1`
-   censoring: `status_num %in% c(0, 2)`

Therefore:

-   right censoring remains censoring
-   remission (`status_num == 2`) is treated as censoring

This rule applies to all four models in this simplified specification.

### 3.2 Interpretation consequence

All reported survival and cumulative risk values in the main tables and
main figures refer to:

-   **transition occurrence**
-   under a convention where remission is treated as censoring

This means that the project is estimating a transition-only event
process, not a remission-aware competing-risk process.

------------------------------------------------------------------------

## 4. Core covariate structure

To keep the framework simple and interpretable, the default covariate
set is:

-   `age_s`
-   `sex_num`

For `merged (site-adjusted)`, add:

-   `site`

No age-by-sex interaction is included in the default simplified project.

For cure models, the same default covariate set is used in both the
incidence and latency components unless explicitly justified otherwise
in a later implementation note.

------------------------------------------------------------------------

## 5. Model set

## 5.1 Model 1: Kaplan–Meier survival estimation

For each dataset view, estimate the overall transition-only survival
function:

$$
\hat S_{KM}(t)
$$

using:

$$
\text{Surv}(\text{days\_followup},\ I(\text{status\_num}=1))
$$

with remission treated as censoring.

The corresponding cumulative transition risk is:

$$
\hat R_{KM}(t) = 1 - \hat S_{KM}(t)
$$

This model is descriptive and unadjusted.

------------------------------------------------------------------------

## 5.2 Model 2: Non-cure parametric survival model (log-normal only)

Assume:

$$
\log T_i = \eta_i + \sigma \varepsilon_i,\qquad \varepsilon_i \sim N(0,1)
$$

with linear predictor:

$$
\eta_i = \beta_0 + \beta_a \, age_{s,i} + \beta_m \, sex_{num,i}
$$

and, for `merged (site-adjusted)`:

$$
\eta_i = \beta_0 + \beta_a \, age_{s,i} + \beta_m \, sex_{num,i} + \beta_s \, site_i
$$

The subject-specific survival and cumulative risk are:

$$
S_i(t)=1-\Phi\left(\frac{\log t-\eta_i}{\sigma}\right), \qquad
R_i(t)=1-S_i(t)
$$

This model assumes the whole population is susceptible.

### 5.2.1 Population-standardized reporting

Because Table 1 is defined on the **overall population** scale,
model-based estimates should be reported as covariate-standardized
averages over the empirical covariate distribution of the corresponding
analysis dataset:

$$
\bar S(t)=\frac{1}{n}\sum_{i=1}^n S_i(t), \qquad
\bar R(t)=1-\bar S(t)
$$

For `merged (site-adjusted)`, this standardization should preserve the
empirical site mix of the merged dataset unless a different
standardization target is explicitly declared.

------------------------------------------------------------------------

## 5.3 Model 3: Frequentist mixture cure model (log-normal latency only)

Let $U_i$ be the latent susceptibility indicator:

-   $U_i=1$: susceptible / uncured
-   $U_i=0$: non-susceptible / cured

### 5.3.1 Incidence component

Define the susceptible probability:

$$
\pi_i = P(U_i=1)
$$

with logistic incidence model:

$$
\operatorname{logit}(\pi_i)=\alpha_0+\alpha_a\,age_{s,i}+\alpha_m\,sex_{num,i}
$$

and, for `merged (site-adjusted)`:

$$
\operatorname{logit}(\pi_i)=\alpha_0+\alpha_a\,age_{s,i}+\alpha_m\,sex_{num,i}+\alpha_s\,site_i
$$

### 5.3.2 Latency component

For susceptible subjects:

$$
\log T_i \mid (U_i=1)=\eta_i+\sigma\varepsilon_i,\qquad \varepsilon_i\sim N(0,1)
$$

with:

$$
\eta_i=\beta_0+\beta_a\,age_{s,i}+\beta_m\,sex_{num,i}
$$

and, for `merged (site-adjusted)`:

$$
\eta_i=\beta_0+\beta_a\,age_{s,i}+\beta_m\,sex_{num,i}+\beta_s\,site_i
$$

The susceptible-only survival is:

$$
S_{u,i}(t)=1-\Phi\left(\frac{\log t-\eta_i}{\sigma}\right)
$$

### 5.3.3 Overall survival and risk

The overall survival is:

$$
S_i(t)=(1-\pi_i)+\pi_i S_{u,i}(t)
$$

and the overall cumulative transition risk is:

$$
R_i(t)=1-S_i(t)=\pi_i\{1-S_{u,i}(t)\}
$$

### 5.3.4 Susceptible-only quantities

The latent susceptible-only quantities are:

$$
S_{u,i}(t), \qquad R_{u,i}(t)=1-S_{u,i}(t)
$$

For dataset-level reporting, use covariate-standardized averages:

$$
\bar S_u(t)=\frac{1}{n}\sum_{i=1}^n S_{u,i}(t), \qquad
\bar R_u(t)=1-\bar S_u(t)
$$

These are the quantities reported in Table 2.

In the current simplified implementation, these susceptible-only
cumulative risk probabilities are explicitly exported on the 1-10 year
grid and can also appear in supplementary comparison CSV files and plots
when they are clearly labeled as `susceptible-only`.

------------------------------------------------------------------------

## 5.4 Model 4: Bayesian mixture cure model (transition-only Stage 8A-type structure only)

The Bayesian model keeps the same transition-only cure decomposition as
Model 3 but estimates it in a Bayesian framework.

### 5.4.1 Scope restriction

This simplified Bayesian block includes only the **main Stage 8A-type
transition-only cure structure**.

The following are excluded:

-   remission-aware competing-risk modeling
-   Stage 8B-type remission-sensitive branch
-   sensitivity analyses
-   prior-family sensitivity branches
-   alternative latency family competition

### 5.4.2 Likelihood

Let:

$$
\delta_i = I(\text{status\_num}_i = 1)
$$

Then the observed-data likelihood contribution is:

$$
L_i = [\pi_i f_u(T_i)]^{\delta_i}\left[(1-\pi_i)+\pi_i S_u(T_i)\right]^{1-\delta_i}
$$

This means:

-   if a subject transitions, they must belong to the susceptible
    subgroup
-   if a subject is censored, they may be either non-susceptible or
    susceptible but not yet transitioned

### 5.4.3 Bayesian targets

The Bayesian branch must report:

-   posterior overall survival $S(t)$
-   posterior overall cumulative risk $R(t)$
-   posterior cure fraction / susceptible fraction
-   posterior susceptible-only survival $S_u(t)$
-   posterior susceptible-only cumulative risk $R_u(t)$

all on the 1–10 year grid.

In plain-language reporting, this means the simplified Bayesian analysis
also includes risk probabilities evaluated for the latent susceptible
subgroup alone, not only the overall-population risk.

------------------------------------------------------------------------

## 6. Tables and integrated exports to be produced

## 6.1 Table 1. Overall survival and risk comparison

For each dataset view and each horizon:

$$
t = 1,2,\ldots,10\ \text{years}
$$

report the following overall-population quantities:

-   Kaplan–Meier estimated survival probability
-   log-normal non-cure estimated survival probability
-   frequentist log-normal mixture cure estimated survival probability
-   Bayesian log-normal mixture cure posterior mean survival probability
-   corresponding cumulative risks $1-S(t)$

If interval estimates or horizon-support governance fields are available
in a harmonized export, they may be appended. In the current simplified
implementation, the core required output is the point-estimate
comparison on the shared 1-10 year grid.

This table must be on the **overall population scale** only.

------------------------------------------------------------------------

## 6.2 Table 2. Susceptible-only quantities from mixture cure models

For each dataset view and each horizon:

$$
t = 1,2,\ldots,10\ \text{years}
$$

report:

-   frequentist mixture cure susceptible-only survival
-   frequentist mixture cure susceptible-only cumulative risk
-   Bayesian mixture cure susceptible-only survival
-   Bayesian mixture cure susceptible-only cumulative risk

These susceptible-only cumulative risk probabilities are part of the
current simplified analysis and are explicitly reported for the two cure
models.

### 6.2.1 Required interpretation rule

Primary tabular reporting should keep Table 1 and Table 2 conceptually
separated.

-   Table 1 = overall population
-   Table 2 = latent susceptible-only subgroup

However, supplementary comparison exports may place these quantities
side by side, provided that:

-   susceptible-only entries are explicitly labeled as
    `susceptible-only`
-   KM and non-cure entries are explicitly described as
    overall/reference values
-   the figure subtitle or caption explains that more than one
    estimation scope is being shown

## 6.3 Integrated horizon-wise risk comparison exports

For each horizon $t=1,\ldots,10$, the simplified workflow may also
export a comparison CSV designed for quick cross-model review.

Each row contains:

-   dataset name
-   model name
-   estimated risk probability

In the current implementation, these integrated comparison CSV files
include:

-   overall risk values for Kaplan-Meier, non-cure, frequentist cure,
    and Bayesian cure models
-   susceptible-only risk values for the frequentist and Bayesian cure
    models
-   a pooled KM reference duplicated into `merged (site-adjusted)`
    comparisons for visual and tabular benchmarking

These integrated CSV exports are convenience comparison files. They do
not erase the conceptual distinction between overall-population and
susceptible-only risk; that distinction must remain explicit in the
model labels.

------------------------------------------------------------------------

## 7. Plot sets and comparison visualizations to be produced

## 7.1 Plot Set 1. Overall survival and cumulative risk curves

For each high-level dataset view, produce:

1.  Overall survival curves\
2.  Overall cumulative risk curves

In the current implementation, the high-level dataset views are:

-   `PNU`
-   `SNU`
-   `Merged`

The `Merged` view is a comparison view that overlays both non-adjusted
and site-adjusted model variants in the same figure.

Overlay the following on the same axis:

-   Kaplan–Meier
-   log-normal non-cure
-   frequentist log-normal mixture cure
-   Bayesian log-normal mixture cure
-   for `Merged`, the corresponding non-adjusted and site-adjusted model
    variants

### 7.1.1 Required visual rules

-   x-axis: years since cohort entry
-   integer-year marks from 1 through 10
-   survival and risk must be clearly separated by title, y-axis, and
    legend
-   all dataset views must use a consistent graphical template
-   late horizons should be discussed in captions or surrounding text as
    increasingly model-dependent, even if a dedicated shaded projection
    band is not drawn

### 7.1.2 Special rule for `merged (site-adjusted)`

Because Kaplan–Meier itself does not adjust for site, the
`merged (site-adjusted)` comparison reuses the pooled KM curve from the
unadjusted merged dataset as a benchmark reference.

In a compact integrated plot, this may appear simply as `KM`, but it
should be interpreted as:

-   **pooled KM reference (unadjusted)**

This preserves comparability while preventing misinterpretation.

------------------------------------------------------------------------

## 7.2 Plot Set 2. Susceptible-only survival and cumulative risk curves

For each dataset view, produce:

1.  Susceptible-only survival curves\
2.  Susceptible-only cumulative risk curves

The current simplified implementation also allows comparison-oriented
figures in which KM and non-cure curves appear together with the
cure-model susceptible-only curves as reference comparators.

The figure titles, subtitles, legends, or captions must explicitly
include:

-   `susceptible-only`
-   or `uncured-only`
-   or an equally clear explanation of which rows/curves are true
    susceptible-only estimates and which rows/curves are
    overall/reference values

### 7.2.1 Curves to include

Because the user explicitly requested the following visual comparison,
Plot Set 2 should include:

-   Kaplan–Meier
-   frequentist non-cure log-normal model
-   frequentist mixture cure model
-   Bayesian mixture cure model

### 7.2.2 Interpretation safeguard for Plot Set 2

Only the frequentist and Bayesian cure models provide true latent
susceptible-only estimates.

Therefore:

-   the non-cure log-normal curve should be labeled as an
    **all-susceptible reference**
-   the KM curve should be labeled as an **overall empirical reference,
    not a latent susceptible-only estimator**

This preserves the requested visual comparison while preventing category
confusion.

### 7.2.3 Additional visual rules

-   x-axis: years since cohort entry
-   integer-year marks from 1 through 10
-   captions or subtitles must explain when overall/reference values and
    susceptible-only values appear together in one comparison figure
-   this plot set should be presented in a section visually separated
    from overall-population plots when possible

## 7.3 Integrated risk heatmaps and horizon-wise comparison plots

In addition to the curve plots, the current simplified analysis may
export:

-   dataset-level risk heatmaps with annotated risk probabilities across
    the 1-10 year grid
-   horizon-specific bar plots for cross-model risk comparison
-   combined comparison figures where overall/reference risks and
    susceptible-only risks appear in separate panels or clearly labeled
    entries

These integrated comparison visuals are supplementary tools for quick
model review. They may place more than one estimation scope in the same
figure, but the scope of each row, bar, or curve must be made explicit
in the labels or caption.

------------------------------------------------------------------------

## 8. Site-adjustment rule

In this simplified project, `site-adjusted` means:

-   in non-cure log-normal models: add `site` to the linear predictor
-   in frequentist mixture cure models: add `site` to both incidence and
    latency components
-   in Bayesian mixture cure models: add `site` to both incidence and
    latency components

No site interaction terms are included in the default simplified
specification.

Any site coefficient must be interpreted as a contextual structural
proxy and not as a treatment-effect parameter.

------------------------------------------------------------------------

## 9. Horizon-support and projection governance

This simplified project keeps the existing support logic.

### 9.1 PNU

-   1 year: `primary_supported`
-   2 years: `sensitivity`
-   3+ years: `projection`

### 9.2 SNU

-   1–2 years: `primary_supported`
-   3–5 years: `secondary`
-   6+ years: `projection`

### 9.3 merged

-   1–2 years: `primary_supported`
-   3–5 years: `secondary`
-   6+ years: `projection`

### 9.4 Horizon evidence class

When a harmonized horizon registry is attached to an output, each
horizon row should also carry:

-   `directly_observed_data_supported`
-   `partly_model_dependent`
-   `mostly_extrapolated`

Preferred mapping:

-   PNU:
    -   1 year → directly observed
    -   2 years → partly model-dependent
    -   3+ years → mostly extrapolated
-   SNU / merged:
    -   1–2 years → directly observed
    -   3–5 years → partly model-dependent
    -   6+ years → mostly extrapolated

In the current simplified implementation, these governance labels are
not yet attached uniformly to every stage-specific CSV export. They
remain reporting rules that can be merged into harmonized comparison
outputs.

### 9.5 Claim restriction

When horizon-governance metadata are attached, each horizon row should
also carry:

-   `primary_claim_allowed`
-   `secondary_or_sensitivity_only`
-   `projection_only`

Late-horizon results must be discussed as projection-dependent and not
as primary empirical claims.

------------------------------------------------------------------------

## 10. Reporting and interpretation rules

## 10.1 Main comparison story

The simplified reporting should answer four linked questions:

1.  How do overall survival and cumulative risk differ across the four
    models?
2.  Does allowing a cure component materially change supported-horizon
    risk estimates?
3.  What latent susceptible-only trajectory is implied by the cure
    models?
4.  How much does the site-adjusted merged analysis differ from the
    site-unadjusted merged analysis?

## 10.2 Primary interpretation window

Primary interpretation should emphasize:

-   `PNU`: 1 year
-   `SNU`: 1–2 years
-   `merged`: 1–2 years

### 10.3 Projection language

For later horizons, wording must explicitly state whether the estimate
is:

-   directly follow-up supported
-   partly model dependent
-   mostly extrapolated

### 10.4 What not to claim

This simplified project must not claim that:

-   the analysis proves that a cured subgroup truly exists
-   remission has no impact, merely because it was censored in the
    simplified design
-   site effect estimates represent treatment effects
-   late-horizon extrapolation is equivalent to direct observation

------------------------------------------------------------------------

## 11. Practical implementation note

This specification is designed to be used immediately in a
“transition-only, remission-as-censoring” workflow. It is intentionally
shorter and lighter than the previous full framework, but it preserves:

-   the existing backbone and variable dictionary
-   the shared time horizon grid
-   the main comparative model logic
-   the distinction between overall-population and susceptible-only
    quantities
-   the current susceptible-only risk outputs from the two cure models
-   the supplementary integrated comparison CSV, bar-plot, and heatmap
    workflow
-   and the supported-horizon interpretation discipline

If a later analysis needs explicit remission competition, that should be
treated as a separate extension rather than being silently folded into
this simplified main specification.

------------------------------------------------------------------------

## 12. Drafting basis for this simplified specification

This simplified document was written to preserve the current project
backbone while simplifying the master framework and the Bayesian Stage 8
design into a transition-only main analysis. It was drafted against the
existing project backbone, data dictionary, and Bayesian Stage 8A
structure, rather than inventing a new schema from scratch.
