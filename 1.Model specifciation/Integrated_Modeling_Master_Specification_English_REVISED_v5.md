# Integrated Modeling Master Specification (English, Revised v5)

**Status:** revised consolidated master specification for the current project (feedback-integrated v5)  
**Scope:** shared backbone, staged workflow, frequentist/non-Bayesian operational rules, Stage 6 screening specification, cross-stage output contracts, interpretation guardrails, and stage-specific reporting/QC rules  
**Separate companion retained:** `Bayesian_Modeling_Specification_Stage8_REVISED_v5.md`

---

## 0. Document role and intended use

This file is the **single integrated non-Bayesian/shared master specification** for the current project.

It is intended to **fully absorb** the project content previously spread across:

- `1.Analysis Objective_🇬🇧ENG.md`
- `2.General Model Specifications_🇬🇧ENG.md`
- `3.Advice from previous modeling results_🇬🇧ENG.md`
- `4.Modeling Framework Breakdown_🇬🇧ENG.md`
- `5.💛Model Specification Framework_🇬🇧ENG.md`
- `staged_workflow_guideline.md`

and to absorb the remaining overlapping operational material from:

- `🟩전체 모델링 검토.txt`
- the non-Bayesian parts of later planning notes
- the Stage 6 frequentist-only screening specification
- the stage-specific reporting/QC prompt, adapted into English reporting rules

Accordingly, this file now serves **simultaneously** as:

1. the **canonical shared framework** for the staged workflow,
2. the **operational stage-by-stage breakdown** for what each stage computes and exports,
3. and the **workflow execution guide** describing round structure, parallelization limits, and carry-forward logic.

This file is **not** intended to replace the separate Bayesian companion document:

- `Bayesian_Modeling_Specification_Stage8_REVISED_v5.md`

That Bayesian document remains the dedicated Stage 8 specification for full Bayesian modeling detail, prior encoding, admissibility, and exported Bayesian outputs.

### 0.1 What this file now replaces

For project interpretation and non-Bayesian/shared modeling rules, this file now replaces the need to keep the following documents as separate governing workflow/specification sources:

- `1.Analysis Objective_🇬🇧ENG.md`
- `2.General Model Specifications_🇬🇧ENG.md`
- `3.Advice from previous modeling results_🇬🇧ENG.md`
- `4.Modeling Framework Breakdown_🇬🇧ENG.md`
- `5.💛Model Specification Framework_🇬🇧ENG.md`
- `staged_workflow_guideline.md`

Their roles are now integrated below in:

- the scientific-target section,
- the canonical common backbone,
- the shared evaluation/output rules,
- the stage-specific operational sections,
- the workflow/parallelization section,
- and the implementation-facing notes.

### 0.2 What remains intentionally separate

The following files should remain separate companion documents:

- `Rules Before Generating R Code_🇬🇧ENG.md` for code-generation behavior
- `3.Data Dictionary_🇬🇧ENG.md` for variable definitions and field semantics
- `Bayesian_Modeling_Specification_Stage8_REVISED_v5.md` for Stage 8 Bayesian detail

### 0.3 Backward-compatibility rule

Because the existing codebase was already written against the current backbone, this integrated file intentionally preserves that backbone unless revision is unavoidable because of a true contradiction.

When a new preferred naming convention is introduced, legacy names may be retained temporarily for compatibility.

### 0.4 Legacy-reference rule

Older scripts, notes, or manifests may still cite the absorbed filenames listed above.

Those legacy citations should now be interpreted as pointing to the corresponding sections of **this integrated master specification**.


## 1. Governing-source hierarchy

Use the following priority order when interpreting or extending the workflow:

1. **this integrated master specification** for all shared and non-Bayesian project rules
2. `Rules Before Generating R Code_🇬🇧ENG.md` for code-generation behavior
3. `3.Data Dictionary_🇬🇧ENG.md` for variable definitions and field semantics
4. `Bayesian_Modeling_Specification_Stage8_REVISED_v5.md` for Stage 8 only
5. stage-specific code, exported metadata, and result registries

### 1.1 Conflict rule

If an older overlapping note conflicts with this integrated master file, follow this integrated file unless a later stage-specific exported artifact explicitly documents a justified exception.

### 1.2 Language rule

This integrated master specification is written in **English**.

A later stage report may be written in another language without changing the scientific or modeling rules recorded here.

### 1.3 Deletion/archiving rule for absorbed documents

Once this file has been adopted as the active non-Bayesian/shared master specification, the absorbed workflow/specification files

- `1.Analysis Objective_🇬🇧ENG.md`
- `2.General Model Specifications_🇬🇧ENG.md`
- `3.Advice from previous modeling results_🇬🇧ENG.md`
- `4.Modeling Framework Breakdown_🇬🇧ENG.md`
- `5.💛Model Specification Framework_🇬🇧ENG.md`
- `staged_workflow_guideline.md`

may be archived or removed **without changing the intended shared non-Bayesian modeling logic**, because their governing content is now integrated here.

This rule does **not** apply to:

- `Rules Before Generating R Code_🇬🇧ENG.md`
- `3.Data Dictionary_🇬🇧ENG.md`
- `Bayesian_Modeling_Specification_Stage8_REVISED_v5.md`



## 2. Core scientific target

The central scientific question is:

> **In future CHR-P research, does allowing cure-like heterogeneity materially change estimated transition risk and reduce unnecessary high-risk classification compared with KM and no-cure models, and does externally informed Bayesian cure modeling add practical value beyond a no-external-information Bayesian fit?**

This project does **not** aim to prove that a cured subgroup exists as an observed fact.

Its main aim is to show whether, relative to KM and no-cure formulations, models that allow a cure-like heterogeneous subgroup change clinically relevant risk estimates enough to matter, and whether Bayesian modeling informed by large external age- and sex-specific incidence evidence improves stabilization, interpretability, and false-positive burden.

The intended top-line message is therefore conditional and comparative:

> **Future CHR-P studies may need modeling frameworks that allow cure-like heterogeneity; when such models are used, Bayesian estimation informed by large external meta-analytic incidence information may be preferable if it materially changes supported-horizon risk estimates or reduces unnecessary high-risk classification.**

Even a modest change can be clinically meaningful. In this project, a 3 percentage-point reduction in unnecessary high-risk classification should be read as **3 fewer false positives per 100 people assessed**, which is clinically relevant at the individual-patient level.

### 2.1 Five integrated project objectives

The project has five linked objectives.

#### 2.1.1 Cohort-level objective

For each dataset:

- `PNU`
- `SNU`
- `merged`

evaluate how the estimated cumulative transition risk changes when cure-like heterogeneity is allowed.

The main quantities of interest are the estimated transition risks at:

- 1 year,
- 2 years,
- ...,
- 10 years after cohort entry

and these should be compared across:

- KM-based risk estimates
- no-cure models
- frequentist mixture cure models
- Bayesian mixture cure models

For Bayesian cure models, the cohort-level objective also includes an explicit **external-information value assessment**:

- compare anchor-informed versus neutral / no-external-information Bayesian fits
- quantify how much risk, false-positive burden, and net benefit change when the external meta-analytic incidence anchor is removed
- show how the anchored incidence pattern is updated by the posterior rather than treating the anchor as fixed truth

#### 2.1.2 Individual-level objective

At pre-specified prediction horizons and risk thresholds, compare:

- KM-based benchmark classifiers
- no-cure models
- cure models

with respect to:

- unnecessary high-risk classification
- false-positive burden
- false positives per 100
- related operating characteristics
- and clinical usefulness

KM must therefore be treated not only as a descriptive curve but also as a **formal horizon-specific benchmark classifier**.

The project is interested in whether cure-aware modeling changes downstream patient labeling enough to matter. Even small absolute differences should therefore be translated into natural-frequency terms whenever possible.

#### 2.1.3 Cohort-comparison objective

Using `PNU`, `SNU`, and `merged`, compare:

- how the estimated risk trajectory changes across cohorts
- how false-positive burden changes across cohorts
- how follow-up maturity and shrinking risk sets alter interpretation
- how external incidence anchoring does or does not stabilize the Bayesian fits
- and how these patterns affect apparent evidence for or against cure-like heterogeneity

#### 2.1.4 Within-CHR-P heterogeneity characterization objective

The project also aims to characterize clinically meaningful heterogeneity within the CHR-P population itself.

This objective is intentionally indirect and triangulated, not class-identifying. Evidence should be read jointly from:

- timing-separation analysis,
- follow-up sufficiency and contradiction checks,
- cure-appropriateness screening,
- latency-family sensitivity,
- remission-aware comparison,
- anchor-informed versus no-external-information Bayesian comparison,
- same-family cure-versus-non-cure fit-gain diagnostics,
- and the cure-model-only supporting decomposition block.

The goal is to evaluate whether a single homogeneous no-cure survival process appears insufficient for the observed CHR-P data.

This objective does not authorize claims that the current data prove a discrete cured subgroup or identify the true number of latent classes.

#### 2.1.5 Future-method recommendation objective

A linked methodological objective is to determine whether future CHR-P survival analyses should consider mixture-cure-type models when supported-horizon risk shifts, false-positive burden changes, and triangulated heterogeneity diagnostics indicate that homogeneous no-cure models are clinically incomplete.

### 2.2 Scientific priorities

The project therefore prioritizes seven linked aims:

1. **cohort-level cumulative transition-risk change**
2. **horizon-specific false-positive burden**
3. **clinical usefulness on the same horizon/threshold grid**
4. **explicit separation of timing difference, follow-up immaturity, and cure-like heterogeneity**
5. **explicit comparison of externally anchored versus no-external-information Bayesian cure modeling**
6. **triangulated indirect assessment of within-CHR-P heterogeneity**
7. **a future-method recommendation about whether mixture-cure-type modeling should be considered in later CHR-P survival studies**

This project is **not** primarily asking which model has the lowest AIC or the highest global fit score.

It is also **not** primarily asking for proof that a cure fraction truly exists.

The primary scientific emphasis is:

- horizon-specific prediction after cohort entry
- threshold-based unnecessary high-risk labeling
- whether allowing cure-like heterogeneity changes those decisions enough to matter clinically
- whether the external meta-analytic incidence anchor materially changes supported-horizon conclusions
- whether the total triangulated evidence is more compatible with clinically meaningful within-CHR-P heterogeneity than with a single homogeneous no-cure survival process
- whether current evidence is strong enough to recommend that future CHR-P survival studies explicitly consider mixture-cure-type models
- and whether any apparent Bayesian stabilization is data-supported rather than merely prior-driven

### 2.3 Cohort context that motivates the staged workflow

The current empirical context is:

- `PNU` has about 2 years of follow-up and shows a consistently earlier-event pattern
- `SNU` includes some patients followed up to about 10 years
- in `SNU`, the risk set becomes progressively smaller after about 5 years, so late-horizon estimation becomes increasingly unstable

A key principle follows from this:

- late-horizon instability must be clearly quantified and reported as a limitation
- but it must **not** automatically replace the main scientific question, which remains the comparison of horizon-specific prediction and false-positive burden from 1 to 10 years after cohort entry

### 2.4 What the project ultimately wants to see

The final integrated analysis should show:

- how much risk estimation changes when a cure structure is allowed
- whether that change reduces false-positive burden at clinically relevant horizons and thresholds
- whether externally informed Bayesian modeling changes supported-horizon risk, false-positive burden, or net benefit compared with a neutral / no-external-information Bayesian fit
- how the incidence anchor shape changes from prior center to posterior after estimation
- how that pattern differs according to follow-up maturity, cohort, and site structure
- whether observed differences reflect true timing difference, limited follow-up, remission handling, external-prior dependence, or evidence for a cure fraction
- whether the total pattern of evidence is more compatible with clinically meaningful within-CHR-P heterogeneity than with a single homogeneous no-cure survival process
- whether the project can responsibly recommend that future CHR-P survival studies consider mixture-cure-type models
- and how much confidence should be placed in late-horizon comparisons given shrinking risk sets and tail instability

## 3. Core interpretation rules that apply across all stages

### 3.1 Three phenomena must never be blurred

The following are analytically distinct and must be reported separately:

1. **timing difference**
2. **follow-up immaturity / late-tail instability**
3. **cure-fraction or cure-like heterogeneity evidence**

A late plateau, tail separation, or long-horizon divergence is **not** by itself evidence of cure.

### 3.1A Triangulation rule for within-CHR-P heterogeneity

No single stage, model family, or diagnostic may by itself establish within-CHR-P heterogeneity.

The heterogeneity argument must instead be triangulated from multiple distinct sources:

- Stage 4 timing separation,
- Stage 6 cure-appropriateness screening and follow-up contradiction checks,
- Stage 7 cure-versus-non-cure risk change and same-family fit-gain diagnostics,
- Stage 8 anchor-informed versus no-external-information Bayesian comparison,
- Stage 9 remission-aware reanalysis,
- and Stage 10 integrated interpretation.

Accordingly, any apparent cure-like signal must be read as indirect evidence whose strength depends on the coherence of these sources rather than on a single statistically favorable result.

### 3.2 Site-effect interpretation rule

There is currently no direct information on how drug treatment was administered.

Therefore, any site effect must be interpreted as a proxy for:

- broader treatment context,
- care pathway,
- referral structure,
- selection,
- or follow-up maturity,

and **not** as a clean treatment effect.

This rule applies to:

- Stage 4 site contrasts,
- Stage 5 site-adjusted non-cure models,
- Stage 7 incidence/latency site placement,
- Stage 8 incidence/latency/remission site coefficients,
- and any merged overall or sex summary that carries a site term.

Accordingly, no merged site coefficient may be written up as a causal clinical treatment effect. It is a structural proxy for contextual differences unless genuine treatment-process variables are later added.

### 3.3 Supported-horizon hierarchy

Interpretation priority is fixed as follows:

- **PNU**
  - 1 year = primary-supported
  - 2 years = sensitivity
  - 3+ years = projection

- **SNU**
  - 1–2 years = primary-supported
  - 3–5 years = secondary
  - 6+ years = projection

- **merged**
  - 1–2 years = primary-supported
  - 3–5 years = secondary
  - 6+ years = projection

Late-horizon instability must always be reported explicitly, but it does **not** automatically erase the 1–10 year comparison goal. It changes the **strength of interpretation**, not the existence of the comparison.

Every horizon-specific row exported from Stage 5, Stage 7, Stage 8, Stage 9, and Stage 10 must also carry a second machine-readable evidence-dependence field:

- `horizon_evidence_class = directly_observed_data_supported`
- `horizon_evidence_class = partly_model_dependent`
- `horizon_evidence_class = mostly_extrapolated`

Preferred mapping is:

- **PNU**
  - 1 year → `directly_observed_data_supported`
  - 2 years → `partly_model_dependent`
  - 3+ years → `mostly_extrapolated`

- **SNU**
  - 1–2 years → `directly_observed_data_supported`
  - 3–5 years → `partly_model_dependent`
  - 6+ years → `mostly_extrapolated`

- **merged**
  - 1–2 years → `directly_observed_data_supported`
  - 3–5 years → `partly_model_dependent`
  - 6+ years → `mostly_extrapolated`

A parallel reporting field `claim_restriction_flag` must then use:

- `primary_claim_allowed`
- `secondary_or_sensitivity_only`
- `projection_only`
- `projection_plus_prior_sensitive` when a Bayesian tail is both projection-dominant and materially prior-sensitive

Rows labeled `mostly_extrapolated` or `projection_only` may not support primary claims. Rows labeled `partly_model_dependent` may support only qualified secondary or sensitivity claims. Only rows labeled `directly_observed_data_supported` are eligible for the strongest manuscript claims, and even then they remain subject to the other stage-specific guardrails.

### 3.4 Stage logic rule

The workflow is staged for scientific reasons, not only for coding convenience:

- Stage 1 = backbone lock
- Stage 2 = follow-up maturity
- Stage 3 = KM benchmark classifier
- Stage 4 = timing-difference separation
- Stage 5 = individualized no-cure comparator
- Stage 6 = cure appropriateness screening
- Stage 7 = frequentist cure block
- Stage 8 = Bayesian cure block
- Stage 9 = remission sensitivity
- Stage 10 = unified comparison / final scientific answer

### 3.5 Stage 8 boundary rule

Stage 8 is **Stage 8 only**.  
It does **not** replace:

- Stage 2 follow-up maturity,
- Stage 4 timing-difference assessment,
- Stage 6 screening,
- Stage 9 remission sensitivity,
- or Stage 10 unified comparison.

### 3.6 Integrated design refinements derived from previous results

The following project refinements are now part of the standing specification:

1. **Revise the merged completeness block**
   - Recalculate `percentage_method`, `CCI`, and `SPT` for merged summaries using the site-specific administrative end date.
   - Pre-revision merged completeness values should not be used as manuscript-ready values.

2. **Always present KM with reverse KM and the numbers-at-risk table**
   - This is especially important for follow-up interpretation and late-tail caution.

3. **Present transition-only and remission-sensitive views side by side where interpretation may change**
   - The main backbone still treats remission as censoring.
   - But the effect of remission handling must remain visible, especially for PNU.

4. **Use an interval-based main site effect when comparing PNU and SNU**
   - The main timing-separation analysis must use:
     - 0 to 1 year
     - 1 to 2 years
     - 2 to 5 years
   - `>5 years` is a tail-diagnostic extension only.
   - A fully time-varying site effect may be added only as a sensitivity/descriptive supplement; it does not replace the interval-based main analysis.

5. **Run common-window restricted analyses**
   - A common 1-year window and a common 2-year window should be used to determine whether cohort differences are already present early.

6. **Run cure-appropriateness screening before making strong cure claims**
   - Use Maller–Zhou, HSU, Xie, and RECeUS / RECeUS-AIC as the Stage 6 screening block.
   - Do not rerun that screening inside Stage 8.

7. **Keep supported horizons primary and expose model dependence explicitly**
   - PNU: 1 year primary-supported / directly observed
   - PNU: 2 years sensitivity / partly model dependent
   - SNU and merged: 1–2 years primary-supported / directly observed
   - SNU and merged: 3–5 years secondary / partly model dependent
   - 6+ years projection / mostly extrapolated

8. **Retain both separate-cohort and merged views in the primary structural analysis set**
   - Separate-cohort fits remain essential for interpretation.
   - Merged results should still be retained in both site-free and site-adjusted form.

9. **Preserve downstream comparability**
   - All later no-cure, cure, and Bayesian outputs must remain on the same horizon/threshold ruler used for KM.

10. **Interpret site cautiously**
    - Site remains a proxy for broader context, care pathway, selection, and follow-up structure until treatment/process variables become available.

11. **For merged overall and merged sex completeness summaries, retain both structural views but prefer site-adjusted reporting**
    - `site_adjusted` is the preferred reporting view.
    - `site_free` is retained only as a sensitivity structure.

12. **Merged site-specific follow-up summaries must be coherent with the single-site cohorts they represent**
    - Merged site panels should inherit support-tier labels from the matched single-site cohort where available.
    - Reverse-KM and horizon-summary quantities for merged site panels should match the corresponding single-site overall summaries up to numerical tolerance.

13. **Carry explicit late-horizon site-composition descriptors and warnings in merged completeness outputs**
    - Site composition should be described for merged overall and sex panels.
    - Warning escalation should be rule-based and should not automatically trigger on primary-supported horizons.

14. **Preserve sparse-tail masking and preferred-reporting metadata as explicit outputs**
    - Stage 2 plots should be generated from exported source-of-truth tables.
    - Preferred-reporting flags and masking flags must remain inspectable in exported metadata.

## 4. Canonical common backbone

This section preserves the common backbone that all main model classes must share.

### 4.1 Datasets

Analyze the following in parallel:

- `PNU`
- `SNU`
- `merged`

Separate-cohort analyses remain essential for interpretation.

The merged dataset stays in the **primary structural analysis set**, not as a discarded afterthought.

### 4.2 Person identifier

Use:

- `site + id`

as the unique person identifier.

### 4.3 Time variables

Use:

- `days_followup` as the default analysis time variable
- `time_year = days_followup / 365.25` for horizon-based reporting and prediction

### 4.4 Main event definition

Main event:

- `status_num == 1` only

This means **transition** is the only event in the common main analysis.

### 4.5 Main censoring definition

Main censoring:

- `status_num %in% c(0, 2)`

Thus, remission is treated as censoring in the main backbone.

A remission-sensitive competing-risk or multi-state reanalysis belongs to **Stage 9**, not to the main backbone.

**Locked risk-scale rule.** The project now uses one shared horizon / threshold comparison ruler but two locked risk scales:

- `transition_only_main`
  - event of interest: `status_num == 1`
  - censoring: `status_num %in% c(0, 2)`
  - remission treated as censoring
  - used by the main backbone, Stage 3, Stage 5, Stage 6, Stage 7, Stage 8A, and the **primary** Stage 10 comparison

- `transition_cif_competing`
  - event of interest: `status_num == 1`
  - competing event: `status_num == 2`
  - right censoring: `status_num == 0`
  - used by Stage 8B, Stage 9, and the **remission-aware supplementary** Stage 10 comparison

Direct model-class comparison is allowed only **within the same risk scale**. Cross-scale comparison may be shown descriptively to quantify how remission handling changes interpretation, but it must never be presented as a fair one-to-one competition among models.


### 4.6 Common horizon grid

Use the same common prediction horizons for all model classes:

- 1, 2, ..., 10 years

### 4.7 Common threshold grid

Use the same pre-specified risk thresholds at each horizon for:

- KM benchmark classifiers
- no-cure models
- frequentist cure models
- Bayesian cure models

No model class may use a custom threshold grid if direct comparison is intended.

### 4.8 Common covariates

Use:

- `age_exact_entry`
- `sex_num`
- `site` where applicable in merged analyses

Create:

- `age_s = (age_exact_entry - mean(age_exact_entry)) / (2 * sd(age_exact_entry))`

### 4.9 Common formula registry

#### PNU / SNU

- Base: `age_s + sex_num`
- Interaction: `age_s + sex_num + age_s:sex_num`

#### merged

- Base: `age_s + sex_num`
- Interaction: `age_s + sex_num + age_s:sex_num`
- Site-added: `age_s + sex_num + site`
- Site + interaction: `age_s + sex_num + age_s:sex_num + site`

### 4.10 Missing-data rule

Assume there are no missing values in the main scripts unless a later explicit exception is introduced.

### 4.11 KM role

KM has **two roles** in this project:

1. cohort-level descriptive baseline
2. formal horizon-specific benchmark classifier for false-positive comparison

Use **overall KM** as the primary benchmark classifier.

Subgroup KM may be used only as a secondary/sensitivity benchmark when:

- subgroup strata are pre-specified, and
- sample size is adequate

KM is a benchmark, not a fully individualized regression model.

### 4.12 Output-unit rule

For every model class, calculate and export both:

- survival
- risk

### 4.13 Output-saving rule

- save fitted model objects as `.rds` when appropriate
- export key numbers as CSV
- write all outputs into one designated `export_path` folder
- do not scatter outputs across multiple subfolders unless absolutely necessary

---

## 5. Shared model menu and shared evaluation ruler

### 5.1 Model menu

#### KM benchmark
- overall KM (primary benchmark classifier)
- subgroup KM (secondary/sensitivity only, if justified)

#### Non-cure MLE
- exponential
- Weibull
- lognormal
- loglogistic
- Cox PH

#### Frequentist mixture cure
- exponential latency
- Weibull latency
- lognormal latency
- loglogistic latency
- semiparametric Cox-latency cure
- AFT-latency cure

#### Bayesian mixture cure
- Stage 8A main transition-only cure branch with exponential / Weibull / lognormal / loglogistic latency families
- Stage 8B remission-sensitive competing-risk cure extension

### 5.2 Cure-model extension hierarchy

AFT-latency cure is now part of the **required structural comparator set**, because AFT-latency and Cox-latency cure models can behave differently in the tail, particularly when follow-up sufficiency differs across covariate regions. Flexible parametric AFT cure ideas are acknowledged as methodologically relevant, but they remain secondary extension paths unless directly implemented in a later stage while preserving the shared comparison ruler.

### 5.3 Common comparison ruler

All direct model-class comparisons must preserve:

- the same event or risk-scale definition,
- the same horizon grid,
- the same threshold grid,
- the same output structure.

The project therefore uses **one shared horizon / threshold ruler** but does **not** authorize pooling across incompatible risk scales.

Every long-format table exported for downstream comparison should therefore carry an explicit `risk_scale` field so that primary transition-only results and remission-aware supplementary results cannot be merged accidentally.

### 5.4 Shared evaluation targets

For each model, horizon, and threshold where applicable, the shared comparison framework should support extraction of:

- survival
- risk
- confidence interval or credible interval
- discrimination summaries appropriate for censored data
- calibration summaries
- Brier-type summaries
- positive classification rate
- false-positive count
- false-positive burden
- false positives per 100
- PPV where meaningful
- net benefit
- risk-set size
- instability marker
- `support_tier`
- `horizon_evidence_class`
- `claim_restriction_flag`
- `risk_scale`

### 5.5 Common outputs to save from all model classes

All model classes should support export of the following common output families where applicable:

- fitted object (`.rds`)
- model registry table
- coefficient table
- fit summary table
- 1–10 year survival table
- 1–10 year risk table
- threshold-based classification table
- false-positive comparison table
- discrimination summary
- calibration summary
- Brier / IBS summary
- clinical-usefulness summary
- stability / risk-set summary
- reporting-governance summary
- plot-ready standardized visualization data

Tables and plots should be reproducible later from the saved fitted objects and exported prediction tables.


---

## 6. Stage-by-stage specification

Each stage below is structured under the same four headings:

1. **What this stage should do**
2. **What quantities, estimates, or calculations are needed**
3. **Reference basis**
4. **What later R implementation must preserve**

---

## 6.1 Stage 1 — Lock the estimand, datasets, and common comparison grid

### 6.1.1 What this stage should do

Freeze one common analytic backbone before any model fitting.

This stage exists to prevent post hoc drift in:

- event definition
- horizon choice
- threshold choice
- formula choice
- interpretation hierarchy

Stage 1 is a **backbone lock** stage, not a model-fitting stage.

### 6.1.2 What quantities, estimates, or calculations are needed

At minimum:

- dataset counts
- event counts
- remission counts
- censoring counts
- person-time summary
- age / sex / site summary
- modeling registry
- formula registry
- horizon registry
- threshold registry
- metadata registry
- saved scaling constants
- analysis datasets prepared for later stages

### 6.1.3 Reference basis

- Sections 2–5 of this integrated master specification
- `2008_A weakly informative default prior distr_Gelman et al..pdf`
- `R1. Gelman A. Scaling regression inputs by dividing by two standard deviations. Stat Med. 2008. (PubMed).pdf`

### 6.1.4 What later R implementation must preserve

- the dataset split registry
- the 1–10-year horizon vector
- the threshold vector
- the saved age-scaling constants
- the formula registry
- the metadata needed by all downstream stages

---

## 6.2 Stage 2 — Describe follow-up maturity before any cure claim

### 6.2.1 What this stage should do

Evaluate follow-up maturity, not cure.

For each dataset:

- estimate the KM curve
- estimate the reverse-KM curve
- produce numbers-at-risk tables
- summarize event / censoring / remission composition over time
- summarize follow-up completeness
- flag late-horizon instability

For merged analyses:

- recompute completeness using the **site-specific administrative end date**
- retain both `site_free` and `site_adjusted` completeness structures for merged overall and merged sex summaries
- treat `site_adjusted` as the **preferred reporting view**
- retain `site_free` as sensitivity only
- use `site_adjusted` completeness for merged site-stratified summaries

Transition-only and remission-sensitive summaries should be presented side by side, especially for PNU.

### 6.2.2 What quantities, estimates, or calculations are needed

At minimum:

- risk set size by horizon
- cumulative event counts by horizon
- cumulative censoring counts by horizon
- reverse-KM follow-up summary
- `percentage_method`
- `CCI`
- `SPT`
- late-horizon instability flags
- merged site-contribution summary by horizon where relevant
- `analysis_structure` fields distinguishing `single_site`, `site_free`, and `site_adjusted`
- `preferred_for_reporting` and `reporting_priority` for merged overall and sex summaries
- sparse-tail masking fields for plot-ready completeness outputs
- site-dominance descriptors and warning fields for merged overall / sex summaries
- support-tier labels for horizons and panels
- validation logs that check merged site panels against the corresponding single-site summaries

These are indicators of tail maturity, not direct evidence of cure.

### 6.2.3 Reference basis

- Sections 2–5 of this integrated master specification
- `R4. Betensky RA. Measures of follow-up in time-to-event studies. Clin Trials. 2015. (PubMed).pdf`
- `R20. Xue et al. New methods for estimating follow-up rate. 2017. (To verify).pdf`

### 6.2.4 What later R implementation must preserve

- site-specific administrative end-date lookup
- follow-up summary tables by dataset and sex
- reverse-KM summaries
- long-format maturity / instability outputs reusable in later reporting
- `support_tier_standard` or an equivalent machine-readable support label
- merged site-panel support labels matched to the corresponding single-site cohort
- explicit preferred-reporting and plot-masking fields
- exported validation summaries and issue rows

## 6.3 Stage 3 — Operationalize KM as the formal horizon-specific benchmark classifier

### 6.3.1 What this stage should do

Use **overall KM** as the primary benchmark classifier.

At each horizon `k`:

- compute `KM risk = 1 - S_KM(k)`
- apply the same threshold rule used for model-based predictions

This gives a simple, transparent, deliberately non-individualized benchmark.

### 6.3.2 What quantities, estimates, or calculations are needed

For each horizon and threshold:

- KM survival
- KM risk
- KM-based positive classification rate
- false-positive count
- false-positive burden
- false positives per 100 where useful
- specificity / PPV where meaningful
- benchmark positive pattern (`all_high_risk`, `all_low_risk`, or mixed) if the implementation exports it

### 6.3.3 Reference basis

- Sections 2–5 of this integrated master specification
- `R11. Vickers AJ, Elkin EB. Decision curve analysis- A novel method for evaluating prediction models. Med Decis Making. 2006..pdf`
- `R12. Vickers AJ, Cronin AM, Elkin EB, Gonen M. Extensions to decision curve analysis, a novel method for evaluating diagnostic tests, prediction models and molecular markers. BMC Med Inform Decis Mak. 2008.pdf`

### 6.3.4 What later R implementation must preserve

- a KM risk table by horizon
- a threshold-crossing function
- a long-format benchmark classification table aligned to later non-cure, cure, and Bayesian outputs

---

## 6.4 Stage 4 — Quantify the PNU–SNU timing difference before interpreting cure

### 6.4.1 What this stage should do

Separate timing difference from cure interpretation using an **interval-based main analysis**.

The purpose is to determine whether PNU’s earlier transition pattern and SNU’s later onset pattern are already visible in short and intermediate horizons where follow-up remains interpretable.

The **main** timing-separation analysis must use a piecewise site effect with intervals:

- 0–1 year
- 1–2 years
- 2–5 years

The `>5 years` segment must be treated as a **projection-dominant tail diagnostic only**. It is not part of the main early/intermediate timing-separation contrast and must not drive the primary interpretation of site timing difference.

This Stage 4 interval backbone does **not** replace the project-wide annual reporting grid used later for horizon-specific prediction, calibration, Brier, decision, and hazard-plausibility summaries. Downstream Stages 7–10 remain anchored to the common annual horizon grid of **1, 2, ..., 10 years from cohort entry**.

A fully time-varying site effect may be added only as a sensitivity or descriptive extension; it does **not** replace the interval-based main analysis.

Also run aligned-window restricted analyses:

- common 1-year restricted analysis
- common 2-year restricted analysis

This stage should answer:

> Is the PNU–SNU difference already present early, or is it mainly created by the longer SNU tail?

### 6.4.2 What quantities, estimates, or calculations are needed

At minimum:

- interval-specific site contrasts for:
  - 0–1 year
  - 1–2 years
  - 2–5 years
- a separate `>5 years` tail-diagnostic contrast
- 1-year restricted-window contrasts
- 2-year restricted-window contrasts
- hazard-pattern summaries
- explicitly labeled outputs:
  - `supported_short_horizon_contrast`
  - `intermediate_horizon_contrast`
  - `tail_diagnostic_contrast`

### 6.4.3 Reference basis

- Sections 2–5 of this integrated master specification
- `R17. Bullement A, et al. Evaluation of survival extrapolation in immuno-oncology models. BMC Med Res Methodol. 2020.pdf`
- `R18. Kearns B, et al. How uncertain is the survival extrapolation? 2020.pdf`

### 6.4.4 What later R implementation must preserve

- split-time or start-stop data structures for piecewise analyses using:
  - 0–1 years
  - 1–2 years
  - 2–5 years
  - optional `>5 years` tail-diagnostic extension
- common 1-year and 2-year restricted datasets
- explicit output labels distinguishing supported early/intermediate contrasts from projection-dominant tail diagnostics

---



## 6.5 Stage 5 — Fit the non-cure model block as the main individualized comparator

### 6.5.1 What this stage should do

Fit the non-cure MLE models as the main individualized no-cure baseline comparator:

- exponential
- Weibull
- lognormal
- loglogistic
- Cox PH

The primary output is not the coefficient table alone but the **subject × horizon risk matrix**.

### 6.5.2 What quantities, estimates, or calculations are needed

At each horizon:

- model family label
- maximized log-likelihood (`logLik`)
- effective parameter count
- AIC / BIC
- convergence status
- individualized `S(k)` and `1 - S(k)`
- time-dependent AUC (prefer cumulative/dynamic for fixed-horizon interpretation)
- Harrell’s C or Uno’s C where appropriate
- calibration summaries
- Brier score / scaled Brier / IBS
- positive classification rate
- false-positive count
- false-positive burden
- PPV / FDP where meaningful
- decision-curve net benefit

Interpretation priority remains:

- PNU: 1 year primary
- SNU and merged: 1–2 years primary
- 5 years secondary
- 10 years projection

Do not choose a single best non-cure family using AIC/BIC alone when tail behavior differs meaningfully.

### 6.5.3 Reference basis

- `R13. McLernon DJ, Giardiello D, Van Calster B, et al. Assessing performance and clinical usefulness in prediction models with survival outcomes. Ann Intern Med. 2023.pdf`
- `R3. Park SY, Park JE, Kim H, Park SH. Review of statistical methods for evaluating the performance of survival or other time-to-event prediction models. Korean J Radiol. 2021. (PubMed).pdf`
- `R2. Heagerty PJ, Lumley T, Pepe MS. Time-dependent ROC curves for censored survival data and a diagnostic marker. Biometrics. 2000. (PubMed).pdf`
- `R16. Beyene KM, Chen D-G. Evaluating predictive accuracy of prognostic models for censored time-to-event data analysis in clinical trials. (Chapter).pdf`
- `R11. Vickers AJ, Elkin EB. Decision curve analysis- A novel method for evaluating prediction models. Med Decis Making. 2006..pdf`
- `R12. Vickers AJ, Cronin AM, Elkin EB, Gonen M. Extensions to decision curve analysis, a novel method for evaluating diagnostic tests, prediction models and molecular markers. BMC Med Inform Decis Mak. 2008.pdf`

### 6.5.4 What later R implementation must preserve

- long-format subject × horizon prediction table as source of truth
- family-level fit summary table with `logLik`, effective parameter count, AIC / BIC, and convergence status
- horizon-specific performance table
- threshold-based classification table
- false-positive burden table
- decision-curve outputs
- supported/secondary/projection labels for every horizon
- joinable inputs for the Stage 7 family-matched cure-versus-non-cure fit-contrast block

---

## 6.6 Stage 6 — Cure appropriateness screening (frequentist screening only)

This section incorporates the full frequentist-only Stage 6 screening specification and is the main substantive update in this consolidated master file.

### 6.6.1 What this stage should do

Stage 6 is the **cure-appropriateness screening** stage.

Its purpose is to determine whether cure modeling is plausible and sufficiently supported by the observed follow-up structure **before** substantive cure-model fitting is interpreted.

Stage 6 is a screening layer only.

It is:

- **not** a final model-comparison stage,
- **not** a replacement for Stage 7 or Stage 8,
- and **not** a Bayesian modeling stage.

Accordingly, Bayesian mixture-cure models are explicitly excluded from Stage 6.

Stage 6 does **not** make the final substantive claim that a cure model is true.

Instead, Stage 6 assigns a **carry-forward cure-model eligibility flag**:

- `supportive`
- `equivocal`
- `unsupportive`

for later interpretation of both frequentist and Bayesian cure fits.

#### Document-placement rule

Within the documentation hierarchy:

- this integrated master file now contains both the shared backbone and the detailed operational Stage 6 rules,
- Stage 8 Bayesian detail remains in the separate Bayesian companion document,
- and any legacy Stage 6 references in older documents should be interpreted through this integrated file.

### 6.6.2 What quantities, estimates, or calculations are needed

#### Backbone inherited by Stage 6

Unless explicitly overridden below, Stage 6 inherits the common backbone:

- datasets: PNU, SNU, merged
- merged structural variants: site-free and site-adjusted
- main event definition: transition only
- remission treated as censoring in the main Stage 6 backbone
- time scale in years from baseline using the common conversion
- prediction horizons: 1–10 years
- threshold grid inherited unchanged from the common backbone
- standardized export structure keyed by dataset variant

#### Dataset variants to screen

Stage 6 must evaluate:

- `PNU`
- `SNU`
- `merged__site_free`
- `merged__site_adjusted`

For the merged site-adjusted variant, raw tail-based quantities are inherited from the same merged cohort.

In practice:

- raw-tail summaries,
- Maller–Zhou quantities,
- Xie quantities,
- and RECeUS quantities

are inherited from the `merged__site_free` branch unless a method itself requires site-adjusted re-estimation.

The HSU working-model branch is rerun for `merged__site_adjusted`, because it is working-model-dependent.

#### Methods included in Stage 6

At minimum, Stage 6 must include the following **frequentist screening tools only**:

- Maller–Zhou `d_n`
- Maller–Zhou `alpha_n` / `q_n` summary quantities
- HSU sup-score test, or a clearly documented pragmatic approximation
- Xie sufficient-follow-up test
- RECeUS with a fixed Weibull family
- RECeUS-AIC across the pre-specified candidate parametric latency families

These are screening tools, not final truth.

Their role is to summarize:

- cure signal,
- follow-up sufficiency,
- and cure-model appropriateness

before substantive cure fitting.

#### Bayesian exclusion rule

The following are forbidden inside Stage 6:

- Bayesian likelihood fitting
- Bayesian priors
- Bayesian posterior summaries
- Bayesian admissibility or retention rules
- any Bayesian rule that changes the Stage 6 screening flag

All Bayesian modeling belongs to Stage 8 only.

### 6.6.3 Stage 6 decision philosophy

Stage 6 should prioritize a **joint cure-appropriateness screen** over pure presence-only or pure follow-up-only screens.

Accordingly:

- **RECeUS-AIC is the primary Stage 6 decision gate**
- **Maller–Zhou `d_n` and HSU are presence-only supporting tests**
- **Xie is a contradiction-oriented follow-up check**
- **`alpha_n` / `q_n` are descriptive tail summaries only**

Presence-only evidence may soften an otherwise negative conclusion to `equivocal`, but it must not by itself create a final `supportive` conclusion.

Likewise, Xie may downgrade an otherwise positive conclusion to `equivocal`, but it must not serve as a stand-alone positive gate.

#### Role table

| Component | Role in Stage 6 | Reporting role |
|---|---|---|
| `RECeUS-AIC` | Primary cure-model eligibility gate | Main flag |
| `RECeUS-Weibull` | Sensitivity gate | Sensitivity row |
| `HSU` family-set result | Presence-only modifier | Supportive/equivocal modifier |
| `Maller–Zhou d_n` | Presence-only modifier | Supportive/equivocal modifier |
| `Xie` | Follow-up contradiction check | Downgrade to equivocal |
| `Maller–Zhou alpha_n / q_n` | Descriptive tail-follow-up summary | Table/figure only |

### 6.6.4 Primary and secondary Stage 6 rules

#### Primary gate

The primary Stage 6 gate is `RECeUS-AIC`.

`RECeUS-AIC` selects the best-supported latency family from the pre-specified candidate family set using AIC and then classifies cure-model appropriateness using the pre-specified RECeUS thresholds.

The fixed-family Weibull RECeUS result is retained as a **sensitivity analysis**, not as the main decision gate.

#### Presence-only support: Maller–Zhou and HSU

The following are treated as presence-only support:

- Maller–Zhou `d_n`
- HSU sup-score statistic

A positive result from at least one of these methods indicates that the data contain some evidence compatible with a cure fraction, but it does **not** by itself establish that cure modeling is sufficiently appropriate for substantive interpretation.

#### HSU must not be Weibull-only

In this project, HSU is **not** restricted to a single Weibull-only working latency model.

The preferred Stage 6 HSU working-family set is:

- `weibull`
- `lognormal`
- `loglogistic`

#### Required HSU outputs

For each dataset variant, Stage 6 should export:

- one family-specific HSU result row for each working family
- one aggregated HSU family-set result row

Recommended family-specific method names:

- `hsu_supscore_weibull_approx`
- `hsu_supscore_lognormal_approx`
- `hsu_supscore_loglogistic_approx`

Preferred aggregated row:

- `hsu_supscore_familyset_approx`

#### Preferred HSU aggregation rule

The preferred HSU aggregate is a **bootstrap-calibrated family-set supremum test** over the prespecified HSU working-family set.

That is, the family-set statistic should be based on the maximum family-specific HSU statistic across:

- Weibull
- lognormal
- loglogistic

with a bootstrap-calibrated p-value for the family-set maximum.

If full family-set bootstrap calibration is unavailable, an acceptable fallback is:

- export all family-specific p-values,
- apply a prespecified multiplicity adjustment,
- and derive the aggregated HSU flag from the adjusted family-set result.

Under either implementation, the official HSU contribution to the final Stage 6 decision is the **aggregated family-set flag**, not an unadjusted single-family minimum p-value.

#### RCT-arm rule not imported here

The RCT-specific extension of RECeUS that evaluates appropriateness by randomized arm and may use an “at least one randomized arm” decision rule is **not** adopted as an operational Stage 6 rule in this project.

This project analyzes **observational site-cohort variants**:

- `PNU`
- `SNU`
- `merged__site_free`
- `merged__site_adjusted`

not randomized treatment arms.

The RCT-arm rule may be cited as an external methodological reference, but it must not be copied directly into the operational Stage 6 classification rule here.

#### Follow-up contradiction check

The Xie sufficient-follow-up test is used as a contradiction-oriented follow-up check.

A contradictory Xie result may downgrade an otherwise supportive Stage 6 conclusion, but Xie must not be used as the sole positive gate.

#### Descriptive tail summaries

The following remain descriptive only:

- `alpha_n`
- `q_n`
- plateau-length summaries
- KM tail survival
- other related extreme-tail summaries

They may inform interpretation, but they must **not** by themselves generate a final supportive Stage 6 conclusion.

### 6.6.5 Operational thresholds

Use the common Stage 6 screening significance level:

- `alpha_screening = 0.05`

#### RECeUS classification thresholds

Use the following project thresholds:

- `supportive` if `cure_fraction_hat > 0.025` and `receus_ratio_hat < 0.05`
- `equivocal` if:
  - `cure_fraction_hat > 0.025` and `receus_ratio_hat < 0.10`, or
  - `cure_fraction_hat > 0.010` and `receus_ratio_hat < 0.10`
- `unsupportive` otherwise

These exact bands should be fixed in code and exported with the Stage 6 metadata registry.

#### Presence-only significance thresholds

Use `alpha_screening = 0.05` for:

- Maller–Zhou `d_n`
- HSU family-specific p-values
- HSU family-set aggregate p-value

#### Xie contradiction threshold

Use `alpha_screening = 0.05` for Xie.

A contradictory Xie result is defined as:

- `p_Xie < alpha_screening`

under the chosen Xie direction convention.

### 6.6.6 Final Stage 6 classification rule

Let:

- `R` = RECeUS-AIC classification (`supportive`, `equivocal`, `unsupportive`)
- `P_MZ` = whether Maller–Zhou `d_n` is positive
- `P_HSU` = whether the aggregated HSU family-set result is positive
- `P` = `P_MZ OR P_HSU`
- `X` = whether Xie gives a contradictory follow-up signal

Then the final Stage 6 classification is:

1. **supportive** if `R = supportive` and `X = FALSE`
2. **equivocal** if `R = supportive` and `X = TRUE`
3. **equivocal** if `R = equivocal`, regardless of `P`
4. **equivocal** if `R = unsupportive` but `P = TRUE`
5. **unsupportive** if `R = unsupportive` and `P = FALSE`

This ensures that:

- RECeUS-AIC remains the primary gate,
- presence-only tests may rescue an otherwise negative screen only to `equivocal`,
- Xie acts as a downgrade signal rather than a stand-alone positive gate,
- descriptive tail summaries do not dominate the final decision,
- HSU is incorporated as a family-set presence-only block rather than a single-family shortcut.

### 6.6.7 Preferred exported naming convention

To avoid implying that Stage 6 is the final substantive verdict, the preferred cross-stage names are:

- `cure_model_eligibility_flag` instead of `final_decision_flag`
- `primary_gate_method`
- `primary_gate_flag`
- `receus_primary_class`
- `presence_modifier_flag`
- `cure_presence_support_flag`
- `followup_contradiction_flag`
- `followup_not_contradicted_flag`
- `descriptive_tail_summary_flag`
- `screening_context`

Recommended values include:

- `primary_gate_method = RECeUS-AIC`
- `screening_context = observational_cohort`

In this revised specification, the paired decomposed fields

- `receus_primary_class`
- `cure_presence_support_flag`
- `followup_not_contradicted_flag`

are **canonical carry-forward exports**, not merely convenience aliases. They are retained alongside the more general semantic fields because later Stage 7 / Stage 8 logic repeatedly uses them directly.

For backward compatibility, legacy outputs may temporarily retain:

- `final_decision_flag`

but the preferred semantic name is `cure_model_eligibility_flag`.

### 6.6.8 Required exports

At minimum, Stage 6 must export:

- a variant registry keyed by `dataset_key`
- a follow-up sufficiency summary table
- a long-format screening-method results table
- a RECeUS candidate-family summary table
- an HSU family-set / familywise summary table, or equivalent familywise rows in the long-format results
- a wide Stage 6 screening summary table
- a carry-forward flag table for later stages
- a metadata registry documenting thresholds, bootstrap settings, candidate families, and implementation notes
- a reusable Stage 6 bundle
- optional compact fitted objects and bootstrap distributions

#### Minimum long-format method names

The long-format screening-method table should include at minimum:

- `maller_zhou_dn`
- `maller_zhou_alpha_n`
- `xie_sufficient_followup`
- `receus_weibull`
- `receus_aic`
- `hsu_supscore_weibull_approx`
- `hsu_supscore_lognormal_approx`
- `hsu_supscore_loglogistic_approx`
- `hsu_supscore_familyset_approx`

#### Minimum carry-forward fields

The carry-forward flag table should preserve at minimum:

- `dataset_key`
- `cure_model_eligibility_flag`
- `primary_gate_method`
- `primary_gate_flag`
- `receus_primary_class`
- `presence_modifier_flag`
- `cure_presence_support_flag`
- `followup_contradiction_flag`
- `followup_not_contradicted_flag`
- `descriptive_tail_summary_flag`
- `supporting_methods`
- `contradicting_methods`
- `screening_note`
- `common_horizon_vector`
- `common_threshold_vector`
- `screening_context`

Legacy loader aliases may additionally be derived where helpful:

- `stage6_final_class = cure_model_eligibility_flag`
- `carry_forward_stage8 = cure_model_eligibility_flag != "unsupportive"`

but these legacy aliases are now secondary to the standardized canonical fields above.

#### Legacy compatibility note

If an existing implementation still exports:

- `final_decision_flag`

and only a limited HSU row such as:

- `hsu_supscore_weibull_approx`

that may be retained for backward compatibility.

However, the preferred forward specification is the multi-family HSU set plus the aggregated family-set row.

### 6.6.9 Relationship to later stages

Every Stage 6 dataset-variant result must be merged forward into later stages by `dataset_key`.

Stage 7 and Stage 8 may use these as:

- interpretation flags
- stratification labels
- exclusion warnings
- projection-strength annotations

They must not reinterpret Stage 6 as a Bayesian model-comparison exercise.

#### Required Stage 8 downstream reading rule

Stage 8 Bayesian cure models are a stabilization and uncertainty layer **conditional on the Stage 6 carry-forward eligibility flag**.

Accordingly:

- dataset variants flagged `supportive` or `equivocal` remain in the primary Bayesian interpretation set
- dataset variants flagged `unsupportive` may still be fitted as sensitivity analyses, but they are not treated as primary substantive evidence for cure-like heterogeneity

Stage 8 may refine parameter estimation and uncertainty characterization, but it must not retroactively redefine Stage 6 screening status.

### 6.6.10 Interpretation boundary

Stage 6 does **not** decide which cure model is scientifically best overall.

It decides only whether cure modeling is:

- sufficiently supported,
- mixed / borderline,
- or insufficiently supported

under the frequentist screening layer defined here.

The direct comparison among:

- KM benchmark classifiers,
- no-cure models,
- frequentist cure models,
- Bayesian cure models

belongs to Stage 10, not to Stage 6.

### 6.6.11 One-paragraph Stage 6 summary

> **Stage 6 performs cure-appropriateness screening on PNU, SNU, merged site-free, and merged site-adjusted variants using frequentist screening tools only. Bayesian mixture-cure models are explicitly excluded from Stage 6 and first appear in Stage 8. Stage 6 uses RECeUS-AIC as the primary cure-model eligibility gate, retains fixed-family Weibull RECeUS as a sensitivity screen, uses Maller–Zhou `d_n` and a multi-family HSU family-set screen as presence-only supporting tests, uses Xie as a contradiction-oriented follow-up check, treats `alpha_n` and `q_n` as descriptive tail summaries only, exports standardized screening and carry-forward tables keyed by dataset variant, and passes these results forward as interpretation flags into both the frequentist cure block and the Bayesian cure block without repeating screening inside Stage 8.**

### 6.6.12 Reference basis

- Sections 2–6 of this integrated master specification
- `R5. Maller RA, Zhou S. Testing for the presence of immune or cured individuals in censored survival data. Biometrics. 1995. (PubMed).pdf`
- `R6. Hsu C-H, Todem D, Kim K. A sup-score test for the cure fraction in mixture models for long-term survivors. Biometrics. 2016. (PubMed).pdf`
- `R7. Xie P, Escobar-Bach M, Van Keilegom I. Testing for sufficient follow-up in censored survival data by using extremes. Biometrical Journal. 2024. (PubMed).pdf`
- `R8. Selukar S, Othus M. RECeUS - Ratio estimation of censored uncured subjects, a different approach for assessing cure model appropriateness in studies with long-term survivors. Stat Med. 2023. (PubMed).pdf`
- `R19. Kouadio C, et al. Detecting cure model appropriateness in randomized clinical trial data. 2025.pdf` (methodological context only; not imported as an operational observational-cohort rule)

### 6.6.13 What later R implementation must preserve

- standardized variant registry
- standardized long-format screening rows
- standardized carry-forward flag table
- compatibility with Stage 7 and Stage 8 by `dataset_key`
- explicit metadata about thresholds, bootstrap settings, and candidate families
- backward-compatible field names where needed

---

## 6.7 Stage 7 — Fit the frequentist cure block

### 6.7.1 What this stage should do

Fit the frequentist cure models after Stage 6 screening is complete.

Include as the **required structural comparator set**:

- parametric mixture cure with exponential latency
- Weibull latency
- lognormal latency
- loglogistic latency
- semiparametric Cox-latency mixture cure
- AFT-latency mixture cure

Cox-latency and AFT-latency are required structural comparators in this revised specification, not optional targeted sensitivities.

More elaborate extensions, such as flexible parametric AFT cure implementations, may be added as targeted secondary sensitivities only if they preserve the shared backbone and exported comparison ruler.

For merged analyses, the cure block must preserve explicit site-placement structures:

- `site_in_neither`
- `site_in_incidence_only`
- `site_in_latency_only`
- `site_in_both`

These site-placement terms are structural decomposition devices for interpreting where cohort differences operate. They must **not** be written up as causal treatment effects.

The goal is not only to estimate cure fractions, but to compare how allowing cure-like heterogeneity changes:

- overall risk,
- susceptible-only survival,
- false-positive burden,
- and clinical usefulness.

This stage is also the main frequentist contributor to the project’s indirect within-CHR-P heterogeneity triangulation and to the future-method recommendation about cure-aware modeling.

### 6.7.2 What quantities, estimates, or calculations are needed

At minimum:

- incidence coefficients
- latency coefficients
- site-placement structure labels for merged fits
- cure fraction estimates
- `susceptible_fraction = 1 - cure_fraction`
- susceptible-only / uncured-only survival and risk on the annual 1–10-year grid
- optional uncured-only mean survival summary (`MSTu`) only when the uncured support is sufficiently identified and the result is not contradicted by follow-up screening; otherwise keep uncured-only reporting on the annual grid and do not force a single extrapolated mean
- overall survival / risk at 1–10 years
- `DeltaRisk_NC_C(t)` versus non-cure models
- threshold-based false-positive burden
- false positives per 100
- PPV / TPR where meaningful
- net benefit
- horizon-specific discrimination summaries whenever estimable
- horizon-specific calibration summaries whenever estimable
- horizon-specific Brier / IBS summaries whenever estimable
- hazard-shape plausibility exports on the common annual horizon grid:
  - `hazard_target`
  - `hazard_1y`
  - `hazard_2y`
  - `hazard_3y`
  - `hazard_4y`
  - `hazard_5y`
  - `hazard_6y`
  - `hazard_7y`
  - `hazard_8y`
  - `hazard_9y`
  - `hazard_10y`
  - `hazard_ratio_10y_vs_1y`
  - `shape_class`
- uncertainty intervals
- `support_tier`
- `horizon_evidence_class`
- `claim_restriction_flag`
- explicit `risk_scale = transition_only_main`

The main inferential weight remains on directly observed or otherwise supported horizons. Hazard-plausibility reporting must therefore follow the same annual 1–10-year grid as the rest of the project rather than a sparse anchor-year subset only.

The uncured-only / non-cure-fraction quantities are useful **supporting decomposition estimands** but remain separate from the main common-scale comparison ruler.

### 6.7.2A Family-matched cure-versus-non-cure likelihood-contrast block

For the fully parametric families that exist in both Stage 5 and Stage 7, the workflow must additionally export one family-matched cure-versus-non-cure fit-contrast block.

Required matched pairs are:

- exponential non-cure MLE versus exponential mixture cure MLE
- Weibull non-cure MLE versus Weibull mixture cure MLE
- lognormal non-cure MLE versus lognormal mixture cure MLE
- loglogistic non-cure MLE versus loglogistic mixture cure MLE

The purpose of this block is supportive rather than decisive. It helps quantify whether allowing a cure component improves within-family fit and whether that improvement coheres with the broader heterogeneity argument.

Minimum exported fields are:

- `dataset_key`
- `family_pair`
- `logLik_noncure`
- `logLik_cure`
- `delta_logLik_cure_minus_noncure`
- `LR_2delta_logLik`
- `delta_AIC_cure_minus_noncure`
- `delta_BIC_cure_minus_noncure`
- `lrt_calibration_status`
- optional `lrt_pvalue_bootstrap`
- `same_family_fit_gain_signal`
- `convergence_pair_flag`

#### Nonregular-inference rule

A naive asymptotic chi-square likelihood-ratio p-value must not be used as the default inferential basis for cure-versus-non-cure comparison.

The cure-null is a boundary / nonregular problem, and the current workflow therefore does not authorize routine textbook chi-square interpretation of `LR_2delta_logLik`.

If a p-value is reported, it must be clearly labeled as one of:

- `bootstrap_calibrated_nonregular_LRT`
- `not_reported_nonregular_problem`

The bootstrap-calibrated version is allowed only as an optional supportive sensitivity for fully parametric same-family pairs.

#### Interpretation rule

The family-matched fit-contrast block is part of the heterogeneity-triangulation evidence set.

It may support statements such as “allowing a cure component materially improved within-family fit,” but it may not by itself prove cure, replace Stage 6 screening, or outweigh supported-horizon risk and false-positive-burden evidence.

### 6.7.3 Reference basis

- Sections 2–6 of this integrated master specification
- `R5. Maller RA, Zhou S. Testing for the presence of immune or cured individuals in censored survival data. Biometrics. 1995. (PubMed).pdf`
- `R6. Hsu C-H, Todem D, Kim K. A sup-score test for the cure fraction in mixture models for long-term survivors. Biometrics. 2016. (PubMed).pdf`
- `R14A. Peng Y, Yu B. The parametric cure model. Cure Models - Methods, Applications, and Implementation. CRC Press. 2021-2022. (Textbook chapter).pdf`
- `R14B. Peng Y, Yu B. The semiparametric and nonparametric cure models. Cure Models - Methods, Applications, and Implementation. CRC Press. 2021-2022. (Textbook chapter).pdf`
- `R9. Cai C, Zou Y, Peng Y, Zhang J. smcure - An R-package for estimating semiparametric mixture cure models. Comput Methods Programs Biomed. 2012. (PubMed).pdf`
- `R10. Parsa M, Van Keilegom I. Accelerated failure time vs Cox proportional hazards mixture cure models. Statistical Papers. 2023. (Springer).pdf`
- `2022_An Accelerated Failure Time Cure Model w_Aida et al..pdf`
- `2025_Flexible Parametric Accelerated Failure_Akynkozhayev et al..pdf`
- `2024_A two-sample comparison of mean survival_Dobler and Musta.pdf` for uncured-only / `MSTu` supporting-estimand justification
- `2024_Cureit  An End-to-End Pipeline f_Whiting et al..pdf` for implementation-facing context only, not as a governing workflow source

### 6.7.4 What later R implementation must preserve

- cure-model registry with family, incidence link, latency type, and site-free/site-adjusted branch
- explicit merged site-placement labels:
  - `site_in_neither`
  - `site_in_incidence_only`
  - `site_in_latency_only`
  - `site_in_both`
- standardized prediction outputs
- Stage 6 carry-forward flag join by `dataset_key`
- convergence diagnostics
- long-format outputs directly joinable into Stage 10
- hazard-shape plausibility table
- family-matched cure-versus-non-cure fit-contrast table with `LR_2delta_logLik`, `delta_AIC_cure_minus_noncure`, `delta_BIC_cure_minus_noncure`, and `lrt_calibration_status`
- cure-model-only supporting decomposition table with `susceptible_fraction`, uncured-only annual summaries, and optional `MSTu`
- explicit `risk_scale` labels
- `support_tier`, `horizon_evidence_class`, and `claim_restriction_flag`


## 6.8 Stage 8 — Bayesian cure block (shared interface to Stage 8A and Stage 8B only)

**Important:** The full Stage 8 mathematical, prior, and export specification remains in the separate companion document:

- `Bayesian_Modeling_Specification_Stage8_REVISED_v5.md`

This integrated master file retains only the shared interface rules needed to keep Stage 8 aligned with the rest of the workflow.

### 6.8.1 What this stage should do

Stage 8 acts as a **substantive-prior, stabilization, uncertainty, and interpretation-governance layer**.

Its role is to:

1. stabilize estimation when MLE cure models are unstable
2. propagate uncertainty more transparently
3. allow carefully specified external substantive prior information
4. export Bayesian predictions on the same downstream comparison grid
5. keep the **Stage 8A transition-only branch** as the main Bayesian branch
6. add the **Stage 8B remission-sensitive competing-risk branch** as a visible supplementary extension rather than a silent rewrite of the backbone
7. require an explicit **anchor-informed versus no-external-information comparison** because integration of external meta-analysis is a central scientific feature of the Bayesian block
8. require a **prior-to-posterior incidence-shape update display** so the external anchor is shown as a starting point, not hidden as if it were the posterior itself
9. export explicit horizon-support, prior-tail, hazard-plausibility, anchor-value, and branch-delta tables so later interpretation cannot overclaim the tail or the prior

### 6.8.2 What quantities, estimates, or calculations are needed

At minimum for each admissible Bayesian fit:

- posterior cure-fraction and susceptible-fraction summaries
- posterior yearly survival / risk or CIF at 1–10 years
- posterior yearly uncured-only / susceptible-only survival and risk on the annual 1–10-year grid
- optional posterior `MSTu` only when support is adequate and follow-up screening does not contradict the required tail interpretation
- posterior `DeltaRisk_NC_C(t)` where defined
- posterior threshold-based false-positive burden
- false positives per 100
- PPV / TPR where estimable
- posterior net benefit / decision curves
- horizon-specific discrimination summaries whenever estimable
- horizon-specific calibration summaries whenever estimable
- horizon-specific Brier / IBS summaries whenever estimable
- hazard-shape plausibility exports on the common annual horizon grid:
  - `hazard_target`
  - `hazard_1y`
  - `hazard_2y`
  - `hazard_3y`
  - `hazard_4y`
  - `hazard_5y`
  - `hazard_6y`
  - `hazard_7y`
  - `hazard_8y`
  - `hazard_9y`
  - `hazard_10y`
  - `hazard_ratio_10y_vs_1y`
  - `shape_class`
- convergence diagnostics
- posterior predictive checks
- prior sensitivity summaries
- explicit `prior_branch` labels (`anchor_informed` or `neutral_no_external_info`)
- explicit `site_prior_family` labels when a site term exists
- `support_tier`
- `horizon_evidence_class`
- `claim_restriction_flag`
- `prior_tail_sensitive`
- explicit `risk_scale` labels
- explicit branch labels (`Stage8A` or `Stage8B`)
- a mandatory anchor-informed versus neutral delta table
- a mandatory prior-to-posterior incidence-shape source table
- a mandatory Stage 8B-exported Stage 8A-versus-Stage 8B delta table for remission-aware change quantification

The default hazard-plausibility export for Stage 8 is the full annual `1, 2, ..., 10 years` grid from cohort entry, not a reduced anchor-year subset only.

### 6.8.3 Shared Stage 8 interface rules

The Stage 8 Bayesian block must preserve the common backbone in:

- datasets
- person identifiers
- common horizon grid
- common threshold grid
- exported output structure

Stage 8A uses the main project scale:

- `risk_scale = transition_only_main`
- event = `status_num == 1`
- censoring = `status_num %in% c(0, 2)`

Stage 8B uses the remission-sensitive competing-risk scale:

- `risk_scale = transition_cif_competing`
- event = `status_num == 1`
- competing event = `status_num == 2`
- right censoring = `status_num == 0`

Stage 8B may not silently substitute for Stage 8A in primary reporting.

Cross-scale Stage 8A-versus-Stage 8B deltas are allowed only as **remission-aware change summaries**. They must not be misread as a within-scale model-ranking contest.

The mandatory no-external-information comparator for any anchored Bayesian fit is the neutral branch labeled:

- `prior_branch = neutral_no_external_info`

Legacy implementation aliases such as `neutral_weakly_informative` may still exist internally, but Stage 10 reporting should prefer the clearer no-external-information label.

When site coefficients are present in Stage 8, the preferred main prior family is:

- `normal_0_1_main`

with mandatory robustness sensitivity under:

- `student_t3_0_1_sensitivity`

Detailed prior encoding remains in the companion specification, but the exported metadata must preserve which site-prior family was actually used.

Any merged site coefficient in incidence, latency, or remission is a structural context proxy, not a causal treatment effect.

### 6.8.4 Stage 8 carry-forward rule

Stage 8 must **not** rerun Stage 6 screening.

It must carry forward Stage 6 screening results as interpretation flags.

The preferred upstream input fields are the standardized Stage 6 carry-forward fields listed in Section 6.6.8, including the decomposed canonical fields:

- `receus_primary_class`
- `cure_presence_support_flag`
- `followup_not_contradicted_flag`

Legacy aliases such as `stage6_final_class` or `carry_forward_stage8` may still exist internally, but they are convenience fields only and must remain semantically equivalent to the preferred standardized inputs.

### 6.8.5 Admissibility and retention rule

Stage 8 should not force a single Bayesian winner.

Instead:

1. exclude non-admissible fits
2. retain all admissible Bayesian fits
3. export all admissible fits to the later unified comparison engine
4. use internal Bayesian fit statistics as supportive evidence only

### 6.8.6 Reference basis

- `Bayesian_Modeling_Specification_Stage8_REVISED_v5.md`
- this integrated master specification for the shared interface rules
- `R14C. Peng Y, Yu B. Bayesian cure model. Cure Models - Methods, Applications, and Implementation. CRC Press. 2021-2022. (Textbook chapter).pdf`
- `2008_A weakly informative default prior distr_Gelman et al..pdf`
- `R1. Gelman A. Scaling regression inputs by dividing by two standard deviations. Stat Med. 2008. (PubMed).pdf`
- `2014_Systematic review and collaborative reca_Van Der Werf et al..pdf` when external age–sex incidence priors are encoded
- `2024_A two-sample comparison of mean survival_Dobler and Musta.pdf` for the auxiliary uncured-only / non-cure-fraction supporting block

### 6.8.7 What later R implementation must preserve

- shared horizon / threshold grid
- dataset compatibility with Stage 10
- Stage 6 carry-forward interpretation
- admissibility diagnostics
- exported long-format posterior prediction tables
- hazard-shape plausibility table
- prior-tail governance table
- anchor-informed versus neutral delta table
- prior-to-posterior incidence-shape source table
- Stage 8A-versus-Stage 8B delta table, exported by Stage 8B
- cure-model-only supporting decomposition table with `susceptible_fraction`, uncured-only annual summaries, and optional `MSTu`
- explicit `prior_branch`, `site_prior_family`, branch labels, and `risk_scale` labels

## 6.9 Stage 9 — Remission sensitivity

### 6.9.1 What this stage should do

Reanalyze remission as a competing event or multi-state sensitivity block.

This stage is especially important for **PNU**, where remission may materially change apparent long-horizon transition risk and may create a false cure-like signal.

Stage 9 is the **frequentist remission-sensitive mirror** of the Bayesian Stage 8B idea. It should stay on the same supplementary risk scale so that remission-sensitive frequentist and Bayesian outputs remain directly comparable.

### 6.9.2 What quantities, estimates, or calculations are needed

At minimum:

- `status_num == 0` reported separately
- `status_num == 2` reported separately
- transition cumulative incidence by horizon
- remission cumulative incidence by horizon
- absolute difference from the main censoring-based transition risk
- change in false-positive burden
- change in net benefit
- CIF-scale discrimination / calibration / Brier summaries where estimable
- `support_tier`
- `horizon_evidence_class`
- `claim_restriction_flag`
- `risk_scale = transition_cif_competing`

If feasible, include a simple multi-state structure:

- entry → transition
- entry → remission

### 6.9.3 Reference basis

- this integrated master specification
- `R12. Vickers AJ, Cronin AM, Elkin EB, Gonen M. Extensions to decision curve analysis, a novel method for evaluating diagnostic tests, prediction models and molecular markers. BMC Med Inform Decis Mak. 2008.pdf`
- `R13. McLernon DJ, Giardiello D, Van Calster B, et al. Assessing performance and clinical usefulness in prediction models with survival outcomes. Ann Intern Med. 2023.pdf`

### 6.9.4 What later R implementation must preserve

- cause-specific or subdistribution analysis code
- horizon-specific cumulative incidence predictions
- mirrored threshold-evaluation framework so the sensitivity results remain directly comparable with the main analysis
- `support_tier`, `horizon_evidence_class`, and `claim_restriction_flag`
- `risk_scale = transition_cif_competing`

## 6.10 Stage 10 — Build one unified comparison engine and one reporting structure


### 6.10.1 What this stage should do

After all model blocks are fitted, Stage 10 must build **five coordinated reporting objects**, not one undifferentiated table:

1. **one main common-scale comparison table** on `transition_only_main`
2. **one remission-aware comparison table** on `transition_cif_competing`
3. **one cure-model-only supporting decomposition block** kept separate from the main common-scale comparison table
4. **one remission-aware delta block** quantifying how explicit remission handling changes apparent cure-like signal and downstream classification burden
5. **one external-anchor value block** quantifying what changes when the external meta-analytic incidence anchor is removed and how the anchored incidence shape is updated by the posterior

The cure-model-only supporting decomposition block is where incidence–latency separation, susceptible-only summaries, optional uncured-only mean summaries, and merged site-placement decomposition are read. It must remain separate from the main common-scale comparison table so that decomposition outputs do not masquerade as the shared primary ruler.

Then answer the scientific question in this order:

1. **Main comparison on `transition_only_main`**  
   At supported horizons, how do KM, no-cure, frequentist cure, and Bayesian Stage 8A differ in false-positive burden and clinical usefulness?

2. **Cohort-level risk shift on `transition_only_main`**  
   How much does estimated cumulative risk change when cure-like heterogeneity is allowed?

3. **External-anchor value block**  
   How much do supported-horizon risk, false-positive burden, and net benefit change when the external incidence anchor is removed, and does the posterior meaningfully update away from the anchor pattern?

4. **Cure-model-only supporting decomposition block**  
   Are the cohort differences operating mainly through incidence, latency, both, or neither, and does the merged site-placement structure point toward `site_in_incidence_only`, `site_in_latency_only`, `site_in_both`, or `site_in_neither`?

5. **PNU–SNU timing difference**  
   Is the site difference already present in the common 1-year and 2-year windows, and how does the 2–5 year interval behave before the tail becomes projection-dominant?

6. **Tail instability and support governance**  
   Which contrasts are directly observed, partly model dependent, or mostly extrapolated?

7. **Remission-aware change block**  
   How much do remission-sensitive frequentist and Bayesian branches change apparent cure-like signal, false-positive burden, and decision-curve conclusions when remission is treated explicitly as a competing event?

8. **Final interpretation bucket**  
   Is the later cure-like signal more compatible with true cure heterogeneity, timing difference, remission handling, external-prior dependence, or projection-dominant tail behavior?

### 6.10.1A Heterogeneity-triangulation block

Stage 10 must additionally produce one heterogeneity-triangulation block.

This block is interpretive rather than competitive. It does not rank models. Its role is to summarize whether the total pattern of evidence is more compatible with clinically meaningful within-CHR-P heterogeneity than with a single homogeneous no-cure survival process.

Minimum fields are:

- `timing_signal`
- `followup_signal`
- `screening_signal`
- `latency_family_signal`
- `remission_signal`
- `anchor_dependence_signal`
- `same_family_fit_gain_signal`
- `uncured_only_support_signal`
- `triangulated_heterogeneity_bucket`
- `future_cure_modeling_recommendation_flag`

Preferred `triangulated_heterogeneity_bucket` values are:

- `indirect_support_for_within_CHR-P_heterogeneity`
- `mixed_but_suggestive`
- `timing_or_followup_explains_most`
- `insufficient_or_tail_driven`

Preferred `future_cure_modeling_recommendation_flag` values are:

- `recommended_to_consider`
- `conditionally_suggestive_needs_larger_data`
- `not_supported_by_current_evidence`

### 6.10.2 What quantities, estimates, or calculations are needed

For each model and horizon, the final unified tables should contain as far as available:

- survival
- risk or transition CIF
- confidence or credible interval
- time-dependent discrimination summaries
- Brier / IBS
- calibration summaries
- positive classification rate
- false-positive count
- false-positive burden
- false positives per 100
- PPV
- TPR
- net benefit
- net reduction in unnecessary intervention
- risk-set size
- instability marker
- `support_tier`
- `horizon_evidence_class`
- `claim_restriction_flag`
- `prior_tail_sensitive` where applicable
- `prior_branch` where applicable
- `site_prior_family` where applicable
- branch label
- risk_scale label

The **cure-model-only supporting decomposition block** must additionally contain where applicable:

- incidence coefficients
- latency coefficients
- cure fraction summaries
- `susceptible_fraction`
- susceptible-only / uncured-only survival and risk on the common annual 1–10-year grid
- optional `MSTu` when support is adequate
- `uncured_mean_support_flag` indicating whether a single-number uncured-only mean is supportable
- merged site-placement labels
- hazard-shape plausibility exports on the common annual horizon grid:
  - `hazard_target`
  - `hazard_1y`
  - `hazard_2y`
  - `hazard_3y`
  - `hazard_4y`
  - `hazard_5y`
  - `hazard_6y`
  - `hazard_7y`
  - `hazard_8y`
  - `hazard_9y`
  - `hazard_10y`
  - `hazard_ratio_10y_vs_1y`
  - `shape_class`

The **external-anchor value block** must include where estimable:

- `delta_risk_anchor_minus_neutral`
- `delta_cure_fraction_anchor_minus_neutral`
- `delta_false_positive_burden_anchor_minus_neutral`
- `delta_FP100_anchor_minus_neutral`
- `delta_NB_anchor_minus_neutral`
- `delta_PPV_anchor_minus_neutral`
- `delta_TPR_anchor_minus_neutral`
- prior-to-posterior incidence-shape source fields needed for plotting the anchor update

The **remission-aware delta block** must include where estimable:

- `delta_risk_8B_minus_8A`
- `delta_cure_fraction_8B_minus_8A`
- `delta_false_positive_burden_8B_minus_8A`
- `delta_FP100_8B_minus_8A`
- `delta_NB_8B_minus_8A`
- `delta_PPV_8B_minus_8A`
- `delta_TPR_8B_minus_8A`

and, where a frequentist remission-aware mirror exists, analogous Stage 9-versus-main delta fields.

The **heterogeneity-triangulation block** must additionally contain:

- `timing_signal`
- `followup_signal`
- `screening_signal`
- `latency_family_signal`
- `remission_signal`
- `anchor_dependence_signal`
- `same_family_fit_gain_signal`
- `uncured_only_support_signal`
- `triangulated_heterogeneity_bucket`
- `future_cure_modeling_recommendation_flag`

### 6.10.3 Four-axis interpretation matrix and derived bucket

Stage 10 must assemble four interpretation axes:

- `axis_timing_difference`:
  - `early_supported_difference`
  - `weak_or_absent_early_difference`
  - `late_only_difference`
- `axis_followup_support`:
  - `directly_observed_data_supported`
  - `partly_model_dependent`
  - `mostly_extrapolated`
- `axis_cure_support`:
  - `supportive`
  - `equivocal`
  - `unsupportive`
- `axis_remission_distortion`:
  - `low`
  - `moderate`
  - `high`

A derived `interpretation_bucket` must then use one of:

- `timing_difference_dominant`
- `remission_distorted_apparent_cure`
- `cure_like_heterogeneity_plausible`
- `tail_only_projection_dominant`
- `mixed_or_indeterminate`

Preferred reading rules are:

- `timing_difference_dominant` when meaningful site differences are already present in the directly observed early windows and cure-specific support is not the dominant explanation.
- `remission_distorted_apparent_cure` when remission-aware branches materially attenuate or otherwise rewrite the apparent cure-like signal or false-positive burden conclusions.
- `cure_like_heterogeneity_plausible` when Stage 6 is supportive or at least not contradicted, supported horizons still show clinically meaningful differences after cure modeling, and remission-aware reanalysis does not erase the main pattern.
- `tail_only_projection_dominant` when divergence is driven mainly by projection-dominant or prior-tail-sensitive horizons.
- `mixed_or_indeterminate` when no single mechanism clearly dominates.


### 6.10.4 Required figures and supporting display objects

The following display objects are mandatory:

- `horizon_support_panel`
  - must show horizon, risk-set support, `support_tier`, `horizon_evidence_class`, and `claim_restriction_flag`
  - must visibly distinguish directly observed horizons from model-dependent and extrapolated horizons
  - must overlay `prior_tail_sensitive` where applicable

- `8A_vs_8B_delta_panel`
  - must be computed and exported by Stage 8B using Stage 8A as the reference branch
  - must show the remission-aware delta quantities by dataset and horizon
  - its purpose is to prevent readers from treating the competing-risk extension as invisible

- `anchor_vs_neutral_delta_panel`
  - must show how supported-horizon risk, false-positive burden, and clinical-usefulness summaries change when the external incidence anchor is removed
  - its purpose is to make the value of the external meta-analytic anchor visible rather than implicit

- `incidence_anchor_update_panel`
  - must show the externally anchored incidence pattern alongside the posterior-updated incidence pattern
  - it may be displayed on the logit or one-year-risk scale, but the chosen scale must be explicit
  - its purpose is to show how the initial anchor changes after estimation rather than letting the anchor disappear into the posterior

- one cure-model-only supporting decomposition table or figure
  - this must remain separate from the main common-scale comparison table
  - it is the correct place for incidence–latency decomposition, `susceptible_fraction`, uncured-only annual summaries, optional `MSTu`, and merged site-placement interpretation

- `heterogeneity_triangulation_panel`
  - must display the synthesized stage-level signals and the resulting `triangulated_heterogeneity_bucket`
  - may be rendered as a structured table or figure, but it must remain separate from the main common-scale comparison table
  - its purpose is to keep the indirect heterogeneity argument explicit, auditable, and distinct from the primary risk-comparison ruler

### 6.10.5 Reference basis

- this integrated master specification
- `R13. McLernon DJ, Giardiello D, Van Calster B, et al. Assessing performance and clinical usefulness in prediction models with survival outcomes. Ann Intern Med. 2023.pdf`
- `R11. Vickers AJ, Elkin EB. Decision curve analysis- A novel method for evaluating prediction models. Med Decis Making. 2006..pdf`
- `R12. Vickers AJ, Cronin AM, Elkin EB, Gonen M. Extensions to decision curve analysis, a novel method for evaluating diagnostic tests, prediction models and molecular markers. BMC Med Inform Decis Mak. 2008.pdf`
- `2014_Systematic review and collaborative reca_Van Der Werf et al..pdf`
- `2024_A two-sample comparison of mean survival_Dobler and Musta.pdf`
- `R21. BMJ 2023 074819. Full text part 1. (To verify).pdf`
- `R22. BMJ 2023 074820. Full text part 2. (To verify).pdf`


### 6.10.6 What later R implementation must preserve

- one long-format prediction table as source of truth
- one long-format classification / performance table
- one reporting metadata table with `support_tier`, `horizon_evidence_class`, `claim_restriction_flag`, `prior_tail_sensitive`, `prior_branch`, and `site_prior_family` where applicable
- one cure-model-only supporting decomposition table
- one remission-aware delta table
- one external-anchor value table
- one heterogeneity-triangulation table
- direct comparability across all model classes **within each risk scale**
- explicit supplementary treatment of remission-sensitive competing-risk outputs
- figure-ready source tables for the `horizon_support_panel`, `8A_vs_8B_delta_panel` (exported by Stage 8B), `anchor_vs_neutral_delta_panel`, `incidence_anchor_update_panel`, and `heterogeneity_triangulation_panel`

## 7. Operational workflow and parallelization rule

This section is the integrated replacement for the former standalone workflow guide and for the operational execution role previously split across the framework and breakdown documents.

The clean operational workflow is:

> **Stage 1 → (Stage 2, 3, 4, 5, 6 in parallel) → (Stage 7, 8, 9 code preparation can be parallel) → Stage 10**

### 7.1 Canonical workflow role now contained in this file

This file now carries three workflow roles internally:

- the **canonical staged framework** role formerly associated with `5.💛Model Specification Framework_🇬🇧ENG.md`
- the **operational breakdown** role formerly associated with `4.Modeling Framework Breakdown_🇬🇧ENG.md`
- the **round-based execution guide** role formerly associated with `staged_workflow_guideline.md`

No separate shared non-Bayesian workflow source is required once this section and the stage-specific sections below are adopted.

### 7.2 Round structure

Recommended operational rounds are:

- **Round A:** Stage 1 only
- **Round B:** Stage 2, 3, 4, 5, and 6 in parallel after Stage 1
- **Round C:** Stage 7, 8, and 9 code preparation and execution after the Round B backbone outputs exist, with Stage 6 carry-forward available for interpretation
- **Round D:** Stage 10 always last

### 7.3 What can be parallelized

After Stage 1:

- Stage 2, 3, 4, 5, and 6 can be prepared and run in parallel
- Stage 7, 8, and 9 code can also be prepared in advance

### 7.4 What is not fully independent

Interpretation is still staged.

In particular:

- Stage 7 should carry forward Stage 6 screening results
- Stage 8 should carry forward Stage 6 screening results
- Stage 8 must not repeat Stage 6 screening
- Stage 8A and Stage 8B are not substitutes for each other
- Stage 7 and Stage 8 should be read in light of Stage 2 and Stage 4 as well, not only Stage 6
- Stage 9 should be interpreted on the supplementary competing-risk scale
- Stage 10 is the only stage allowed to answer the final overall scientific question

### 7.5 Practical code-design rule

A practical implementation pattern is:

- Stage 7 scripts accept Stage 6 screening flag CSV as an optional input
- Stage 8 scripts accept Stage 6 screening flag CSV as an optional input
- Stage 7 / 8 code may be prepared before Stage 6 finishes, but final interpretation must still use the carried-forward Stage 6 flags
- Stage 8 outputs must include explicit `branch` and `risk_scale` fields so they can join safely with later tables

### 7.6 Stage-prompt architecture rule

When requesting stage-specific code or stage-specific analyses, the shared assumptions do **not** need to be restated from external framework documents.

Instead, this integrated file should be treated as containing:

- the shared backbone,
- the stage purpose,
- the required quantities,
- the required exports,
- the carry-forward rules,
- and the workflow ordering constraints.

Separate prompts may still refer to `Rules Before Generating R Code_🇬🇧ENG.md`, `3.Data Dictionary_🇬🇧ENG.md`, and the Bayesian Stage 8 companion, but not to now-absorbed framework/workflow files.

### 7.7 Deletion implication

Because the canonical framework role, operational breakdown role, and round-based workflow role are now all recorded inside this document, the former shared non-Bayesian framework/workflow documents may be archived or removed after project adoption of this file.

---

## 8. Recommended reporting hierarchy for the final manuscript

### 8.1 Primary findings

Use these as the main scientific claims:

- **PNU:** 1-year horizon only as the safest primary result
- **SNU:** 1-year and 2-year horizons as primary
- **merged:** 1-year and 2-year horizons as primary, with both site-free and site-adjusted versions shown
- main target: false-positive burden and clinical usefulness
- main event: transition only
- remission sensitivity shown alongside, not hidden

### 8.2 Secondary findings

Use these as important supporting results:

- 5-year results for SNU and merged, when they remain only partly model dependent rather than projection-dominant
- cohort-level risk shift when cure-like heterogeneity is allowed
- anchor-informed versus neutral / no-external-information Bayesian comparison on supported horizons
- the prior-to-posterior incidence-shape update panel
- discrimination, calibration, and Brier / IBS summaries
- timing pattern of the site effect
- separate incidence and latency interpretation
- the cure-model-only supporting decomposition block
- the heterogeneity-triangulation block
- the future cure-modeling recommendation flag

### 8.3 Sensitivity findings

Use these as important but more tentative results:

- 10-year projections
- any row labeled `projection_only` or `projection_plus_prior_sensitive`
- subgroup KM benchmarks
- remission competing-risk / multi-state analyses
- merged site-free versus site-adjusted contrasts
- flexible parametric AFT-cure or other extension-path sensitivities
- site-prior family sensitivities beyond the main retained Bayesian fits
- any extra extrapolation-sensitivity summaries

---

## 9. Stage-specific reporting and QC rules (adapted into English)

This section incorporates the useful reporting and QC logic from the stage-report prompt, but converts it into a reusable English reporting standard.

### 9.1 First rule: detect the stage before interpretation

Before writing interpretation:

- identify which stage the files belong to
- cluster mixed files by stage
- do not collapse multiple stages into one undifferentiated conclusion

### 9.2 Source-of-truth priority

When available, use this source priority:

1. export manifest
2. metadata registry
3. validation summary / validation issue rows / quality flags
4. machine-readable CSV tables
5. RDS bundle summaries
6. PDF plots

CSV or other machine-readable tables are the source of truth for exact numbers.

Plots are for visual corroboration, not for inventing precise numbers.

### 9.3 Mandatory QC before interpretation

QC must be performed before narrative interpretation.

Check:

- file completeness
- stage coherence
- obvious required-file omissions
- duplicate layered exports versus true duplication
- cross-file consistency in labels, model names, families, horizons, thresholds, counts, and support flags
- backbone consistency
- validation summaries and quality flags
- fitted-model diagnostics where relevant
- plot-versus-table consistency
- support strength at 1, 2, 5, and 10 years

### 9.4 Do-not-misclassify rules

Do **not** automatically call the following computational failures:

- Stage 2 sparse-tail or masking flags
- Stage 3 all-high-risk or all-low-risk KM patterns
- Stage 4 intentionally suppressed or bootstrap-unstable late-tail rows
- Stage 6 mixed evidence or equivocal screening decisions
- Stage 6 inherited raw-tail metrics for `merged__site_adjusted`
- Stage 8 retaining multiple admissible Bayesian fits rather than forcing one winner

Call them failures only if they contradict the stated stage logic or create impossible outputs.

### 9.5 Major-failure versus no-major-failure branching

#### If a major computational or structural failure exists

The report should clearly state:

- what is broken
- whether the defect is fatal or limited
- which claims become insecure
- whether repair requires:
  - re-export,
  - rerun current stage,
  - rerun from earlier stage,
  - partial code revision,
  - or only in extreme cases, major rewrite

#### If no major failure exists

The report should clearly state:

- no major computational/data-structure failure was identified
- nonfatal technical cautions
- interpretive cautions
- and only then proceed to full interpretation

### 9.6 Evidence-first reporting discipline

Do not jump directly from a number to interpretation.

Before each major interpretation, use an evidence block with the following structure:

#### Evidence block [number]. [title]

- Source file(s):
- Rows / filters used:
- Key columns used:
- Why this table matters:
- Claim-relevant complete subtable:

Then explain in this order:

1. what the numbers directly say
2. simple-language explanation
3. why it matters for CHR-P
4. what can be claimed now
5. what cannot be claimed now

### 9.7 CHR-P interpretation lens

Every stage report should explicitly ask:

- How would this change CHR-P transition-risk communication?
- Does it reduce or increase unnecessary high-risk labeling?
- Does it alter monitoring intensity or follow-up duration?
- Could PNU vs SNU differences reflect timing, selection, care-pathway, or follow-up structure rather than biology?
- Could remission handling distort apparent long-term transition risk?
- Could an apparent cure-like signal instead be follow-up immaturity or tail instability?


### 9.8 Interpretation rules for claims

Do not reduce interpretation to:

- “model X won”
- “lowest AIC”
- “statistically superior”

Instead interpret:

- disease process
- timing difference
- follow-up maturity
- remission handling
- latent heterogeneity
- external-prior dependence
- cohort/site context
- care pathway proxy effects
- assumption plausibility

### 9.8.1 Allowed claims

Acceptable claims include statements such as:

- at directly observed or otherwise supported horizons, allowing cure-like heterogeneity changed estimated transition risk and altered false-positive burden or clinical usefulness by a stated amount
- remission-aware competing-risk analysis attenuated, amplified, or otherwise changed the apparent cure-like signal by a stated amount
- anchor-informed versus neutral / no-external-information Bayesian comparison did or did not materially change supported-horizon conclusions
- the posterior-updated incidence pattern remained close to, or moved away from, the external anchor in a stated way
- merged site-placement decomposition suggests that cohort difference is operating mainly through incidence, latency, both, or neither
- late-horizon results are better treated as partly model dependent or projection-dominant rather than directly observed facts
- the observed CHR-P data are more compatible with clinically meaningful heterogeneity than with a single homogeneous no-cure survival process
- future CHR-P studies should consider cure-aware modeling if supported-horizon risk and false-positive burden change materially relative to no-cure comparators
- future CHR-P survival studies should consider mixture-cure-type modeling when supported-horizon shifts and the heterogeneity-triangulation block jointly support that recommendation

### 9.8.2 Forbidden claims

The following are forbidden:

- “the Bayesian cure model proved cure”
- “the present study proved a discrete cured subgroup within CHR-P”
- “the present study identified the true latent subgroup structure of CHR-P”
- “the external meta-analysis prior proved the true cure fraction”
- treating `projection_only` or `mostly_extrapolated` horizons as observed facts
- using `projection_plus_prior_sensitive` horizons for primary or headline claims
- treating posterior stabilization alone as proof of a discrete cured subgroup
- treating a subject-level posterior cure probability as a stand-alone patient-level decision number in the current project
- writing a merged site coefficient or site-placement term as a causal treatment effect
- presenting latent-class or frailty interpretations as if they were established by the current baseline workflow
- presenting a cross-scale Stage 8A-versus-Stage 8B delta as if it were a fair within-scale model-ranking contest
- presenting the anchor-informed branch as automatically correct simply because it uses a larger external data source

### 9.9 Natural-frequency rule

When horizons are sufficiently supported, translate important absolute risks into simple language such as:

- “about X out of 100”

If the horizon is weakly supported, explicitly say it is better treated as a model-based projection.

### 9.10 Carry-forward rule for stage reports

Every stage report should explicitly state:

- what this stage securely contributes to later stages
- what later stages still need to resolve
- what should be carried forward into Stage 10
- what must not yet be treated as final

---

## 10. Appendix: implementation-facing notes

### 10.1 Support-tier and evidence-dependence standardization

A machine-readable support field should use:

- `primary_supported`
- `secondary`
- `sensitivity`
- `projection`

A parallel machine-readable evidence-dependence field should use:

- `directly_observed_data_supported`
- `partly_model_dependent`
- `mostly_extrapolated`

A machine-readable reporting-governance field should use:

- `primary_claim_allowed`
- `secondary_or_sensitivity_only`
- `projection_only`
- `projection_plus_prior_sensitive`

Human-readable labels may render these as `primary-supported`, `secondary`, `sensitivity`, `projection`, `directly observed`, `partly model dependent`, `mostly extrapolated`, and `projection + prior-sensitive`.

The standardized meaning must remain unchanged even if display labels vary.

### 10.2 Preferred Stage 2 reporting rule

For merged overall and merged sex completeness outputs:

- `site_adjusted` is the preferred reporting view
- `site_free` is retained only as sensitivity
- merged site-stratified summaries use `site_adjusted`
- merged site panels should inherit support tiers from the matched single-site cohort when available

### 10.3 Preferred Stage 6 carry-forward mapping for Stage 7 / Stage 8

Preferred upstream Stage 6 fields are those listed in Section 6.6.8.

In this revised specification, the following paired fields are part of the **canonical cross-stage export**:

- `primary_gate_flag`
- `receus_primary_class`
- `presence_modifier_flag`
- `cure_presence_support_flag`
- `followup_contradiction_flag`
- `followup_not_contradicted_flag`

The decomposed fields are retained deliberately because later Stage 7 / Stage 8 logic uses them directly.

Legacy convenience aliases may still be derived where needed:

- `stage6_final_class = cure_model_eligibility_flag`
- `carry_forward_stage8 = cure_model_eligibility_flag != "unsupportive"`

but the preferred cross-stage contract is now the canonical field set in Section 6.6.8, not the legacy aliases.

### 10.4 Frequentist cure incidence-link default

For the incidence part of frequentist cure models, default to:

- **logit link**

Alternative links (for example complementary log-log or probit) may be treated as sensitivity analyses rather than defaults unless a strong substantive reason exists.

### 10.5 Tail identifiability and flexible AFT note

In semiparametric cure models, tail identifiability can be fragile.  
A practical tail restriction such as a zero-tail completion may be needed after the last uncensored time to avoid identifiability problems.

This is especially relevant when comparing PNU and SNU, because follow-up maturity differs strongly by cohort.

Recent AFT-cure and flexible parametric AFT-cure developments are recognized as methodologically relevant. In this revised specification, standard AFT-latency cure is a required structural comparator in Stage 7, while more flexible parametric AFT-cure forms remain extension paths unless they are implemented in a way that preserves the shared Stage 1 backbone and Stage 10 comparison ruler.


### 10.6 Auxiliary uncured-only / non-cure-fraction decomposition note

Recent cure-model comparison work supports the idea that, in the presence of a cure fraction, it can be informative to compare the cure fractions together with the survival of the uncured rather than relying only on overall survival.

In this project, that literature justifies a **cure-model-only supporting decomposition block** that may include:

- `susceptible_fraction = 1 - cure_fraction`
- uncured-only / susceptible-only survival and risk on the annual 1–10-year grid
- optional mean survival time of the uncured (`MSTu`) when the uncured support is sufficiently identified and the result is not contradicted by follow-up screening

However, these estimands remain **secondary/supporting** because they condition on latent uncured status and therefore should not replace the main common-scale comparison ruler centered on horizon-specific transition risk, false-positive burden, and clinical usefulness.

If a full `MSTu` would require aggressive tail extrapolation or is not supportable from the observed follow-up, the workflow should keep uncured-only reporting on the annual grid and avoid forcing a single-number extrapolated summary.

### 10.7 Out-of-scope but acknowledged extensions

The following literatures are recognized as informative but are **not** adopted as baseline comparators in the current staged workflow:

- partially observed cure-status methods, because the current project does not directly observe verified transition-nonsusceptibility labels
- latent-class survival heterogeneity models, because they would change the interpretation of the current shared comparator set
- frailty-based cure extensions, because they add another heterogeneity layer beyond the current primary covariate/site structure

These may motivate future sensitivity analyses, but should not silently replace the current baseline architecture.

#### 10.7A Higher-order heterogeneity-model boundary

Given the current cohort sizes and rapidly shrinking late risk sets, higher-order heterogeneity models such as latent-class survival models, mixture-of-experts survival models, and frailty-based cure extensions are acknowledged as future methodological directions only.

They are not baseline comparators in the current project and must not replace the current KM / no-cure / frequentist cure / Bayesian cure architecture.

Their use would require substantially larger cohorts, stronger tail support, and a more stable identification basis than is currently available in the present data.

#### 10.7B Frailty-specific caution

Frailty is a useful conceptual language for unobserved heterogeneity, but in univariate survival data it can be difficult to distinguish from time-varying covariate effects or other departures from proportional hazards.

Accordingly, univariate frailty fits in the present project should be regarded as speculative background explanations rather than baseline workflow components.

Frailty-based models may therefore be discussed as future directions or limited methodological curiosities, but they must not be presented as established explanations of within-CHR-P heterogeneity in the current data.

### 10.8 Integrated replacement rule for former shared workflow documents

The following documents are now treated as **fully absorbed** into this file for shared non-Bayesian project use:

- `4.Modeling Framework Breakdown_🇬🇧ENG.md`
- `5.💛Model Specification Framework_🇬🇧ENG.md`
- `staged_workflow_guideline.md`

No separate shared non-Bayesian workflow source is required once this revised file is adopted.


## 11. One-paragraph integrated summary

> **This project uses one common horizon / threshold comparison ruler across PNU, SNU, and merged analyses, but direct model comparison is locked to a declared risk scale. The main scale is `transition_only_main`, on which transition is the event of interest and remission is treated as censoring; the supplementary remission-sensitive scale is `transition_cif_competing`, on which remission is modeled explicitly as a competing event. The workflow is staged: first lock the backbone, then describe follow-up maturity, operationalize KM as a formal benchmark, separate timing difference from cure interpretation using the interval-based `0–1`, `1–2`, `2–5` main structure with `>5` as tail diagnostic, fit the non-cure individualized comparator, screen cure appropriateness using frequentist tools only, fit the frequentist cure block with both Cox-latency and AFT-latency structural comparators, fit the Bayesian cure block as a separate stabilization and uncertainty layer without repeating Stage 6 screening, perform remission-sensitive reanalysis, and answer the final scientific question only in the unified comparison stage. The project’s aim is not to prove cure as an observed fact, but to test whether allowing cure-like heterogeneity changes supported-horizon risk estimates, false-positive burden, and clinical usefulness enough to matter, whether externally informed Bayesian modeling adds value over a no-external-information Bayesian fit, and whether the total triangulated evidence is more compatible with clinically meaningful within-CHR-P heterogeneity than with a single homogeneous no-cure survival process. Stage 6 is a frequentist cure-model eligibility gate centered on RECeUS-AIC, and its outputs are carried forward into Stages 7 and 8 using canonical semantic and decomposed flags rather than being overwritten. Stage 7 additionally exports family-matched cure-versus-non-cure fit-gain diagnostics, but naive chi-square LRT interpretation is not authorized as the default inferential rule because those contrasts remain supportive rather than decisive in the current workflow. Throughout the workflow, timing difference, follow-up immaturity, remission handling, external-prior dependence, and cure-like heterogeneity must remain analytically distinct; every horizon-specific row must carry support-tier, evidence-dependence, and claim-restriction fields; merged overall and sex follow-up completeness must prefer the site-adjusted structure; Stage 8A remains the main Bayesian branch; Stage 8B and Stage 9 remain remission-sensitive supplementary branches; the anchor-informed versus neutral Bayesian comparison and the prior-to-posterior incidence-shape update are mandatory; cure-model decomposition and the heterogeneity-triangulation block must remain separate from the main common-scale comparison table; higher-order latent-class or frailty explanations remain future directions only; and projection-dominant or prior-sensitive tails must be reported explicitly without being overclaimed as observed fact.**
