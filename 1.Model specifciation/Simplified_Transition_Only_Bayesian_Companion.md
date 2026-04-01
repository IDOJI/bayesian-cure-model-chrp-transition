# Simplified Bayesian Companion Specification for the Transition-Only Cure Modeling Project

## 1. Role of this companion document

This document is the Bayesian companion for the simplified transition-only cure modeling project.

It is intentionally narrower than the previous Stage 8 specification. It keeps the **Stage 8A-type Bayesian mixture cure model** as the only Bayesian branch and removes:

- the Stage 8B remission-aware competing-risk extension
- all sensitivity-analysis branches
- the previous multi-branch Bayesian comparison structure

The purpose of this companion is to define a practical Bayesian cure-model block that can be run directly inside the simplified project while preserving:

- the existing data backbone
- the transition-only main risk scale
- prior rationale based on external incidence evidence
- posterior inference
- and mandatory convergence review using trace plots

---

## 2. Fixed backbone inherited unchanged

This Bayesian companion inherits the shared project backbone without changing variable definitions.

### 2.1 Datasets

Fit and report the Bayesian model for:

- `PNU`
- `SNU`
- `merged`
- `merged (site-adjusted)`

### 2.2 Identifier and time variables

Use:

- `site + id` as the person identifier
- `days_followup` as the main time variable
- \[
  \text{time\_year} = \frac{\text{days\_followup}}{365.25}
  \]
  for reporting and prediction

### 2.3 Outcome coding retained from the data dictionary

Observed outcome coding remains:

- `status_num == 1` → transition
- `status_num == 2` → remission
- `status_num == 0` → right censoring

### 2.4 Simplified Bayesian event rule

This companion uses only the transition-only main scale:

**Label:** `transition_only_main`

- event: `status_num == 1`
- censoring: `status_num %in% c(0, 2)`

Therefore, remission is treated as censoring in this Bayesian block.

---

## 3. Purpose of the simplified Bayesian model

The Bayesian model serves four functions:

1. stabilize cure-model estimation when frequentist fits are numerically fragile  
2. propagate uncertainty more transparently than a single plug-in estimate  
3. allow external age/sex incidence evidence to inform the **incidence** component  
4. produce posterior horizon-specific risk summaries directly comparable to the frequentist outputs  

This model is not intended to identify cure by itself. It is a stabilization and uncertainty layer inside the broader workflow.

---

## 4. Model definition

## 4.1 Latent susceptibility structure

Let \(U_i\) denote latent transition susceptibility:

- \(U_i=1\): susceptible / uncured
- \(U_i=0\): non-susceptible / cured

Define:

\[
\pi_i = P(U_i=1)
\]

### 4.1.1 Plain-language interpretation

The model assumes the study population is a mixture of:

- subjects who are effectively not susceptible to transition
- subjects who are susceptible and, if followed long enough, can transition according to a latency distribution

Because remission is censored here, this is a **transition-only cure model**, not a remission-aware competing-risk model.

---

## 4.2 Incidence component

Model the susceptible probability with logistic regression:

\[
\operatorname{logit}(\pi_i)=\alpha_0+\alpha_a\,age_{s,i}+\alpha_m\,sex_{num,i}
\]

For `merged (site-adjusted)`:

\[
\operatorname{logit}(\pi_i)=\alpha_0+\alpha_a\,age_{s,i}+\alpha_m\,sex_{num,i}+\alpha_s\,site_i
\]

The default covariates for the incidence component are:

- `age_s`
- `sex_num`
- `site` for merged site-adjusted analysis

---

## 4.3 Latency component

Use a **log-normal AFT latency model** only.

For susceptible subjects:

\[
\log T_i \mid (U_i=1)=\eta_i+\sigma\varepsilon_i,\qquad \varepsilon_i \sim N(0,1)
\]

with

\[
\eta_i=\beta_0+\beta_a\,age_{s,i}+\beta_m\,sex_{num,i}
\]

and, for `merged (site-adjusted)`:

\[
\eta_i=\beta_0+\beta_a\,age_{s,i}+\beta_m\,sex_{num,i}+\beta_s\,site_i
\]

The susceptible-only survival is:

\[
S_{u,i}(t)=1-\Phi\left(\frac{\log t-\eta_i}{\sigma}\right)
\]

and the susceptible-only cumulative risk is:

\[
R_{u,i}(t)=1-S_{u,i}(t)
\]

---

## 4.4 Overall survival and cumulative risk

The overall survival is:

\[
S_i(t)=(1-\pi_i)+\pi_i S_{u,i}(t)
\]

The overall cumulative transition risk is:

\[
R_i(t)=1-S_i(t)=\pi_i\{1-S_{u,i}(t)\}
\]

### 4.4.1 Population-standardized reporting

For dataset-level reporting, posterior subject-level predictions should be averaged over the empirical covariate distribution of the corresponding analysis dataset:

\[
\bar S(t)=\frac{1}{n}\sum_{i=1}^{n} S_i(t), \qquad
\bar R(t)=1-\bar S(t)
\]

For susceptible-only summaries:

\[
\bar S_u(t)=\frac{1}{n}\sum_{i=1}^{n} S_{u,i}(t), \qquad
\bar R_u(t)=1-\bar S_u(t)
\]

These standardized summaries are what should populate the final Bayesian result tables.

---

## 4.5 Observed-data likelihood

Let:

\[
\delta_i = I(\text{status\_num}_i = 1)
\]

Then the observed-data likelihood contribution is:

\[
L_i = [\pi_i f_u(T_i)]^{\delta_i}\left[(1-\pi_i)+\pi_i S_u(T_i)\right]^{1-\delta_i}
\]

### 4.5.1 Plain-language reading of the likelihood

This says:

- if a subject transitions, they must belong to the susceptible subgroup and their event time follows the susceptible latency distribution
- if a subject is censored, they may either be truly non-susceptible or be susceptible but still event-free at censoring time

That is the essential statistical logic of the simplified Bayesian cure model.

---

## 5. Prior structure

## 5.1 General prior principle

External substantive information is allowed to enter only through the **incidence** component.

The latency component remains weakly informative and internally estimated from the project data.

This preserves the original project logic that the external meta-analytic information is used to guide who is likely to belong to the susceptible subgroup, not to force the timing pattern of transition among susceptible subjects.

---

## 5.2 Incidence priors informed by external meta-analytic evidence

Use the project’s external age/sex incidence evidence as a **prior-centering device** for the incidence coefficients.

### 5.2.1 Interpretation rule

The external evidence is used as:

- prior or calibration anchor
- starting structure for incidence heterogeneity
- not as fixed truth
- not as proof of the true cure fraction

The posterior must be allowed to move away from the prior center if the current data support it.

### 5.2.2 Practical prior statement

A practical implementation should:

- map the external age/sex incidence pattern to the current design matrix
- center the incidence priors on that mapped pattern
- use moderate variance so the current dataset can still update the prior

For the simplified project, no neutral-vs-anchor sensitivity branch is required. Only the anchored main prior is fit.

---

## 5.3 Default weakly informative priors for the remaining parameters

Recommended defaults:

### 5.3.1 Incidence component

- incidence intercept:  
  \[
  \alpha_0 \sim N(0, 2.5^2)
  \]
  unless replaced by an explicitly anchored center

- incidence non-site coefficients around the anchored center:
  \[
  \alpha_j \sim N(\mu_j^{anchor}, 1^2)
  \]

- site coefficient in merged site-adjusted incidence model:
  \[
  \alpha_s \sim N(0, 1^2)
  \]

### 5.3.2 Latency component

- latency intercept:
  \[
  \beta_0 \sim N(0, 2.5^2)
  \]

- latency non-site coefficients:
  \[
  \beta_j \sim N(0, 1^2)
  \]

- site coefficient in merged site-adjusted latency model:
  \[
  \beta_s \sim N(0, 1^2)
  \]

- log-normal scale prior:
  \[
  \log(\sigma) \sim N(0, 0.5^2)
  \]

These defaults are consistent with the project’s continuous-covariate scaling discipline and are intended to be weakly informative rather than strongly shrinkage-dominant.

---

## 6. MCMC estimation and convergence rules

## 6.1 Required sampler behavior

The Bayesian fit must use:

- multiple chains
- a warm-up / burn-in period
- posterior draws sufficient for stable horizon-specific posterior summaries

The exact number of iterations can be implementation-specific, but the diagnostics below are mandatory.

## 6.2 Mandatory convergence review

Every fitted Bayesian model must be checked with:

- trace plots
- \(\hat R\)
- effective sample size
- divergence / sampler pathology review where relevant

### 6.2.1 Trace plots are mandatory

Trace plots are not optional in this simplified project. They must be exported and reviewed for key parameters, including:

- incidence intercept
- major incidence coefficients
- latency intercept
- major latency coefficients
- log-normal scale parameter
- representative posterior cure fraction
- representative posterior 1-year and 2-year risk summaries

### 6.2.2 Minimum admissibility standard

A Bayesian fit should not be treated as admissible if any of the following remain unacceptable:

- severe chain non-mixing
- clearly unstable trace plots
- materially inflated \(\hat R\)
- inadequate effective sample size
- obvious divergence or degeneracy

---

## 7. Mandatory Bayesian outputs

## 7.1 Posterior summary table

For every fitted model, export a posterior summary table including:

- incidence coefficients
- latency coefficients
- log-normal scale parameter
- posterior cure fraction / susceptible fraction summaries
- posterior mean, SD, median, 95% CrI

---

## 7.2 Horizon-specific prediction tables

### 7.2.1 Overall-population table contribution

For each dataset view and each horizon \(t=1,\dots,10\) years, export:

- posterior mean overall survival
- posterior lower / upper 95% CrI overall survival
- posterior mean overall cumulative risk
- posterior lower / upper 95% CrI overall cumulative risk
- `support_tier`
- `horizon_evidence_class`
- `claim_restriction_flag`

These values populate the Bayesian columns of Table 1.

### 7.2.2 Susceptible-only table contribution

For each dataset view and each horizon \(t=1,\dots,10\) years, export:

- posterior mean susceptible-only survival
- posterior lower / upper 95% CrI susceptible-only survival
- posterior mean susceptible-only cumulative risk
- posterior lower / upper 95% CrI susceptible-only cumulative risk
- `support_tier`
- `horizon_evidence_class`
- `claim_restriction_flag`

These values populate the Bayesian columns of Table 2.

---

## 7.3 Trace plot output

Export a trace-plot panel for the key parameters listed above.

This should be treated as a required model diagnostic figure and not as supplementary-only decoration.

---

## 7.4 Prior rationale note

Every Bayesian fit should be accompanied by a short structured note documenting:

- which external incidence source was used
- how it was mapped into the incidence prior center
- why that prior is treated as a centering device rather than fixed truth

This note may be short, but it must be explicit.

---

## 8. Relationship to the non-Bayesian workflow

The Bayesian block plugs into the same simplified reporting framework as the non-Bayesian models.

### 8.1 Shared comparison framework

The Bayesian results must remain directly comparable to:

- Kaplan–Meier
- frequentist non-cure log-normal
- frequentist mixture cure log-normal

### 8.2 Same horizon grid

Use the same 1–10 year grid for all Bayesian outputs.

### 8.3 Same late-horizon governance

The Bayesian branch does not escape the project’s follow-up maturity limits. Late-horizon posterior stability is still subject to:

- direct-support vs model-dependent distinction
- projection labeling
- claim restriction

A stable posterior tail is not, by itself, evidence that the tail is empirically well supported.

---

## 9. Interpretation rules specific to the simplified Bayesian branch

## 9.1 Allowed claims

Allowed statements include:

- the Bayesian cure model changed supported-horizon transition risk estimates by a stated amount relative to the non-cure model
- the Bayesian cure model implied a stated posterior susceptible fraction and corresponding susceptible-only trajectory
- the posterior incidence pattern was informed by, but not fixed to, the external age/sex incidence evidence

## 9.2 Claims that are not allowed

Do not claim that:

- the Bayesian model proves that a cured subgroup exists
- posterior stability alone proves cure
- the external prior proves the true cure fraction
- site coefficients are treatment effects
- remission is unimportant simply because it was censored in the simplified design

---

## 10. Simplification rule

This companion intentionally defines a narrower Bayesian project than the previous Stage 8 framework.

Specifically:

- only one Bayesian cure family is fit: log-normal mixture cure
- only one event rule is used: transition with remission censored
- only one Bayesian branch is retained: Stage 8A-type transition-only cure modeling
- no Stage 8B competing-risk remission-aware branch is fit
- no sensitivity analyses are fit
- no latency-family competition is fit inside the Bayesian block

This is deliberate. The point of the simplified Bayesian companion is to provide a stable, interpretable Bayesian counterpart to the frequentist log-normal mixture cure model, not to preserve the full earlier Bayesian model-selection universe.

---

## 11. Drafting basis for this simplified Bayesian companion

This simplified Bayesian companion was drafted by retaining the main structural logic of the existing Stage 8A Bayesian cure design while removing the Stage 8B remission-aware branch and all sensitivity-analysis branches. It also preserves the current data backbone and the project-wide support-tier governance.
