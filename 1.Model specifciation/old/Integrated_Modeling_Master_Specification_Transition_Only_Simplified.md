# Integrated Modeling Master Specification — Simplified Transition-Only Cure Project

**Status:** new simplified companion specification for the transition-only sub-project  
**Purpose:** define a lighter and more interpretable survival modeling framework that keeps the existing data dictionary and backbone, but removes remission-state modeling and unnecessary stage complexity  
**Companion retained:** `Bayesian_Modeling_Specification_Transition_Only_Simplified.md`

---

## 1. Project role

This document defines a **simplified transition-only cure modeling project** nested inside the current CHR-P repository.

It is intended to be:

- shorter than the current full master specification,
- operationally compatible with the same data backbone,
- focused on a smaller set of models,
- and directly usable for a cleaner comparative report.

This simplified project is a **sub-project**, not a silent rewrite of the full staged workflow.

---

## 2. Scientific target

The simplified question is:

> When remission is treated as censoring and the analysis is reduced to transition versus censoring only, how do risk estimates change across a minimal set of interpretable survival models?

The project therefore focuses on:

1. a transparent overall survival/risk comparison,
2. a clear separation between overall-population and susceptible-only estimands,
3. a side-by-side comparison of site-unadjusted and site-adjusted reporting,
4. and a Bayesian cure model that keeps the external-anchor concept without the remission-aware extension.

---

## 3. Locked event definition

The event definition is fixed as:

- `event = (status_num == 1)` = transition
- `status_num == 2` = remission, treated as censoring
- `status_num == 0` = right censoring

Therefore the operational risk scale remains:

- `transition_only_main`

This simplified project does **not** fit a remission-aware model and does **not** use a competing-risk transition-vs-remission estimand.

---

## 4. Shared data backbone

The project inherits the same data dictionary and backbone definitions as the current repository.

### 4.1 Unit and key

- one row = one subject
- person identifier = `site + id`

### 4.2 Time variable

Use:

$$
T_i = \frac{\text{days\_followup}_i}{365.25}
$$

for reporting on the yearly grid:

$$
t \in \{1, 2, \ldots, 10\}\ \text{years}
$$

### 4.3 Core covariates

The minimal covariate backbone is:

- `age_exact_entry`
- `sex_num`
- `site` when a merged site-adjusted analysis is requested

Retain the stored Stage 1 scaled age definition:

$$
\text{age\_s} = \frac{\text{age\_exact\_entry} - \operatorname{mean}(\text{age\_exact\_entry})}{2\operatorname{sd}(\text{age\_exact\_entry})}
$$

### 4.4 Dataset versions to report

Report the following four dataset versions:

- `PNU`
- `SNU`
- `merged`
- `merged (site-adjusted)`

Interpretation rule:

- `PNU` and `SNU` are single-site cohorts, so site adjustment is not applicable there.
- `merged` means merged without site control.
- `merged (site-adjusted)` means merged with an explicit site term retained in the selected model.

---

## 5. Retained model set

Only the following four analyses belong to the simplified project.

### 5.1 Kaplan-Meier survival estimation

For each dataset version, estimate:

$$
\widehat{S}_{KM}(t)
$$

and define cumulative risk as:

$$
\widehat{R}_{KM}(t) = 1 - \widehat{S}_{KM}(t)
$$

This is the nonparametric benchmark.

### 5.2 Non-cure parametric survival model

Retain only the **log-normal** no-cure model.

In AFT form:

$$
\log(T_i) = \eta_i + \sigma \varepsilon_i,\qquad \varepsilon_i \sim N(0,1)
$$

with

$$
\eta_i = \beta_0 + \beta_1 \text{age\_s}_i + \beta_2 \text{sex\_num}_i
$$

and, for the merged site-adjusted version,

$$
\eta_i = \beta_0 + \beta_1 \text{age\_s}_i + \beta_2 \text{sex\_num}_i + \beta_3 \text{site}_i
$$

The reported estimands are:

- overall survival probability
- overall cumulative risk

### 5.3 Frequentist mixture cure model

Retain only the **log-normal** mixture cure model.

Let:

- $U_i = 1$ indicate that subject $i$ is susceptible (uncured)
- $\pi_i = P(U_i = 1 \mid x_i)$

Then:

$$
S_i(t) = (1 - \pi_i) + \pi_i S_u(t \mid x_i)
$$

and

$$
R_i(t) = 1 - S_i(t)
$$

For the susceptible-only subgroup:

$$
S_{u,i}(t) = S_u(t \mid x_i), \qquad R_{u,i}(t) = 1 - S_u(t \mid x_i)
$$

The simplified project reports both:

- overall-population quantities
- susceptible-only quantities

For the merged site-adjusted cure model, the main simplified definition is:

- site retained in both the incidence and latency parts

This avoids having multiple competing definitions of “site-adjusted” in the simplified report.

### 5.4 Bayesian mixture cure model

Retain only the **Stage 8A-type transition-only Bayesian mixture cure model**, with:

- remission excluded from the likelihood as an event type,
- remission treated as censoring only,
- log-normal latency only,
- the external meta-analytic anchor concept retained,
- no Stage 8B content,
- and no sensitivity-analysis branch in the main simplified report.

The Bayesian details are governed by the separate simplified Bayesian specification.

---

## 6. Reporting structure

### 6.1 Table 1 — overall survival and risk comparison

For each dataset version and each retained model, report at:

- 1, 2, 3, ..., 10 years after cohort entry

the following:

- estimated survival probability
- estimated cumulative risk

This table is strictly an **overall population** table.

### 6.2 Table 2 — susceptible-only quantities

For:

- frequentist mixture cure model
- Bayesian mixture cure model

report at:

- 1, 2, 3, ..., 10 years after cohort entry

the following:

- susceptible-only survival probability
- susceptible-only cumulative risk

Also carry:

- cure fraction where available,
- susceptible fraction where available,
- and late-horizon support flags.

### 6.3 Separation rule

Overall-population quantities and susceptible-only quantities must never be merged into one unlabeled table.

The simplified project must keep them explicitly separate.

---

## 7. Plotting rules

### 7.1 Plot Set 1 — overall curves

Create:

1. overall survival curves  
2. overall cumulative risk curves

Compare on the same axis:

- Kaplan-Meier
- non-cure log-normal
- frequentist mixture cure log-normal
- Bayesian mixture cure log-normal

### 7.2 Plot Set 2 — susceptible-only curves

Create:

1. susceptible-only survival curves  
2. susceptible-only cumulative risk curves

Main rule:

- the main susceptible-only plot set should include **cure models only**, because KM and no-cure models do not define a cure-fraction decomposition.

If reference overlays are ever shown in a later extension, they must be explicitly labeled as overall-population references and not allowed to visually merge with the susceptible-only estimand.

### 7.3 Late-horizon rule

Late horizons must be flagged as extrapolation/projection whenever appropriate.

At minimum:

- projection-dominant regions should be visually shaded or clearly annotated,
- and the table rows should carry `support_tier`, `horizon_evidence_class`, and `claim_restriction_flag`.

---

## 8. Simplified workflow

The simplified sub-project should keep the workflow short.

### 8.1 Minimal operational structure

1. load locked upstream backbone-compatible data or upstream model outputs  
2. retain only the four target model families  
3. standardize reporting to the four dataset versions  
4. export Table 1, Table 2, and the two plot sets  
5. export Bayesian posterior summary and trace-plot material

### 8.2 What is intentionally removed

The simplified project does not include:

- remission-aware modeling,
- competing-risk remission extensions,
- sensitivity-analysis branches,
- broad latency-family sweeps,
- or unnecessary stage proliferation.

---

## 9. Interpretation guardrails

### 9.1 What this project can say

This project can compare how estimates change under a lighter transition-only cure-modeling frame.

### 9.2 What this project cannot say

This project cannot by itself prove the true existence of a cured subgroup as an observed biological fact.

It is still a comparative modeling framework.

### 9.3 Projection rule

Rows marked as projection-dominant should not be used for primary claims without explicit caveats.

---

## 10. Output contract

The simplified project should export, at minimum:

- Table 1 overall comparison
- Table 2 susceptible-only comparison
- overall plot set
- susceptible-only plot set
- Bayesian posterior summary
- Bayesian trace plot
- a short prior-rationale note
- a model-selection registry

This is the required output contract for the simplified transition-only cure sub-project.
