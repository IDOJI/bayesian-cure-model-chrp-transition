# Bayesian Modeling Specification — Simplified Transition-Only Cure Project

**Status:** simplified Bayesian companion specification for the transition-only sub-project  
**Purpose:** keep the scientific core of Stage 8A while removing remission-aware modeling and all sensitivity-analysis branches from the main simplified report  
**Scope:** Bayesian mixture cure model only

---

## 1. Role of this document

This file is the Bayesian companion for the simplified transition-only sub-project.

It keeps the main Stage 8A ideas:

- transition-only event definition,
- mixture cure structure,
- log-normal latency,
- external meta-analytic anchor,
- posterior inference,
- and trace-plot-based convergence checking.

It removes from the main simplified project:

- Stage 8B remission-aware content
- remission competing-risk likelihoods
- neutral/no-external-information reporting branch
- site-prior sensitivity branch
- all other sensitivity analyses

---

## 2. Observed data and estimand

For subject $i$, define:

$$
T_i = \max\left(\frac{\text{days\_followup}_i}{365.25}, 10^{-8}\right)
$$

and

$$
D_i = I(\text{status\_num}_i = 1)
$$

with:

- `status_num == 1` = transition event
- `status_num == 2` = remission, treated as censoring
- `status_num == 0` = right censoring

Therefore the Bayesian model targets the same main scale:

- `transition_only_main`

---

## 3. Bayesian mixture cure structure

Let:

- $U_i = 1$ mean susceptible (uncured)
- $\pi_i = P(U_i = 1 \mid x_i)$

### 3.1 Incidence part

Use a logistic incidence model:

$$
\operatorname{logit}(\pi_i) = \alpha_0 + x_i^\top \beta
$$

with the simplified covariate set:

- `age_s`
- `sex_num`
- and `site` only when the merged site-adjusted model is being fit

### 3.2 Latency part

Conditionally on susceptibility, use a log-normal latency model:

$$
\log(T_i \mid U_i = 1, x_i) \sim N(\mu_i, \sigma^2)
$$

with

$$
\mu_i = \gamma_0 + z_i^\top \gamma
$$

The simplified main site-adjusted merged model retains site in both:

- the incidence component
- and the latency component

### 3.3 Overall survival

The overall survival function is:

$$
S_i(t) = (1 - \pi_i) + \pi_i S_u(t \mid x_i)
$$

and overall risk is:

$$
R_i(t) = 1 - S_i(t)
$$

### 3.4 Susceptible-only survival

For the uncured subgroup:

$$
S_{u,i}(t) = S_u(t \mid x_i)
$$

and

$$
R_{u,i}(t) = 1 - S_u(t \mid x_i)
$$

These susceptible-only quantities must be reported separately from the overall-population quantities.

---

## 4. Prior rationale

### 4.1 Core principle

The simplified Bayesian model keeps the original project’s key idea that large external meta-analytic transition incidence information can stabilize the incidence side of the cure model.

### 4.2 Simplified main prior strategy

Retain only the **anchor-informed** main prior strategy.

Operationally, the external evidence is used as:

- a prior or calibration anchor for the incidence component,
- especially for age- and sex-structured transition incidence patterns,
- not as fixed truth,
- and not as a remission model.

### 4.3 Site-prior rule

For the merged site-adjusted model, retain only the main site prior family:

- `normal_0_1_main`

Do not include the Student-t site prior sensitivity branch in the main simplified report.

### 4.4 Excluded branches

The simplified main Bayesian report excludes:

- neutral/no-external-information branch
- site-prior sensitivity branch
- remission-aware extension
- all Stage 8B content

---

## 5. Posterior inference

Posterior inference should produce, on the yearly grid `1:10`:

- overall survival probability
- overall cumulative risk
- susceptible-only survival probability
- susceptible-only cumulative risk
- cure fraction / susceptible fraction where available

and a parameter-level posterior summary for the main retained models.

The posterior summary should include, at minimum:

- posterior mean
- posterior SD
- 2.5% quantile
- median
- 97.5% quantile

---

## 6. Convergence and diagnostics

### 6.1 Mandatory trace plot

Trace plots are mandatory for the simplified Bayesian project.

They should be used to assess:

- chain mixing
- lack of visible drift
- lack of persistent chain separation

### 6.2 Additional required checks

Also report:

- `Rhat`
- effective sample size
- divergences
- tree-depth problems when relevant

### 6.3 Diagnostic interpretation rule

If the fit is not diagnostically acceptable, the model must not be silently reported as the main Bayesian result.

---

## 7. Reporting outputs

The simplified Bayesian block should export:

1. posterior summary table  
2. trace plot  
3. overall-population yearly estimates  
4. susceptible-only yearly estimates  
5. short prior-rationale note explaining how the external anchor was used

---

## 8. Interpretation limits

The simplified Bayesian model is a stabilization layer, not proof of a biologically observed cured class.

Its role is to show whether an anchor-informed Bayesian mixture cure model changes:

- risk estimates,
- late-horizon behavior,
- and susceptible-only interpretation

relative to the lighter transition-only framework.

---

## 9. Locked exclusions

The simplified Bayesian project must not include:

- remission-aware likelihood terms
- remission CIF modeling
- Stage 8B export structure
- sensitivity-analysis panels
- multi-branch prior comparison as a main reporting target

This document defines the full Bayesian scope for the simplified transition-only cure sub-project.
