# Bayesian Modeling Specification — Stage 8A main branch and Stage 8B / Stage 8-CR extension (Revised v5)

**Status:** enriched companion specification for the current project (feedback-integrated v5)  
**Purpose:** to restate the Stage 8 Bayesian design in a more defensible and literature-anchored form, with explicit competing-risk formulas, plain-language interpretation, a stronger Stage 8 export contract, and a clearer separation between standard methodology and project-specific extensions.  
**Scope:** Stage 8 only. This document does not replace Stage 4 timing separation, Stage 6 cure-appropriateness screening, Stage 9 remission-sensitive frequentist analysis, or Stage 10 unified interpretation.

---

## 1. Why this enriched version exists

The current project already has the right high-level direction:

- **Stage 8A** is the main Bayesian transition-only cure branch.
- **Stage 8B / Stage 8-CR** is a remission-sensitive Bayesian competing-risk extension.
- The project now uses **one shared comparison ruler** but **two locked risk scales**.

What was still worth strengthening is not the basic architecture but the **defensibility of the modeling narrative**:

1. the formulas for the competing-risk likelihood should be written out explicitly;
2. each formula should be translated into ordinary language so that readers can see what it means scientifically;
3. the document should say clearly which parts are **standard cure-model or survival-model practice** and which parts are **project-specific design choices**;
4. the document should identify the literature basis for each major design decision;
5. the document should state the main identification and interpretation limits instead of letting readers assume more than the model can support.

That is the purpose of this enriched companion.

---

## 2. Relationship to the master framework and locked boundaries

This document is subordinate to the current integrated framework and inherits the project’s locked staged architecture.

### 2.1 What remains unchanged

The following remain unchanged in substance:

- datasets: **PNU**, **SNU**, **merged**;
- person identifier: **site + id**;
- main time variable: **days_followup** and reporting variable $$
\text{time\_year} = \frac{\text{days\_followup}}{365.25}
$$
- common horizon grid: $$1, 2, \ldots, 10 \text{ years}$$
- common threshold grid;
- Stage 6 carry-forward logic;
- the rule that **Stage 8A** is the main Bayesian branch;
- the rule that **Stage 8B** is an additional remission-sensitive extension rather than a silent rewrite of the backbone.

### 2.2 What Stage 8 must not do

Stage 8 must not:

- rerun Stage 6 screening;
- quietly redefine the project’s main event rule;
- change the common horizon grid;
- change the common threshold grid;
- collapse the difference between the transition-only estimand and the remission-sensitive competing-risk estimand;
- answer the final scientific question by itself.


### 2.2A What Stage 8 cannot identify by itself

Stage 8 is a stabilization and uncertainty layer, not a stand-alone identification layer.

It may regularize a weakly identified incidence–latency decomposition, but it cannot by itself determine whether an apparent cure-like signal is driven by:

- timing difference,
- follow-up immaturity or late-tail instability,
- remission handling,
- or true cure-like heterogeneity.

Those distinctions remain the joint task of Stages 2, 4, 9, and 10.

### 2.3 Two locked risk scales

The project now uses two locked risk scales.

#### Scale A

**Label:** `transition_only_main`

- event of interest: $\text{status\_num} = 1$
- censoring: $\text{status\_num} \in \{0, 2\}$
- remission is treated as censoring
- this is the main operational scale for Stage 8A

#### Scale B

**Label:** `transition_cif_competing`

- event of interest: $\text{status\_num} = 1$
- competing event: $\text{status\_num} = 2$
- right censoring: $\text{status\_num} = 0$
- this is the remission-sensitive competing-risk scale for Stage 8B

The two scales answer different questions and should never be mixed as if they were identical.

---

## 3. Shared backbone used by both Stage 8A and Stage 8B

Unless stated otherwise, both branches inherit the same backbone.

### 3.1 Datasets and identifiers

Analyze:

- `PNU`
- `SNU`
- `merged`

Use:

- `site + id` as the person identifier
- `days_followup` as the main time variable
- $$
\text{time\_year} = \frac{\text{days\_followup}}{365.25}
$$
for reporting and prediction

### 3.2 Observed-state coding

Observed outcome coding is fixed as:

- $\text{status\_num} = 1 \rightarrow \text{transition}$
- $\text{status\_num} = 2 \rightarrow \text{remission}$
- $\text{status\_num} = 0 \rightarrow \text{right censoring}$

### 3.3 Common horizons and thresholds

Both Stage 8A and Stage 8B must export on the same ruler:

- horizons: $$1, 2, \ldots, 10 \text{ years}$$
- the same pre-specified threshold grid as the rest of the workflow

### 3.4 Shared core covariates

The common covariate backbone is:

- `age_exact_entry`
- `sex_num`
- `site` where applicable in merged analyses

Use the stored Stage 1 scaled age definition:

$$
\text{age\_s} = \frac{\text{age\_exact\_entry} - \operatorname{mean}(\text{age\_exact\_entry})}{2 \cdot \operatorname{sd}(\text{age\_exact\_entry})}
$$

### 3.5 Stage 6 carry-forward fields

Both Stage 8 branches must carry forward, not recreate, the Stage 6 interpretation fields. At minimum use the standardized canonical fields:

- `cure_model_eligibility_flag`
- `primary_gate_method`
- `primary_gate_flag`
- `receus_primary_class`
- `presence_modifier_flag`
- `cure_presence_support_flag`
- `followup_contradiction_flag`
- `followup_not_contradicted_flag`
- `screening_note`

Optional legacy loader aliases may additionally exist for implementation convenience:

- `stage6_final_class = cure_model_eligibility_flag`
- `carry_forward_stage8 = cure_model_eligibility_flag != "unsupportive"`


### 3.6 Horizon-support governance inherited from the master specification

Every horizon-specific Stage 8 export row must also carry:

- `support_tier`
- `horizon_evidence_class`
- `claim_restriction_flag`

The standardized support-tier values are:

- `primary_supported`
- `secondary`
- `sensitivity`
- `projection`

The required evidence-dependence classes are:

- `directly_observed_data_supported`
- `partly_model_dependent`
- `mostly_extrapolated`

The required reporting-governance classes are:

- `primary_claim_allowed`
- `secondary_or_sensitivity_only`
- `projection_only`
- `projection_plus_prior_sensitive`

Preferred horizon mapping is inherited directly from the master specification:

- **PNU**
  - 1 year → `support_tier = primary_supported`, `directly_observed_data_supported`
  - 2 years → `support_tier = sensitivity`, `partly_model_dependent`
  - 3+ years → `support_tier = projection`, `mostly_extrapolated`

- **SNU**
  - 1–2 years → `support_tier = primary_supported`, `directly_observed_data_supported`
  - 3–5 years → `support_tier = secondary`, `partly_model_dependent`
  - 6+ years → `support_tier = projection`, `mostly_extrapolated`

- **merged**
  - 1–2 years → `support_tier = primary_supported`, `directly_observed_data_supported`
  - 3–5 years → `support_tier = secondary`, `partly_model_dependent`
  - 6+ years → `support_tier = projection`, `mostly_extrapolated`

Rows labeled `projection_only` or `projection_plus_prior_sensitive` may not support primary claims. This restriction is part of the Stage 8 export contract itself, not just a later writing suggestion.

## 4. Notation used throughout

For subject $i$, define:

- $$
T_i = \max\left(\frac{\text{days\_followup}_i}{365.25}, 1\mathrm{e}{-8}\right)
$$
- $$
D_i = \text{status\_num}_i
$$

Also define the merged-site indicator:

- $$
s_i = I(\text{site}_i = \text{"SNU"})
$$

Thus:

- $s_i = 0$ corresponds to PNU
- $s_i = 1$ corresponds to SNU

### 4.1 Incidence-part covariates

Let:

- $$
\text{age}_i = \text{age\_exact\_entry}_i
$$
- $$
z_i = \text{sex\_num}_i
$$

with project coding:

- $\text{sex\_num} = 0 \rightarrow \text{Female}$
- $\text{sex\_num} = 1 \rightarrow \text{Male}$

Define age-band indicators:

- $$
x_{20,i} = I(20 \le \text{age}_i < 30)
$$
- $$
x_{30,i} = I(\text{age}_i \ge 30)
$$

Reference stratum:

- Female, age < 20 years

Incidence design vector:

- $$
X_i^{inc} = \left( z_i, x_{20,i}, x_{30,i}, z_i x_{20,i}, z_i x_{30,i} \right)^{T}
$$

### 4.2 Latency-part covariates

Let:

- $$
a_i = \text{age\_s,i}
$$
- $$
az_i = a_i z_i
$$

This follows the project-wide rule of scaling continuous covariates on a **2-SD scale** while keeping binary predictors on the natural 0/1 scale.

---

## 5. Stage 8A — Main Bayesian transition-only cure branch

### 5.1 Purpose

Stage 8A is the main Bayesian cure branch. Its role is to:

- stabilize estimation when frequentist cure fits are numerically fragile,
- propagate uncertainty more transparently,
- allow a carefully specified external prior in the incidence component only,
- export posterior risk and threshold summaries on the same comparison ruler used later by Stage 10.

### 5.2 Estimand and event rule

Stage 8A keeps the current project backbone:

- event: $\text{status\_num} = 1$
- censoring: $\text{status\_num} \in \{0, 2\}$

Its risk scale is:

$$
\text{risk\_scale} = \text{transition\_only\_main}
$$

Its main estimand is the transition-only mixture-cure risk.

### 5.3 Mixture-cure decomposition

Let $U_i$ denote latent transition susceptibility:

- $U_i = 1$ means transition-susceptible
- $U_i = 0$ means transition non-susceptible

Define:

- $$
\pi_i = P(U_i = 1)
$$

Then the population survival is:

$$
S_{pop,i}(t) = (1 - \pi_i) + \pi_i S_u(t \mid i)
$$

and the transition-only cumulative risk is:

$$
R_i(t) = 1 - S_{pop,i}(t) = \pi_i \{1 - S_u(t \mid i)\}
$$

### 5.4 Plain-language interpretation

This means Stage 8A assumes that the population is a mixture of:

- people who are effectively **not susceptible to transition**, and
- people who are susceptible, but among whom transition happens with some latency distribution.

Because remission is absorbed into censoring on this branch, Stage 8A does **not** estimate the real-world probability of transition under explicit competition from remission. It estimates transition risk on the project’s main transition-only scale.


### 5.5 Incidence component

#### 5.5.1 External prior basis

Only the incidence component uses external substantive information. The current external basis is the age-band by sex incidence pattern derived from van der Werf–type schizophrenia incidence evidence.

This external information is used as **prior centering**, not as hard truth and not as a prior assertion about the true cure fraction.

#### 5.5.1A How the external incidence template is converted into a logistic prior center

Let $g$ index the external age-by-sex cells.

- Let $r_g^{(10k)}$ denote the reported annual incidence rate per 10 000 person-years.
- Convert to a person-year rate:
  $$
  r_g = \frac{r_g^{(10k)}}{10000}
  $$
- Convert that annual rate to a one-year risk:
  $$
  q_g = 1 - \exp(-r_g)
  $$
- Move to the logistic scale:
  $$
  y_g = \operatorname{logit}(q_g)
  $$

Take the reference cell $g_0 =$ Female, age $<20$ years. Then:

$$
\alpha_{GP}^{vdW} = y_{g_0}
$$

The current project-specific prior mean vector $\mu^{vdW}$ is a **contrast vector on the current incidence design matrix**, not a free-floating list of external numbers. It is defined so that for each non-reference age-by-sex cell represented by the current design matrix,

$$
X_g^{inc}\mu^{vdW} = y_g - y_{g_0}
$$

When the external template is richer than the current project design matrix, the collapsed contrast vector must be precomputed once and then frozen as a design artifact. It must not be re-derived ad hoc during posterior fitting.

The prior scale vector $s_0^{vdW}$ likewise belongs to the frozen external-centering artifact and must be stored alongside the mapping metadata used to collapse the external template into the current incidence design.

#### 5.5.2 Anchor-informed branch

Main incidence model:

$$
\operatorname{logit}(\pi_i) = \alpha_{GP}^{vdW} + \Delta + X_i^{inc}\beta + \beta_s^{inc} s_i \, I(\text{incidence-site-branch-on})
$$

with priors:

- $\Delta \sim \mathcal{N}(0, 3.5^2)$
- $\beta_j \sim \mathcal{N}(\mu_j^{vdW}, s_{0,j}^2)$ for the five age-band by sex coefficients
- when included in merged incidence-site branches: $\beta_s^{inc} \sim \mathcal{N}(0, 1^2)$ in the main site-prior family
- mandatory site-prior sensitivity for merged incidence-site branches: $\beta_s^{inc} \sim t_{3}(0, 1)$

Female <20 anchor:

$$
\alpha_{GP}^{vdW} = \operatorname{logit}\left(1 - \exp\left(-\frac{0.69}{10000}\right)\right) = -9.581369553169
$$

Current prior mean vector:

$$
\mu^{vdW} = \left(0.419871845822, 0.907608052926, 0.586202561451, 0.466865123863, 0.037997248763\right)^{T}
$$

Current prior scale vector:

$$
s_{0}^{vdW} = \left(0.132789397422, 0.173731076538, 0.191221553945, 0.270393197518, 0.302838606651\right)^{T}
$$

#### 5.5.3 Neutral no-external-information branch

A mandatory prior-sensitivity branch retains the same design matrix but removes the external pattern-centering of the incidence coefficients **and** the externally centered baseline intercept. This is the project’s required **no-external-information comparator** for the anchored Bayesian analysis.

$$
\operatorname{logit}(\pi_i) = \alpha_0 + X_i^{inc}\beta + \beta_s^{inc} s_i \, I(\text{incidence-site-branch-on})
$$

with priors:

- $\alpha_0 \sim \mathcal{N}(0, 3.5^2)$
- $\beta_j \sim \mathcal{N}(0, 2^2)$ for the five age-band by sex coefficients
- $\beta_s^{inc} \sim \mathcal{N}(0, 1^2)$ when included in the main site-prior family
- mandatory site-prior sensitivity for merged incidence-site branches: $\beta_s^{inc} \sim t_{3}(0, 1)$

Preferred exported label:

- `prior_branch = neutral_no_external_info`

Legacy implementation alias allowed:

- `neutral_weakly_informative`

### 5.6 Why two prior branches are required

This is necessary because the project is not only asking whether Bayesian stabilization is numerically useful, but also whether any cure-like signal is robust to the incidence prior.

The external meta-analysis is one of the core substantive motivations for the Bayesian block, so the comparison between:

- `anchor_informed`
- `neutral_no_external_info`

is not an optional side note. It is part of the main Stage 8 scientific question.

If the signal is strong only under the externally anchored incidence prior and weak under the no-external-information branch, the signal is more prior-dependent than data-driven. That comparison must therefore be exported explicitly and shown in Stage 10.

### 5.7 Latency component

The latency component should stay on the common backbone as closely as possible.

#### 5.7.1 PNU and SNU structural branches

- `L0`: $$
\mu_i^{lat} = \gamma_0 + \gamma_a a_i + \gamma_m z_i
$$
- `L1`: $$
\mu_i^{lat} = \gamma_0 + \gamma_a a_i + \gamma_m z_i + \gamma_{am} a_i z_i
$$

#### 5.7.2 Merged structural branches

Use:

- incidence-site switch: `I0` or `I1`
- latency-site / latency-interaction branches:
  - `L0S0`
  - `L1S0`
  - `L0S1`
  - `L1S1`

Merged latency submodels:

- `L0S0`: $$
\mu_i^{lat} = \gamma_0 + \gamma_a a_i + \gamma_m z_i
$$
- `L1S0`: $$
\mu_i^{lat} = \gamma_0 + \gamma_a a_i + \gamma_m z_i + \gamma_{am} a_i z_i
$$
- `L0S1`: $$
\mu_i^{lat} = \gamma_0 + \gamma_a a_i + \gamma_m z_i + \gamma_s s_i
$$
- `L1S1`: $$
\mu_i^{lat} = \gamma_0 + \gamma_a a_i + \gamma_m z_i + \gamma_{am} a_i z_i + \gamma_s s_i
$$

#### 5.7.3 Site-placement interpretation in merged analyses

Merged outputs should preserve four site-placement structures:

- `site_in_neither`
- `site_in_incidence_only`
- `site_in_latency_only`
- `site_in_both`

This is scientifically important because the merged analysis is trying to separate whether PNU/SNU difference is driven by:

- who is susceptible,
- when susceptible subjects transition,
- or both.

Any merged site coefficient or site-placement label is a structural decomposition term only. It must **not** be described as a clean clinical treatment effect, because site here stands in for care pathway, referral, selection, and follow-up maturity differences.

### 5.8 Latency families

The main Stage 8A family set remains:

- exponential
- Weibull
- log-normal
- log-logistic

These are retained because they cover the main parametric shapes already used in the project, permit direct comparison with Stage 7, and remain computationally feasible in a Bayesian workflow.


### 5.9 Priors for latency coefficients

Use weakly informative priors consistent with the Stage 1 scaling discipline:

- $\gamma_0 \sim \mathcal{N}(0, 2.5^2)$
- non-site latency regression coefficients $\gamma_j \sim \mathcal{N}(0, 1^2)$
- when included in merged latency-site branches: $\gamma_s \sim \mathcal{N}(0, 1^2)$ in the main site-prior family
- mandatory site-prior sensitivity for merged latency-site branches: $\gamma_s \sim t_{3}(0, 1)$

Family-specific priors can follow the current companion convention, for example:

- Weibull shape on log scale: $\rho_W \sim \mathcal{N}(0, 0.35^2)$
- log-normal scale: $\log(\sigma_{LN}) \sim \mathcal{N}(0, 0.50^2)$
- log-logistic shape parametrization: $\psi_{LL} \sim \mathcal{N}(0, 0.50^2)$

This unified site-prior policy is deliberate: when site terms appear in incidence, latency, or remission, the preferred main prior is $\mathcal{N}(0,1^2)$ and the robustness sensitivity is $t_3(0,1)$.

### 5.10 Likelihood for Stage 8A

The observed-data likelihood contribution is:

$$
L_i = [\pi_i f_u(T_i \mid i)]^{\delta_i} \left[(1-\pi_i) + \pi_i S_u(T_i \mid i)\right]^{1-\delta_i}
$$

where:

- $$
\delta_i = I(\text{status\_num}_i = 1)
$$

### 5.11 Plain-language reading of the Stage 8A likelihood

This likelihood says:

- if a subject transitions, they must belong to the susceptible subgroup and their transition time must come from the susceptible latency density;
- if a subject is censored, they may be either truly non-susceptible or susceptible but not yet transitioned by their censoring time.

That is the core statistical logic of the main Bayesian cure branch.

---

## 6. Stage 8B / Stage 8-CR — Bayesian remission-sensitive competing-risk cure extension

### 6.1 Purpose

Stage 8B exists because remission can materially change the clinical meaning of long-horizon transition risk, especially when remission is common and is not merely administrative censoring. It is therefore a remission-sensitive extension that should be compared side by side with Stage 8A, not substituted for it.

### 6.2 Interpretive meaning of the cure fraction under Stage 8B

Under Stage 8B, the cure fraction must be interpreted as:

> **non-susceptibility to transition**

It must **not** be described as:

- immunity from both transition and remission,
- or guaranteed all-event-free survival.

This distinction is essential. A subject may be non-susceptible to transition and still experience remission.

### 6.3 Observed-state rule

Stage 8B uses all three observed states explicitly:

- $D_i = 1$ → transition
- $D_i = 2$ → remission
- $D_i = 0$ → right censoring

### 6.4 Latent structure

Introduce latent transition susceptibility:

- $U_i \in \{0,1\}$
- $U_i = 1$ means transition-susceptible
- $U_i = 0$ means transition non-susceptible

Let:

- $T_{1i}$ = latent transition time
- $T_{2i}$ = latent remission time
- $C_i$ = right-censoring time

Structural rule:

- if $U_i = 0$, then $T_{1i} = \infty$
- remission is still allowed whether $U_i = 0$ or $U_i = 1$

Observed time:

- $$
Y_i = \min(T_{1i}, T_{2i}, C_i)
$$

### 6.5 Transition incidence component in Stage 8B

To maximize comparability, Stage 8B reuses the same incidence architecture as Stage 8A:

- same incidence design matrix $X_i^{inc}$
- same `anchor_informed` prior branch
- same `neutral_weakly_informative` prior branch
- same merged incidence-site option

### 6.6 Transition latency component in Stage 8B

To maximize comparability, Stage 8B reuses the same transition latency structure and family set as Stage 8A:

- same PNU/SNU branches `L0` and `L1`
- same merged branches `L0S0`, `L1S0`, `L0S1`, `L1S1`
- same site-placement labels
- same transition families: exponential, Weibull, log-normal, log-logistic

### 6.7 Remission component

#### 6.7.1 Principle

The remission submodel should be simpler than the transition-latency grid. The purpose is remission-sensitive realism, not building a second fully saturated survival-family universe.

#### 6.7.2 Recommended primary remission model

Use a piecewise-exponential remission hazard:

$$
h_{2i}(t) = \exp(\eta_i^{rem}) \rho_j, \quad t \in [\kappa_{j-1}, \kappa_j)
$$

Default intervals:

- $[0,1)$
- $[1,2)$
- $[2,5)$
- $[5,10)$
- $[10,\infty)$

#### 6.7.3 Remission linear predictor

PNU / SNU:

- `R0`: $$
\eta_i^{rem} = \xi_0 + \xi_a a_i + \xi_m z_i
$$
- `R1`: $$
\eta_i^{rem} = \xi_0 + \xi_a a_i + \xi_m z_i + \xi_{am} a_i z_i
$$

Merged:

- `R0S0`: $$
\eta_i^{rem} = \xi_0 + \xi_a a_i + \xi_m z_i
$$
- `R1S0`: $$
\eta_i^{rem} = \xi_0 + \xi_a a_i + \xi_m z_i + \xi_{am} a_i z_i
$$
- `R0S1`: $$
\eta_i^{rem} = \xi_0 + \xi_a a_i + \xi_m z_i + \xi_s s_i
$$
- `R1S1`: $$
\eta_i^{rem} = \xi_0 + \xi_a a_i + \xi_m z_i + \xi_{am} a_i z_i + \xi_s s_i
$$


#### 6.7.4 Recommended priors for remission coefficients

- $\xi_0 \sim \mathcal{N}(0, 2.5^2)$
- non-site remission regression coefficients $\xi_j \sim \mathcal{N}(0, 1.5^2)$
- when included in merged remission-site branches: $\xi_s \sim \mathcal{N}(0, 1^2)$ in the main site-prior family
- mandatory site-prior sensitivity for merged remission-site branches: $\xi_s \sim t_{3}(0, 1)$
- $\log(\rho_j) \sim \mathcal{N}(0, 1.5^2)$ for piecewise baseline remission hazards

### 6.8 Practical identifiability assumption

The minimal practical assumption for Stage 8B is:

> conditional on observed covariates, remission is modeled independently of the latent transition susceptibility process except through observed competition in the likelihood.

This is a **pragmatic identification choice**, not a literal biological claim that remission and transition are mechanistically independent.

### 6.9 Likelihood under the competing-risk cure extension

Define:

- $S_{1u,i}(t)$ = transition survival for susceptible subjects
- $f_{1u,i}(t)$ = transition density for susceptible subjects
- $h_{2i}(t)$ = remission cause-specific hazard
- $$
G_{2i}(t) = \exp\left\{-\int_0^t h_{2i}(u)\,du\right\}
$$
= remission-free survival
- $$
M_{1i}(t) = (1-\pi_i) + \pi_i S_{1u,i}(t)
$$
= transition-free mixture survival

#### 6.9.1 Transition event contribution

If $D_i = 1$:

$$
L_i^{(1)} = \pi_i f_{1u,i}(Y_i) G_{2i}(Y_i)
$$

**Meaning:** the subject must be transition-susceptible, must transition at $Y_i$, and remission must not have happened first.

#### 6.9.2 Remission event contribution

If $D_i = 2$:

$$
L_i^{(2)} = h_{2i}(Y_i) G_{2i}(Y_i) M_{1i}(Y_i)
$$

**Meaning:** remission occurs at $Y_i$, remission had not happened earlier, and transition had not already occurred by $Y_i$. Transition-free survival here is a mixture because some subjects are non-susceptible to transition and others are susceptible but still event-free.

#### 6.9.3 Right censoring contribution

If $D_i = 0$:

$$
L_i^{(0)} = G_{2i}(Y_i) M_{1i}(Y_i)
$$

**Meaning:** by censoring time $Y_i$, neither remission nor transition has yet occurred.

### 6.10 Implied population functions

The remission-sensitive branch naturally implies three core functions.

#### 6.10.1 Transition CIF

$$
F_{1i}(t) = \int_0^t G_{2i}(u) \pi_i f_{1u,i}(u)\,du
$$

This is the cumulative probability that transition has happened by time $t$, explicitly accounting for remission as a competing outcome.

#### 6.10.2 Remission CIF

$$
F_{2i}(t) = \int_0^t h_{2i}(u) G_{2i}(u) M_{1i}(u)\,du
$$

This is the cumulative probability that remission has happened by time $t$.

#### 6.10.3 All-event-free survival

$$
S_i^{all}(t) = G_{2i}(t) M_{1i}(t)
$$

This is the probability that by time $t$, neither transition nor remission has occurred.

### 6.11 Plain-language interpretation of the three functions

These three quantities partition the subject’s future into:

- already transitioned,
- already remitted,
- neither yet happened.

This is why the Stage 8B branch is more appropriate than Stage 8A when the scientific question is the **real-world probability of transition under explicit competition from remission**.

### 6.12 Prediction target under Stage 8B

The primary risk quantity is:

- $$
\text{transition\_CIF}(t)
$$

not the Stage 8A censoring-based transition-only risk.

This point is non-negotiable. Once remission is modeled as a true competing event, prediction, thresholding, and decision-curve calculations must be based on the transition CIF.

### 6.13 Threshold-based clinical-usefulness outputs

For Stage 8B, threshold-based outputs should be based on $$
\text{transition\_CIF}(t)
$$
including:

- positive classification rate
- false-positive count
- false positives per 100
- PPV where meaningful
- TPR where estimable
- false-positive burden
- net benefit

This is necessary because using a censoring-based event probability after explicitly introducing remission as a competing event would mix scales and produce an incoherent comparison.

---


## 7. Mandatory Stage 8 export contract, diagnostics, and interpretation guardrails

### 7.1 Export tables required from every admissible Stage 8 fit

At minimum, Stage 8 must export:

- one long-format posterior prediction table
- one horizon-level performance / classification table
- one prior-sensitivity comparison table
- one anchor-informed versus neutral / no-external-information delta table
- one prior-to-posterior incidence-shape update table
- one hazard-shape plausibility table
- one Stage 8A-versus-Stage 8B delta table
- one uncured-only / non-cure-fraction supporting decomposition table
- one diagnostics / admissibility table
- figure-ready source tables for the `horizon_support_panel`, the `8A_vs_8B_delta_panel`, the `anchor_vs_neutral_delta_panel`, the `incidence_anchor_update_panel`, and the `uncured_only_decomposition_panel`

Every horizon-level row must preserve at minimum:

- `dataset_key`
- `branch`
- `risk_scale`
- `prior_branch`
- `site_prior_family` where a site term is present
- `horizon`
- `threshold` where applicable
- `support_tier`
- `horizon_evidence_class`
- `claim_restriction_flag`
- `prior_tail_sensitive`
- `admissibility_flag`

The horizon-level performance / classification table should include where estimable:

- posterior risk, survival, or CIF
- posterior cure fraction summaries
- posterior `susceptible_fraction = 1 - cure_fraction`
- posterior uncured / susceptible-only survival and risk on the annual 1–10-year grid
- posterior `DeltaRisk_NC_C(t)` where defined
- posterior threshold-based false-positive burden
- false positives per 100
- PPV
- TPR
- net benefit
- discrimination summaries
- calibration summaries
- Brier / IBS summaries

### 7.2 Horizon-support governance inside Stage 8

Stage 8 must restate the master specification’s horizon-support governance directly in its own exports.

Only rows labeled `directly_observed_data_supported` and `primary_claim_allowed` are eligible for the strongest claims.

Rows labeled `partly_model_dependent` may support only qualified secondary or sensitivity claims.

Rows labeled `mostly_extrapolated`, `projection_only`, or `projection_plus_prior_sensitive` may not be used as observed-fact claims.

### 7.3 Prior-tail action rule

The project uses two incidence-prior branches for Stage 8:

- `anchor_informed`
- `neutral_no_external_info`

Legacy implementation alias allowed for the second branch:

- `neutral_weakly_informative`

That design is not only descriptive. It creates an explicit action rule.

For each dataset, branch, family, and horizon, Stage 8 must compare the two prior branches and export prior-sensitivity contrasts for at least:

- risk or CIF
- cure fraction
- false-positive burden
- false positives per 100
- net benefit

Set:

- `prior_tail_sensitive = TRUE`

whenever prespecified materiality thresholds are exceeded **or** a clinically relevant threshold-crossing / sign-of-net-benefit conclusion changes between the two prior branches.

The exact numerical materiality thresholds should be fixed in the implementation metadata registry, but the action rule itself is fixed here.

If:

- `prior_tail_sensitive = TRUE`
- and `horizon_evidence_class = mostly_extrapolated`

then the row must be relabeled:

- `claim_restriction_flag = projection_plus_prior_sensitive`

Such rows may be discussed only as projection-sensitive sensitivity results and may not support primary or headline claims.

### 7.3A How to distinguish genuine stabilization from prior-driven illusion

A visually stable posterior cure fraction is not sufficient evidence of genuine cure-like structure.

Stage 8 should be interpreted as genuine stabilization only when:

1. posterior quantities move meaningfully away from the prior center;
2. supported-horizon risk, delta-risk, and threshold-based outputs stabilize together;
3. directional conclusions remain coherent across latency families;
4. the signal persists after alternative explanations such as site/context and remission handling are examined.

### 7.4 Hazard-shape plausibility export

The project should not evaluate tail behavior only by which family “wins” on a fit index.

Every retained Stage 8 fit must therefore export, on the common annual horizon grid:

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

These are not merely descriptive extras. They are interpretation guardrails so that a numerically admissible but biologically implausible tail is visible in the exported results.

### 7.5 Mandatory anchor-versus-neutral delta exports

Because integration of external meta-analysis is a central scientific feature of Stage 8, every retained primary Stage 8A fit and its Stage 8B mirror must be compared under both prior branches.

At minimum, the following horizon-specific delta fields must be exported wherever estimable:

- `delta_risk_anchor_minus_neutral`
- `delta_cure_fraction_anchor_minus_neutral`
- `delta_false_positive_burden_anchor_minus_neutral`
- `delta_FP100_anchor_minus_neutral`
- `delta_NB_anchor_minus_neutral`
- `delta_PPV_anchor_minus_neutral`
- `delta_TPR_anchor_minus_neutral`

The point of this comparison is not to declare the anchored branch the winner. The point is to show what the external meta-analytic anchor contributes and whether supported-horizon conclusions remain coherent after the anchor is removed.

### 7.5A Mandatory prior-to-posterior incidence-shape update exports

Stage 8 must also export figure-ready source data showing how the initial anchored incidence shape changes after estimation.

At minimum, the incidence-shape update table must contain:

- `dataset_key`
- `branch`
- `risk_scale`
- `prior_branch`
- `retained_fit_id`
- `age_sex_anchor_cell`
- `external_incidence_rate_per10k`
- `external_one_year_risk`
- `prior_center_logit`
- `posterior_mean_logit`
- `posterior_lower_logit`
- `posterior_upper_logit`
- `posterior_mean_one_year_risk`
- `posterior_lower_one_year_risk`
- `posterior_upper_one_year_risk`
- `posterior_minus_prior_logit`
- `posterior_minus_prior_risk`

This panel is mandatory because otherwise the project’s defining external anchor can disappear into the posterior without readers seeing how much updating actually occurred.

### 7.5B Auxiliary uncured-only / non-cure-fraction exports

Every retained cure-model fit should also export a cure-model-only supporting decomposition table containing:

- `cure_fraction`
- `susceptible_fraction`
- `uncured_survival(t)` on the annual 1–10-year grid
- `uncured_risk(t)` on the annual 1–10-year grid
- optional `MSTu`
- `uncured_mean_support_flag`

These are supporting decomposition estimands, not part of the main shared cross-model ruler.

If a full `MSTu` would require aggressive tail extrapolation or is not supportable after follow-up screening, Stage 8 should avoid forcing a single-number uncured-only mean and instead keep the uncured-only outputs on the annual horizon grid only.

### 7.6 Mandatory Stage 8A versus Stage 8B delta exports

The project’s scientific question for the competing-risk extension is not “which branch is better?” but “how much does remission-aware modeling rewrite the apparent cure-like signal and downstream classification burden?”

Therefore, Stage 8 must export horizon-specific deltas for at least:

- `delta_risk_8B_minus_8A`
- `delta_cure_fraction_8B_minus_8A`
- `delta_false_positive_burden_8B_minus_8A`
- `delta_FP100_8B_minus_8A`
- `delta_NB_8B_minus_8A`
- `delta_PPV_8B_minus_8A`
- `delta_TPR_8B_minus_8A`

### 7.7 Shared Bayesian diagnostics

Neither Stage 8A nor Stage 8B should be considered admissible if any of the following are unacceptable:

- serious divergence problems
- poor chain mixing
- unacceptable `Rhat`
- inadequate effective sample size
- gross posterior-predictive contradiction in the observed support
- clear prior or posterior degeneracy

### 7.8 Additional Stage 8B coherence check

For each exported horizon $k$, require:

- $0 \le F_1(k) \le 1$
- $0 \le F_2(k) \le 1$
- $0 \le S_{\text{all}}(k) \le 1$
- $\left|F_1(k) + F_2(k) + S_{\text{all}}(k) - 1\right|$ stays within a small numerical tolerance

### 7.9 Retention rule

Stage 8 should not force one Bayesian winner. The rule is:

1. exclude non-admissible fits;
2. retain all admissible fits;
3. export all admissible fits into the later unified comparison engine;
4. treat internal Bayesian fit statistics as supportive evidence only.

### 7.10 Allowed claims and forbidden claims

Allowed claims include statements such as:

- at supported horizons, allowing cure-like heterogeneity changed estimated transition risk and altered false-positive burden or clinical usefulness by a stated amount;
- anchor-informed and neutral / no-external-information branches gave similar or different supported-horizon conclusions by a stated amount;
- the posterior-updated incidence pattern remained close to, or moved away from, the external anchor in a stated way;
- remission-aware competing-risk analysis changed the apparent cure-like signal or downstream classification burden by a stated amount.

Forbidden claims include:

- “the Bayesian cure model proved cure”;
- “posterior stabilization alone confirms a discrete cured subgroup”;
- “the external meta-analysis prior proved the true cure fraction”;
- treating projection-dominant horizons as observed fact;
- treating a subject-level posterior cure probability as a stand-alone patient decision number in the current project;
- describing a merged site coefficient as a causal treatment effect;
- describing the Stage 8A-versus-Stage 8B delta table as a fair within-scale model-ranking competition;
- presenting the anchor-informed branch as automatically correct simply because it uses a larger external data source.

### 7.11 Required figure-ready source tables

Stage 8 must export figure-ready source tables for:

- `horizon_support_panel`
  - horizon
  - `support_tier`
  - `horizon_evidence_class`
  - `claim_restriction_flag`
  - `prior_tail_sensitive`
  - risk-set support where available

- `8A_vs_8B_delta_panel`
  - horizon
  - dataset
  - `delta_risk_8B_minus_8A`
  - `delta_cure_fraction_8B_minus_8A`
  - `delta_false_positive_burden_8B_minus_8A`
  - `delta_FP100_8B_minus_8A`
  - `delta_NB_8B_minus_8A`
  - `delta_PPV_8B_minus_8A`
  - `delta_TPR_8B_minus_8A`

- `anchor_vs_neutral_delta_panel`
  - horizon
  - dataset
  - `branch`
  - `risk_scale`
  - `delta_risk_anchor_minus_neutral`
  - `delta_cure_fraction_anchor_minus_neutral`
  - `delta_false_positive_burden_anchor_minus_neutral`
  - `delta_FP100_anchor_minus_neutral`
  - `delta_NB_anchor_minus_neutral`
  - `delta_PPV_anchor_minus_neutral`
  - `delta_TPR_anchor_minus_neutral`
  - `support_tier`
  - `claim_restriction_flag`

- `incidence_anchor_update_panel`
  - `dataset_key`
  - `branch`
  - `prior_branch`
  - `age_sex_anchor_cell`
  - `external_incidence_rate_per10k`
  - `external_one_year_risk`
  - `prior_center_logit`
  - `posterior_mean_logit`
  - `posterior_lower_logit`
  - `posterior_upper_logit`
  - `posterior_mean_one_year_risk`
  - `posterior_lower_one_year_risk`
  - `posterior_upper_one_year_risk`
  - `posterior_minus_prior_logit`
  - `posterior_minus_prior_risk`

- `uncured_only_decomposition_panel`
  - `dataset_key`
  - `branch`
  - `risk_scale`
  - `cure_fraction`
  - `susceptible_fraction`
  - annual `uncured_survival(t)` / `uncured_risk(t)` fields in long format
  - optional `MSTu`
  - `uncured_mean_support_flag`

These source tables are mandatory because otherwise a 10-year projection can too easily be read with the same visual weight as a directly observed 1-year estimate, and the defining external incidence anchor can disappear from view.


## 8. What is standard methodology and what is project-specific

### 8.1 Standard or literature-aligned components

The following are standard or literature-aligned:

- mixture cure decomposition into incidence and latency;
- Bayesian weakly informative priors for regression coefficients after sensible scaling;
- use of 2-SD scaling for continuous covariates to improve coefficient interpretability and prior specification;
- use of parametric latency families such as Weibull, log-normal, and log-logistic;
- the AFT-versus-Cox discussion about latency support and insufficient follow-up;
- use of CIF-based decision curves and clinical usefulness under competing risks;
- piecewise-exponential hazards as a practical baseline model for a competing event process;
- separate reporting of cure-fraction-type summaries and uncured-only survival summaries as supporting decomposition outputs when a cure model is used.

### 8.2 Project-specific but methodologically defensible components

The following are project-specific design choices:

- the explicit split into **Stage 8A** and **Stage 8B** rather than replacing the main branch;
- the dual-estimand architecture with `transition_only_main` and `transition_cif_competing`;
- the decision to keep the incidence prior anchored to the van der Werf age–sex pattern while adding a mandatory neutral / no-external-information branch;
- the decision to require an explicit prior-to-posterior incidence-shape update panel rather than hiding the anchor inside the posterior;
- the use of the same shared horizon and threshold grid across KM, no-cure, frequentist cure, Bayesian cure, and Bayesian competing-risk cure;
- the choice to interpret the Stage 8B cure fraction specifically as **non-susceptibility to transition**;
- the decision to keep the remission submodel simpler than the transition family grid;
- the rule that uncured-only / non-cure-fraction summaries remain supporting decomposition estimands rather than the main common-scale ruler;
- the decision to keep frailty-based and latent-class heterogeneity explanations outside baseline Stage 8 interpretation because the current data structure is not a stable basis for treating them as baseline comparators.

### 8.3 How this should be described honestly

This model should therefore be described as:

> a **literature-informed, project-specific Bayesian extension** that combines standard mixture-cure structure, weakly informative prior design, explicit external-incidence prior centering, paired no-external-information comparison, and competing-risk prediction logic while preserving compatibility with the project’s shared comparison ruler.

It should **not** be described as:

- a direct one-to-one implementation of a single canonical published competing-risk cure model,
- an off-the-shelf textbook model copied unchanged,
- a latent-class survival model,
- or a frailty-based cure model.

That distinction improves, rather than weakens, the credibility of the specification because it makes the modeling choices transparent.


## 9. Literature-grounded justification by design choice

### 9.1 Cure-model decomposition and Bayesian cure modeling

The overall cure-model decomposition, and the idea of Bayesian cure modeling, are grounded in the Peng and Yu cure-model texts:

- **The parametric cure model**
- **The semiparametric and nonparametric cure models**
- **Bayesian cure model**

These texts support the basic incidence–latency factorization and the general Bayesian cure-model framework.

### 9.2 Incidence prior centering and age–sex structure

The use of age-band by sex incidence structure as external centering comes from the van der Werf incidence synthesis. In the present project it is used **only** to stabilize and inform the incidence part, not to dictate the cure fraction as fixed truth.

The conversion path is deliberately explicit:

1. external annual incidence rates are taken on a per-10 000 person-year scale;
2. they are converted to one-year risks;
3. those risks are moved to the logit scale;
4. the resulting values are represented on the project’s incidence design matrix as a frozen contrast vector.

That last step matters. The stored vector $\mu^{vdW}$ is not just a list of epidemiologic rates copied into the model; it is a project-specific contrast representation of the external template under the current incidence design matrix.

### 9.3 Weakly informative priors and predictor scaling

The use of weakly informative priors and 2-SD scaling is aligned with Gelman’s work on:

- weakly informative priors for logistic and related regression models;
- scaling regression inputs by dividing continuous predictors by two standard deviations.

These references justify both the prior scale logic and the scaling discipline used in the companion.

### 9.4 Why AFT structure matters when follow-up is uneven

The AFT-versus-Cox cure discussion in Parsa and Van Keilegom supports the concern that Cox-latency cure models impose a common support structure that may be unrealistic when follow-up sufficiency varies across covariate regions, while AFT-style structures can be more realistic in such settings.

Disagreement across latency families should be interpreted first as evidence of latency-side heterogeneity, tail identifiability limits, or both.

It should not be treated as direct evidence of a stable discrete cured subgroup unless the separation also appears in supported horizons and remains coherent under prior sensitivity, site adjustment, and remission-sensitive comparison.

### 9.5 Why Stage 8B must use CIF for prediction and decision curves

The need to base competing-risk clinical usefulness on the CIF rather than a censoring-based event probability is aligned with the decision-curve work of Vickers and colleagues and with the project’s own dual-estimand framework.

### 9.6 Why Stage 8B is supplementary rather than replacing Stage 8A

The integrated dual-estimand framework fixes this rule: the remission-sensitive branch should be added in a way that preserves the main transition-only backbone, rather than silently overwriting it.

### 9.7 Why uncured-only / non-cure-fraction summaries are supporting rather than primary

Recent cure-model comparison work argues that, in the presence of a cure fraction, it can be more informative to compare cure fractions together with the survival of the uncured rather than relying only on overall survival, including through mean survival time of the uncured (`MSTu`).

That supports the current project’s cure-model-only supporting decomposition block.

However, because such estimands condition on uncured status, they should remain supporting rather than replacing the main common-scale comparison ruler.

### 9.8 Why frailty-based and latent-class explanations remain background only

The current Stage 8 companion is not a latent-class survival model and not a frailty-based cure model.

In the present project, such higher-order heterogeneity formulations are acknowledged only as future directions because the current univariate survival setting and weak late-tail support are not a stable basis for treating them as baseline comparators.

Accordingly:

- do not describe Stage 8 as identifying latent classes;
- do not describe posterior regularization as evidence for a frailty distribution;
- and do not use Stage 8 outputs to back-fit a frailty or latent-class narrative that is stronger than the observed data can support.

## 10. Minimum reporting language this specification supports

### 10.1 Methods-facing language for Stage 8A

> We fitted Bayesian mixture-cure models on the project’s main transition-only scale, in which transition was the event of interest and remission was treated as censoring. The incidence component used a prespecified age-band by sex prior-centering strategy informed by external incidence patterns only as a regularizing prior on the incidence component, not as a prior assertion about the true cure fraction or as fixed clinical truth. The latency component remained aligned with the common project covariate backbone. Weakly informative priors were used throughout, and both anchor-informed and neutral / no-external-information prior branches were fitted as prespecified prior sensitivities.

### 10.2 Methods-facing language for Stage 8B

> We additionally fitted a remission-sensitive Bayesian competing-risk cure extension in which transition remained the event of interest, remission was modeled explicitly as a competing event, and the cure fraction was interpreted as non-susceptibility to transition rather than immunity from all events. The transition incidence component retained the same paired anchor-informed versus no-external-information prior structure used in the main Bayesian branch, while remission was modeled using a simpler piecewise-exponential submodel. Threshold-based prediction and decision-curve quantities on this branch were based on the transition cumulative incidence function.

### 10.3 Discussion-facing caution language

> The remission-sensitive Bayesian branch should be interpreted as a literature-informed project extension rather than a direct off-the-shelf canonical competing-risk cure model. Its main value is not to replace the transition-only analysis, but to quantify how much the apparent cure-like signal and downstream false-positive burden change once remission is treated as an explicit competing event.

> Late-horizon Bayesian outputs must be read together with the exported `support_tier`, `horizon_evidence_class`, `claim_restriction_flag`, and `prior_tail_sensitive` fields. Projection-dominant or prior-sensitive tails are not headline evidence.

> Posterior stabilization of the cure fraction, by itself, should not be interpreted as confirmation of a discrete cured subgroup. In this project, such stabilization must be read jointly with supported-horizon risk shifts, family sensitivity, Stage 6 carry-forward screening, anchor-versus-neutral comparison, and the remission-sensitive comparison.

> The current Bayesian companion should also not be over-read as latent-class identification or frailty identification. Those interpretations remain outside the present baseline workflow.

## 11. One-line lock

**Stage 8 should now be understood as a dual-branch Bayesian companion: Stage 8A remains the main transition-only cure model on the shared project backbone, while Stage 8B adds a remission-sensitive competing-risk extension whose core likelihood, prediction target, and clinical-usefulness outputs are explicitly defined on the transition CIF scale and whose interpretation is anchored to non-susceptibility to transition rather than immunity from all events; every exported horizon-specific row must also carry explicit horizon-support, claim-restriction, prior-tail, anchor-versus-neutral, and remission-aware delta metadata; the external incidence anchor must be shown both before and after posterior updating; uncured-only / non-cure-fraction summaries must remain in a separate supporting decomposition block rather than being overclaimed as the main ruler; and neither latent-class nor frailty-based explanations are established by this baseline Stage 8 workflow.**

## 12. Reference basis used in this enriched companion

Core project architecture and branch logic:

- `Integrated_Modeling_Master_Specification_English_REVISED_v5.md`
- `Bayesian_Modeling_Specification_Stage8_REVISED_v5.md`

Methodological references used to justify the modeling choices:

- `R14A. Peng Y, Yu B. The parametric cure model. Cure Models - Methods, Applications, and Implementation.`
- `R14B. Peng Y, Yu B. The semiparametric and nonparametric cure models. Cure Models - Methods, Applications, and Implementation.`
- `R14C. Peng Y, Yu B. Bayesian cure model. Cure Models - Methods, Applications, and Implementation.`
- `2014_Systematic review and collaborative reca_Van Der Werf et al..pdf`
- `2008_A weakly informative default prior distr_Gelman et al..pdf`
- `R1. Gelman A. Scaling regression inputs by dividing by two standard deviations. Stat Med. 2008.`
- `R10. Parsa M, Van Keilegom I. Accelerated failure time vs Cox proportional hazards mixture cure models. Statistical Papers. 2023.`
- `R12. Vickers AJ, Cronin AM, Elkin EB, Gonen M. Extensions to decision curve analysis ... 2008.`
- `R13. McLernon DJ et al. Assessing performance and clinical usefulness in prediction models with survival outcomes.`
- `2025_Flexible Parametric Accelerated Failure_Akynkozhayev et al..pdf`
- `2022_An Accelerated Failure Time Cure Model w_Aida et al..pdf`
- `2020_A tutorial on frailty models_Balan and Putter.pdf`
- `2024_A latent class Cox model for heterogeneo_Pei et al..pdf`
- `2024_A two-sample comparison of mean survival_Dobler and Musta.pdf`

The last five references are not required to define the core Stage 8 architecture, but they are useful as extended background when discussing flexible parametric AFT structures, frailty and latent heterogeneity background, and uncured-only supporting summaries.

