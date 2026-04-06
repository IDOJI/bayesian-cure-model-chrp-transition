
# CHR Cohorts Modeling Specification

## Final version based on the 6-block workflow

## 1. Scope and governing principle

### 1.1 Purpose

본 분석의 목적은 CHR 코호트에서 장기 transition probability를 추정하고 비교하는 것이다. 분석은 다음 **6개 블록**으로 구성한다.

1. **KM, MLE, Bayesian probability comparison**
    
2. **Timing-difference separation**
    
3. **Short vs long follow-up direct comparison**
    
4. **Follow-up adequacy / tail stability assessment**
    
5. **Cure-appropriateness screening**

6. **Cure decomposition stability / latency fragility diagnostics**
    

### 1.2 Governing principle

최종 과학적 해석은 위 여섯 블록을 **순차적으로 읽은 뒤 종합**해서 제시한다.  
어느 한 블록도 단독으로 cure-like heterogeneity를 확정하지 않는다.  
특히 **timing difference**, **follow-up immaturity / tail instability**, **cure-like heterogeneity**는 분리 해석한다.

### 1.3 Relationship to prior source files

기존 source file과 integrated master specification은 **참고 문헌과 해석 경계 설정용 reference**로만 사용한다.  
현재 분석 운영은 이 문서의 6-block workflow를 따른다.

---

## 2. Data backbone

### 2.1 Datasets

분석 대상 데이터셋은 다음 세 개다.

- `PNU`
    
- `SNU`
    
- `merged`
    

### 2.2 Time scale

기본 time scale은 cohort entry 이후 follow-up time이다.

- `time = days_followup`
    
- 필요 시 `time_year = days_followup / 365.25`
    

### 2.3 Main event definition

기본 분석의 main risk scale은 **transition-only**로 둔다.

- event of interest: `status_num == 1`
    
- censoring in the main analysis: `status_num %in% c(0, 2)`
    

즉 기본 spec에서는 transition 이외의 outcome은 모두 censoring으로 처리한다.

### 2.4 Core covariates

공통 기본 공변량은 다음으로 둔다.

- `age_s`
    
- `sex_num`
    
- `site` (merged 및 timing-difference contexts)
    

### 2.5 Interpretation rule for site

site effect는 causal treatment effect로 해석하지 않는다.  
site는 referral structure, care pathway, selection, follow-up maturity를 반영하는 **context proxy**로만 읽는다.

---

## 3. Cross-block hierarchy

### 3.1 Main results block

다음 두 블록은 **main results**에 해당한다.

- **Block 1: KM, MLE, Bayesian probability comparison**
    
- **Block 3: Short vs long follow-up direct comparison**
    

이 두 블록이 본문의 핵심 메시지를 만든다.

### 3.2 Supporting interpretation block

다음 네 블록은 **additional analyses for interpretation**에 해당한다.

- **Block 2: Timing-difference separation**
    
- **Block 4: Follow-up adequacy / tail stability assessment**
    
- **Block 5: Cure-appropriateness screening**

- **Block 6: Cure decomposition stability / latency fragility diagnostics**
    

이 네 블록은 main results를 해석하고 방어하는 구조적 보조분석이다.

### 3.3 Final synthesis rule

최종 결론은 Block 1–6을 모두 본 뒤에만 제시한다.  
어느 한 블록도 단독으로 final scientific conclusion을 담당하지 않는다.

---

## 4. Block 1 — KM, MLE, and Bayesian probability comparison

### 4.1 Objective

이 블록의 목적은 **모델별 장기 transition probability 차이를 직접 보여주는 것**이다.

즉, 어떤 모델이 맞는가를 먼저 선언하기보다,

- KM
    
- non-cure MLE
    
- frequentist cure MLE
    
- Bayesian cure
    

가 장기 확률을 어떻게 다르게 제시하는지를 비교한다.

### 4.2 KM benchmark

KM은 descriptive benchmark로 포함한다.

반드시 제시:

- dataset별 overall KM
    
- numbers at risk
    
- 95% CI where feasible
    
- common horizons에서의 `S_KM(t)`와 `1 - S_KM(t)`
    

Betensky도 KM, 95% CI, numbers at risk를 기본 display로 둔 뒤 follow-up stability를 추가 논의한다.

### 4.3 Non-cure MLE comparator

비교용 non-cure parametric model은 **prespecified lognormal no-cure model**로 고정한다.

이 workflow에서는 latency-side parametric assumption을 block 간에 일관되게 유지하기 위해, non-cure comparator에도 lognormal specification만 사용한다.
    

### 4.4 Frequentist cure model

기본 frequentist cure model은 **prespecified lognormal AFT mixture cure model**으로 둔다.

- incidence submodel: logistic
    
- latency submodel: AFT
    
- latency family: lognormal
    

Peng & Yu는 mixture cure model의 기본 incidence–latency decomposition을 설명하고, incidence link로 logit을 가장 흔한 선택으로 두며, parametric cure chapter에서 lognormal AFT mixture cure model을 직접 예시로 적합한다.

### 4.5 Why AFT only, not Cox

latency 구조에서는 Cox를 main model로 사용하지 않고 **AFT만 채택**한다.

Parsa & Van Keilegom은 Cox mixture cure model이 공통 cure-threshold/support를 강하게 가정할 수 있고, follow-up이 일부 covariate region에서 불충분할 때 AFT mixture cure model이 더 현실적일 수 있다고 설명한다. 두 모델은 충분한 follow-up 영역에서는 비슷할 수 있지만, 부족한 영역 밖에서는 크게 달라질 수 있다. 본 연구의 핵심 질문이 short vs long follow-up이므로, AFT-only latency choice를 main specification으로 채택한다.

### 4.6 Bayesian model

Bayesian cure model은 병렬 probability comparison block에 포함한다.

다만 이 기본 spec에서는 Bayesian prior, anchor, posterior workflow를 다시 정의하지 않는다.  
Bayesian modeling의 세부 구현은 **별도 Bayesian specification**만 따른다. Bayesian companion은 external anchor와 posterior-updated anchor를 별도로 유지할 것을 요구한다.

### 4.7 Required outputs

이 블록의 산출물은 다음이다.

- dataset별 KM curve
    
- dataset별 non-cure MLE fitted curves
    
- dataset별 frequentist cure fitted curves
    
- dataset별 Bayesian cure fitted curves
    
- common horizons에서의 survival / risk probability table
    
- model-by-dataset probability difference summary
    

### 4.8 Allowed interpretation

이 블록에서 허용되는 해석은 다음이다.

- 모델에 따라 장기 probability 추정이 얼마나 달라지는가
    
- no-cure vs cure-aware 계열이 tail probability를 어떻게 다르게 그리는가
    
- Bayesian model이 uncertainty를 어떻게 반영하는가
    

### 4.9 Forbidden interpretation

이 블록만으로 cure appropriateness를 확정하지 않는다.

---

## 5. Block 2 — Timing-difference separation

### 5.1 Objective

이 블록의 목적은 **PNU–SNU 차이가 early/intermediate timing difference인지**, 아니면 later tail structure와 얽혀 있는지를 분리하는 것이다. 기존 integrated spec의 Stage 4도 공식 목적을 이렇게 둔다.

### 5.2 Standard backbone

표준 출력은 다음으로 고정한다.

- `0–1 year` interval-specific site contrast
    
- `1–2 years` interval-specific site contrast
    
- `2–5 years` interval-specific site contrast
    
- `>5 years` tail diagnostic only
    
- `1-year restricted analysis`
    
- `2-year restricted analysis`
    

`>5 years`는 **tail diagnostic only**이며, primary timing interpretation을 주도하면 안 된다.

### 5.3 Optional extension

필요 시 `2–3-year-only slice`를 추가할 수 있다.  
다만 이는 **optional sensitivity / descriptive extension**이며, 표준 backbone을 대체하지 않는다.

### 5.4 Role in the full argument

이 블록은 **independent cure evidence가 아니다.**  
역할은 short follow-up 해석이 단순한 early timing difference의 결과인지 분리하는 것이다.

### 5.5 Required outputs

- interval-specific site contrasts
    
- 1-year restricted analysis summary
    
- 2-year restricted analysis summary
    
- `>5 years` tail-diagnostic summary
    
- concise timing-pattern interpretation
    

### 5.6 Allowed interpretation

- 초기/중간 구간에서 이미 cohort timing 차이가 존재하는가
    
- later long-tail difference 이전에 short-window 차이가 있었는가
    

### 5.7 Forbidden interpretation

- Stage 4 alone = cure evidence
    
- Stage 4 alone = heterogeneity proof
    

---

## 6. Block 3 — Short vs long follow-up direct comparison

### 6.1 Objective

이 블록은 **2–3년 추적만 썼을 때와 장기 추적을 썼을 때 probability와 cure-like signal이 실제로 얼마나 달라지는지 직접 보여주는 핵심 분석**이다.

### 6.2 Constructed datasets

기본 비교는 다음을 사용한다.

- `SNU_full`
    
- `SNU_truncated_2y`
    
- `SNU_truncated_3y` (optional)
    
- `merged_full`
    
- `merged_truncated_2y`
    
- `merged_truncated_3y` (optional)
    

### 6.3 Refit rule

각 truncated dataset과 full dataset에는 **동일한 모델 세트**를 다시 적합한다.

- KM
    
- non-cure MLE lognormal comparator
    
- frequentist lognormal AFT mixture cure
    
- Bayesian cure
    

### 6.4 Required outputs

- full vs truncated KM
    
- full vs truncated model-based probability curves
    
- common-horizon probability tables
    
- divergence summary between full and truncated analyses
    
- descriptive summary of whether cure-like signal weakens or disappears under truncation
    

### 6.5 Role in the full argument

이 블록은 short-vs-long follow-up narrative의 **직접 근거**다.

### 6.6 Allowed interpretation

- 장기 추적이 있으면 probability 추정이 얼마나 달라지는가
    
- 2–3년 추적만으로는 long-tail information이 얼마나 잃어버려지는가
    

### 6.7 Forbidden interpretation

이 블록만으로 timing difference와 cure-like heterogeneity를 분리하지 않는다.  
그 분리는 Block 2, 4, 5가 보완한다.

---

## 7. Block 4 — Follow-up adequacy / tail stability assessment

### 7.1 Objective

이 블록의 목적은 **late tail이 얼마나 해석 가능한지**를 보여주는 것이다.  
follow-up maturity는 tail maturity indicator이지 direct evidence of cure가 아니다.

### 7.2 Parallel reporting rule

이 블록에서는 어떤 하나를 primary로 고정하지 않는다.  
다음 두 축을 **병렬로 같이 제시**한다.

### 7.3 Reverse-KM / censoring-based summaries

포함:

- reverse KM
    
- time-to-censoring summaries
    
- numbers at risk
    
- cumulative censoring summaries
    

### 7.4 Betensky-style stability displays

포함:

- KM upper/lower limits
    
- difference curve
    
- normalized area under the difference curve
    
- partial difference curves where useful
    

Betensky는 흔히 보고되는 follow-up measures만으로는 KM stability를 직접 전달하기 어렵고, upper/lower limits와 difference curves가 stability를 더 직접적으로 보여준다고 제안한다.

### 7.5 Required outputs

- reverse-KM summary table
    
- numbers-at-risk table
    
- tail instability flags
    
- Betensky-style stability figures
    
- cohort별 late-tail narrative summary
    

### 7.6 Allowed interpretation

- late tail support가 충분한가
    
- plateau처럼 보이는 구간이 얼마나 안정적인가
    
- 추가 follow-up이 있으면 KM tail이 얼마나 달라질 수 있는가
    

### 7.7 Forbidden interpretation

이 블록만으로 cure fraction의 존재를 주장하지 않는다.

---

## 8. Block 5 — Cure-appropriateness screening

### 8.1 Objective

이 블록의 목적은 **cure model이 현재 데이터 구조에서 합리적인가를 screening하는 것**이다.  
이 블록은 main modeling의 앞에서 절대적 gatekeeper처럼 쓰기보다, **Block 1–4 결과를 해석하는 additional-analysis block**으로 둔다.

### 8.2 Screening hierarchy

screening hierarchy는 다음으로 고정한다.

1. **Primary screening method:** RECeUS / RECeUS-AIC
    
2. **Secondary contradiction-oriented follow-up test:** Xie sufficient-follow-up test
    
3. **Older contextual reference:** Maller–Zhou family of methods
    

### 8.3 Primary screening method: RECeUS / RECeUS-AIC

RECeUS는 기존 방법들이 cure presence 또는 sufficient follow-up 중 한쪽만 다루는 한계를 정리한 뒤, `π̂_n`과 `r̂_n`을 함께 사용해 cure-model appropriateness를 평가하도록 제안된 pragmatic method다.

default threshold 조합은 보통 다음처럼 쓴다.

- `π̂_n > 0.025`
    
- `r̂_n < 0.05`
    

하지만 이는 **screening thresholds**로 읽어야 하며, universal theorem처럼 고정하지 않는다.

### 8.4 Secondary contradiction-oriented test: Xie

Xie sufficient-follow-up test는 **secondary follow-up contradiction check**로 포함한다.

역할:

- follow-up sufficiency를 좀 더 직접적으로 확인하는 formal check
    
- stand-alone positive gate로는 쓰지 않음
    
- RECeUS supportive result에 대해 contradiction 여부를 표시하는 보조 검정으로 사용
    

### 8.5 Older contextual reference: Maller–Zhou

Maller–Zhou family의 earlier methods는 **historical/contextual reference**로만 둔다.  
primary adjudication rule로 사용하지 않는다.

### 8.6 Required outputs

#### 8.6.1 RECeUS outputs

- `π̂_n`-related quantity
    
- `r̂_n`
    
- RECeUS screening note
    
- RECeUS-AIC result if implemented
    
- threshold note
    

#### 8.6.2 Xie outputs

- test statistic
    
- p-value or bootstrap p-value
    
- contradiction flag
    
- interpretation note
    

#### 8.6.3 Maller–Zhou outputs

- statistic if computed
    
- contextual note
    

### 8.7 Allowed interpretation

- cure-aware modeling을 screening 관점에서 볼 만한가
    
- follow-up adequacy에 대한 contradiction signal이 있는가
    
- older contextual methods와 최신 pragmatic screening이 어떻게 비교되는가
    

### 8.8 Forbidden interpretation

- screening = cure proof
    
- RECeUS or Xie alone = final scientific conclusion
    

---

## 9. Block 6 — Cure decomposition stability / latency fragility diagnostics

### 9.1 Objective

이 블록의 목적은 **cure fraction–latency decomposition이 현재 데이터에서 얼마나 불안정한지 직접 점검하는 것**이다.  
이 블록은 새 headline estimate를 만드는 primary modeling block이 아니라, **Block 1과 Block 3의 cure-aware 결과를 해석하는 guardrail block**으로 둔다.

### 9.2 Core diagnostics

Block 6는 다음 세 분석을 기본 진단 세트로 고정한다.

1. **Leave-one-out / leave-last-k-out influence analysis**
   
2. **Latency family sensitivity**
   
3. **Joint uncertainty between cure fraction and latency scale**

### 9.3 Required outputs

#### 9.3.1 Influence analysis outputs

- dataset별 baseline decomposition summary
   
- late censored / late event removal registry
   
- leave-one-out 영향표
   
- leave-last-k-out 영향표
   
- top influential cases summary

#### 9.3.2 Family sensitivity outputs

- latency family별 fit summary
   
- AIC / log-likelihood comparison
   
- overall risk, cure fraction, latency median comparison
   
- plot-ready family sensitivity figure

#### 9.3.3 Joint uncertainty outputs

- posterior draw table for cure fraction and latency scale
   
- cure fraction vs latency sigma summary
   
- cure fraction vs uncured latency summary
   
- scatter plots or equivalent pairwise uncertainty display

### 9.4 Allowed interpretation

- 어떤 tail case가 decomposition을 크게 흔드는가
   
- family choice에 따라 cure fraction과 latency가 얼마나 달라지는가
   
- posterior 상에서 cure fraction과 latency scale이 trade-off ridge를 이루는가
   
- current data에서 cure decomposition을 strong claim으로 써도 되는가

### 9.5 Forbidden interpretation

- 영향 큰 사례 = data error라고 단정
   
- family sensitivity가 크더라도 하나의 family 결과만 headline으로 고정
   
- posterior trade-off를 discrete cured subgroup의 존재 증명으로 과해석

---

## 10. Cross-block interpretation rules

### 10.1 Separation rule

다음 세 가지는 끝까지 분리해서 쓴다.

- timing difference
    
- follow-up immaturity / tail instability
    
- cure-like heterogeneity

- cure decomposition fragility / latency instability
    

### 10.2 No single-block rule

어느 한 블록도 단독으로 final scientific conclusion을 담당하지 않는다.

### 10.3 Nonregular inference caution

cure vs non-cure 비교에서 naive textbook chi-square LRT를 기본 inferential basis로 쓰지 않는다.  
integrated master specification도 cure-null을 boundary / nonregular problem으로 다루며, bootstrap-calibrated nonregular LRT 또는 descriptive fit contrast만 허용한다.

### 10.4 Bayesian interpretation boundary

Bayesian posterior stabilization만으로 discrete cured subgroup를 확정했다고 말하지 않는다.  
Bayesian 해석은 별도 Bayesian specification의 경계 안에서만 수행한다.

---

## 11. Final synthesis

### 11.1 Question order

최종 결과 취합은 아래 질문 순서로 진행한다.

1. **모델별 probability가 얼마나 다른가?** → Block 1
    
2. **짧은 추적과 긴 추적에서 결과가 얼마나 달라지는가?** → Block 3
    
3. **그 차이가 early/intermediate timing difference 때문인가?** → Block 2
    
4. **late tail이 얼마나 성숙하고 안정적인가?** → Block 4
    
5. **cure model을 screening 관점에서 볼 만한가?** → Block 5

6. **cure fraction–latency decomposition이 얼마나 fragile한가?** → Block 6
    

### 11.2 Final claim logic

최종 문장은 위 여섯 축을 연결해서 쓴다.  
즉,

- Block 1에서 probability divergence를 보이고,
    
- Block 3에서 short vs long follow-up 차이를 직접 보여주고,
    
- Block 2에서 단순 timing difference만으로 환원되지 않는지 확인하고,
    
- Block 4에서 tail maturity를 점검하고,
    
- Block 5에서 cure-model appropriateness screening을 덧붙인 뒤,

- Block 6에서 decomposition fragility와 latency instability를 확인한 뒤,
    

그 결과를 종합해 해석한다.

### 11.3 One-sentence summary

이 spec은 기존 stage numbering을 직접 운영 단위로 쓰지 않고,  
**probability comparison → timing separation → short-vs-long direct comparison → follow-up stability → cure-appropriateness screening → decomposition stability diagnostics**  
의 6-block hierarchy로 분석을 구성한다.
