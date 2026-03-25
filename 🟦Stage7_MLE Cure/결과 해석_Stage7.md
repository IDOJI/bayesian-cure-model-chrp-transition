## 0. QC summary

결론부터 말하면, **업로드된 Stage 7 파일들만 기준으로 볼 때 치명적 계산 실패는 확인되지 않았습니다.**  
즉, “이 결과는 전부 깨졌으니 해석 불가” 수준의 구조적 붕괴는 보이지 않았습니다.

다만, **해석의 신뢰도를 낮추는 중요한 주의점**은 분명히 있습니다.

1. **업로드 세트가 manifest 대비 일부 누락**되어 있습니다.  
    `stage7_subject_predictions.csv`, 시각화 PDF, 일부 bootstrap 상세 파일들이 없어서, 표-그림 완전 대조와 subject-level prediction QC는 끝까지 못 했습니다. 이건 **치명적 오류는 아니지만 QC 범위를 제한하는 문제**입니다.
    
2. **Stage 6 screening flag가 Stage 7에 조인되지 않았습니다.**  
    업로드된 Stage 7 테이블들에서 `stage6_flag_joined = FALSE`가 전부이며, `stage6__screening_note`도 “Stage 6 carry-forward was not joined in this refresh run.”으로 동일합니다. 즉, Stage 7 해석 시 반드시 붙어 있어야 할 **cure appropriateness 맥락**이 빠져 있습니다. 이는 **계산 실패라기보다 해석 체인의 단절**입니다. 프로젝트 master spec은 Stage 6 결과를 Stage 7/8로 carry-forward 하라고 명시합니다.
    
3. **SNU의 일부 cure model은 optimizer 비수렴 경고가 남아 있습니다.**  
    `fit_ok=TRUE`라도 괜찮다고 끝내면 안 되는데, 실제로 SNU의 2개 cure fit은 `convergence_code = 1`, `has_converged_solution = FALSE` 상태에서 “best available finite fit retained”로 남아 있습니다. 이것은 **치명적 실패는 아니지만 수치적 불안정성**입니다.
    
4. **interaction branch의 incidence 계수는 여러 모델에서 경계값/과대불안정 양상**을 보입니다.  
    특히 SNU·merged의 `cure_lnorm`, `cure_llogis` 계열, 그리고 PNU interaction에서 계수와 신뢰구간이 비정상적으로 커집니다. 이건 **모형이 임상 현실을 잘 포착했다기보다, 정보량 대비 파라미터가 과해졌을 가능성**이 큽니다.
    
5. **PNU는 remission이 실제 competing process인데, Stage 7 main backbone은 remission을 censoring으로 처리합니다.**  
    프로젝트 backbone 자체가 `status_num==1`만 event, `status_num in {0,2}`를 censoring으로 두도록 고정되어 있습니다.  
    따라서 PNU에서 보이는 late plateau나 cure-like signal은, **진짜 cure라기보다 remission 처리 방식의 산물일 수 있습니다.** 이건 프로젝트 master spec도 Stage 9에서 따로 재점검하라고 강조합니다.
    
6. **PNU의 horizon support label이 현재 master spec과 일부 어긋납니다.**  
    현재 integrated spec은 PNU 1년 primary, 2년 sensitivity, 3년 이상 projection을 요구하는데, 업로드된 Stage 7 risk summary에서는 PNU 2–5년이 `secondary`로 찍혀 있습니다.  
    이것은 **숫자 계산 실패는 아니고 라벨링/해석 규칙 불일치**입니다.
    

따라서 분기 규칙상 저는 **“major computational/data-structure failure 없음”**으로 판단합니다.  
하지만 **해석은 반드시 보수적으로**, 특히

- PNU의 장기 cure 해석,
    
- interaction model의 incidence coefficient 해석,
    
- SNU 일부 cure fit의 세부 계수 해석
    

에서는 강한 단정이 금물입니다.

---

## 0. Evidence tables

아래 표들은 **업로드된 CSV들을 직접 대조해서 재구성한 QC 표**입니다.  
정확한 수치의 source of truth는 업로드된 Stage 7 CSV들입니다.  
해석 규칙은 project master spec과 data dictionary를 따랐습니다.

---

### 0-1. Data structure / cohort summary

**Source:** `@_stage7_fit_registry.csv`, `@_stage7_ipcw_registry.csv`, `@_stage7_output_manifest.csv`

```text
[업로드된 Stage 7 핵심 파일 구조]

fit_registry                 91 rows
model_registry               12 rows
risk_summary                910 rows
delta_risk                  400 rows
threshold_metrics          3880 rows
ipcw_registry              5240 rows
coefficients                530 rows
bootstrap_model_qc           91 rows
replicate_overview         1000 rows
replicate_failure_by_model    7 rows
stage6_carry_forward_registry 0 rows
output_manifest              19 rows
```

```text
[dataset별 분석 대상 수와 event 구조]

PNU     n=54   transition=12   right_censor=29   remission=13
SNU     n=208  transition=38   right_censor=170  remission=0
merged  n=262  transition=50   right_censor=199  remission=13
```

```text
[모델 row 구조]

PNU:
  benchmark 1
  base 11
  interaction 11

SNU:
  benchmark 1
  base 11
  interaction 11

merged:
  benchmark 1
  base_sitefree 11
  interaction_sitefree 11
  base_siteadjusted 11
  interaction_siteadjusted 11
```

```text
[ipcw_registry 구조 확인]

총 row = 5240
= (PNU 54 + SNU 208 + merged 262) × 10 horizons
```

#### 숫자가 직접 말하는 것

- row 수는 모델 registry, risk summary, IPCW registry가 **서로 산술적으로 잘 맞습니다**.
    
- `threshold_metrics 3880 rows`는 겉보기엔 많아 보이지만, 모델 rows + `treat_all`/`treat_none` benchmark strategy rows가 같이 들어간 결과라 **이상치가 아닙니다**.
    

#### 쉬운 말로 풀면

- 파일이 뒤섞여 중복 저장되었거나, 같은 사람을 여러 번 잘못 쌓아 올린 흔적은 보이지 않습니다.
    
- “행 수가 많다 = 파일이 망가졌다”가 아니라, **모델 수 × horizon × threshold** 구조라서 그렇게 나온 겁니다.
    

#### 그래서 중요한 이유

- 이런 구조가 맞아야, 뒤에 나오는 risk table, delta risk, threshold metrics를 서로 비교할 수 있습니다.
    
- 구조가 틀리면 이후 모든 해석이 무의미해지는데, 여기서는 그 수준의 붕괴는 없습니다.
    

#### 지금 주장할 수 있는 것

- **Stage 7 export layering은 대체로 정상**입니다.
    
- **analytic sample size와 row expansion 규칙은 일관적**입니다.
    

#### 지금 주장하면 안 되는 것

- 업로드되지 않은 `subject_predictions`와 plot PDF까지 완벽히 일치한다고는 아직 못 말합니다.
    
- 전체 Stage 7 산출물이 100% 완전 QC되었다고는 못 말합니다.
    

---

### 0-2. Event and follow-up summary

**Source:** `@_stage7_fit_registry.csv`, `3.Data Dictionary_🇬🇧ENG.md`, `Integrated_Modeling_Master_Specification_English.md`

프로젝트 backbone상 main event는 `status_num == 1`인 transition이고, `status_num in {0,2}`는 censoring입니다. 또한 PNU에만 remission competing-risk 구조가 존재하고, SNU에는 remission이 구조적으로 없습니다.

```text
[main analysis event rule]

event_main   = transition only
censor_main  = right censoring + remission
```

```text
[PNU와 SNU의 outcome 구조 차이]

PNU: transition 12, remission 13
SNU: transition 38, remission 0
```

#### 숫자가 직접 말하는 것

- PNU에서는 transition과 remission 수가 거의 비슷합니다.
    
- SNU에서는 remission이 아예 없습니다.
    

#### 쉬운 말로 풀면

- PNU에서는 “나빠지는 사람”만 있는 게 아니라 “호전/다른 경로로 빠지는 사람”도 꽤 있습니다.
    
- 그런데 main analysis는 remission을 사건이 아니라 censoring으로 다룹니다.
    

#### 그래서 중요한 이유

- PNU에서 late risk가 낮아 보이거나 plateau가 생기면, 그게 **진짜 cure인지**, 아니면 **remission을 censoring 처리한 결과인지** 구별이 어려워집니다.
    
- SNU와 PNU를 같은 event rule로 묶으면 비교는 쉬워지지만, **자연경과 차이와 outcome 정의 차이**가 섞일 수 있습니다.
    

#### 지금 주장할 수 있는 것

- PNU의 cure-like signal은 **반드시 remission sensitivity와 분리해서 읽어야 합니다**.
    
- SNU에서 remission이 없다는 점은 PNU와 구조적으로 다릅니다.
    

#### 지금 주장하면 안 되는 것

- PNU의 plateau를 바로 “cured subgroup 존재”라고 단정하면 안 됩니다.
    
- PNU-SNU 차이를 바로 생물학적 차이라고 해석하면 안 됩니다.
    

---

### 0-3. Follow-up support by horizon

**Source:** `@_stage7_risk_summary.csv`, `Integrated_Modeling_Master_Specification_English.md`

현재 master spec의 supported-horizon rule은 다음과 같습니다.  
PNU는 1년 primary, 2년 sensitivity, 3년 이상 projection; SNU/merged는 1–2년 primary, 3–5년 secondary, 6년 이상 projection입니다.

그런데 업로드된 `risk_summary`의 실제 `horizon_support` 라벨은 다음과 같습니다.

```text
[risk_summary 내부 label]

PNU:
  1 year = primary_supported
  2-5 year = secondary
  6-10 year = projection

SNU:
  1-2 year = primary_supported
  3-5 year = secondary
  6-10 year = projection

merged:
  1-2 year = primary_supported
  3-5 year = secondary
  6-10 year = projection
```

#### 숫자가 직접 말하는 것

- SNU/merged는 spec과 맞습니다.
    
- **PNU 2년이 sensitivity가 아니라 secondary로 라벨링**되어 있습니다.
    

#### 쉬운 말로 풀면

- 숫자 자체보다 “이 숫자를 어느 정도 믿어야 하는지” 붙여 주는 딱지가 조금 어긋나 있습니다.
    

#### 그래서 중요한 이유

- 같은 2년 결과라도 “민감도 확인용”인지 “2차 해석 가능”인지는 해석 강도가 다릅니다.
    
- 특히 follow-up이 짧은 PNU에서는 이 차이가 큽니다.
    

#### 지금 주장할 수 있는 것

- **PNU 2년 이후 해석은 현재 spec 기준으로 더 보수적이어야 합니다.**
    
- 업로드 Stage 7 라벨을 그대로 manuscript language로 쓰면 곤란합니다.
    

#### 지금 주장하면 안 되는 것

- PNU 2–5년 결과를 SNU 2–5년과 같은 해석 무게로 다뤄서는 안 됩니다.
    

---

### 0-4. Model comparison table(s)

**Source:** `@_stage7_risk_summary.csv`, `@_stage7_delta_risk.csv`

아래는 **해석에 꼭 필요한 비교 block**만 뽑은 것입니다.  
모델 ranking이 아니라 **어떤 과학적 패턴이 반복되는지**를 보기 위한 표입니다.

#### 0-4A. PNU base: horizon별 risk 비교 block

```text
[PNU / base / selected horizons]

1 year
  no-cure range:   0.171 ~ 0.223
  cure range:      0.165 ~ 0.243
  no-cure minus cure delta range: about -0.020 ~ +0.002 수준의 작은 차이

5 years
  no-cure range:   대략 0.271 ~ 0.593
  cure range:      대략 0.333 ~ 0.388
  delta range:     음수도 있고 양수도 있음

10 years
  no-cure range:   0.308 ~ 0.865
  cure range:      0.366 ~ 0.440
  delta range:     -0.058 ~ +0.495
```

#### 숫자가 직접 말하는 것

- **1년에서는 cure/no-cure risk 차이가 매우 작습니다.**
    
- **10년에서는 모델에 따라 차이가 폭발적으로 커집니다.**
    

#### 쉬운 말로 풀면

- 짧은 구간에서는 cure를 넣든 안 넣든 예측이 거의 비슷합니다.
    
- 긴 구간으로 갈수록 “모형이 어떻게 꼬리를 그리느냐”에 따라 결과가 완전히 달라집니다.
    

#### 그래서 중요한 이유

- PNU에서 cure fraction이 있어 보이더라도, **관찰로 직접 뒷받침되는 짧은 구간에서는 아직 큰 임상적 차이로 번지지 않습니다**.
    
- 반대로 장기 risk 차이는 대부분 **모형 투영의 산물**일 가능성이 큽니다.
    

#### 지금 주장할 수 있는 것

- PNU에서 “cure model이 단기 임상의사결정을 크게 바꾼다”는 증거는 약합니다.
    
- PNU 장기 cure 해석은 **projection-driven**입니다.
    

#### 지금 주장하면 안 되는 것

- PNU 10년 risk 수치를 “100명 중 몇 명”처럼 직접 관찰치처럼 말하면 안 됩니다.
    

---

#### 0-4B. SNU base: horizon별 risk 비교 block

```text
[SNU / base / selected horizons]

1 year
  no-cure range:   0.051 ~ 0.084
  cure range:      0.070 ~ 0.084

2 years
  no-cure range:   대체로 0.11 ~ 0.16
  cure range:      대체로 0.11 ~ 0.16

5 years
  no-cure range:   0.208 ~ 0.227
  cure range:      0.208 ~ 0.227 정도로 거의 겹침

10 years
  no-cure range:   0.307 ~ 0.404
  cure range:      0.304 ~ 0.351
```

#### 숫자가 직접 말하는 것

- SNU에서는 1–5년 risk가 cure/no-cure 사이에 **거의 비슷**합니다.
    
- 10년에서도 차이는 일부 가족에서만 조금 커집니다.
    

#### 쉬운 말로 풀면

- SNU에서는 cure model을 넣어도 **관찰 범위 안의 절대위험 예측은 거의 안 바뀝니다**.
    

#### 그래서 중요한 이유

- 이것은 “cure model이 무의미하다”가 아니라,  
    **현재 데이터가 단기·중기 위험 예측에서 no-cure와 cure를 뚜렷이 갈라놓지 못한다**는 뜻입니다.
    
- 즉, SNU는 “cure fraction 추정”보다 “관찰 범위 위험 예측의 안정성”이 더 눈에 띕니다.
    

#### 지금 주장할 수 있는 것

- SNU의 1–5년 주장은 상대적으로 robust합니다.
    
- cure model의 핵심 차이는 **장기 tail 해석**에 더 가깝습니다.
    

#### 지금 주장하면 안 되는 것

- SNU에서 cure fraction의 정확한 크기가 확정되었다고 말하면 안 됩니다.
    

---

#### 0-4C. SNU base: cure fraction comparison block

```text
[SNU / base / cure fraction by family]

cure_lnorm          0.0375
cure_llogis         0.0383
cure_weibull        0.4129
cure_exp            0.6004
cure_coxlatency     0.6486
cure_aft_sensitivity 0.6887
```

#### 숫자가 직접 말하는 것

- 같은 SNU 데이터에서 추정 cure fraction이 **약 4%에서 약 69%까지** 뜁니다.
    

#### 쉬운 말로 풀면

- “얼마나 많은 사람이 사실상 장기 non-susceptible인가?”라는 질문에 대해, 모델 가족만 바꿔도 답이 크게 달라집니다.
    

#### 그래서 중요한 이유

- 이것은 **잠재 이질성의 존재 가능성**과 **그 크기의 식별 가능성**을 구분해야 한다는 뜻입니다.
    
- 즉, “cure-like heterogeneity가 전혀 없다”까지는 아니지만, “몇 %가 cured다”는 숫자는 아직 아주 불안정합니다.
    

#### 지금 주장할 수 있는 것

- SNU에서 cure fraction의 **존재 가능성은 논의할 수 있어도, 크기 추정은 family-sensitive**합니다.
    
- supported horizon risk가 거의 안 바뀌는 점을 함께 보면, cure fraction 숫자는 주로 tail parameterization에 좌우됩니다.
    

#### 지금 주장하면 안 되는 것

- “SNU의 cure fraction은 65%다” 같은 단정.
    

---

#### 0-4D. merged base_sitefree vs base_siteadjusted block

```text
[merged / selected comparison]

base_sitefree cure fraction examples:
  cure_lnorm       0.084
  cure_llogis      0.142
  cure_weibull     0.529
  cure_coxlatency  0.624
  cure_exp         0.667

base_siteadjusted:
  family별 수치는 변하지만,
  site adjustment 후에도 cure/no-cure 단기 risk 차이는 크지 않음
  그러나 site coefficient는 여러 latency model에서 일관되게 큼
```

#### 숫자가 직접 말하는 것

- merged에서도 cure fraction은 family-sensitive입니다.
    
- site를 넣어도 단기 risk 예측은 크게 달라지지 않지만, **site effect 자체는 남습니다**.
    

#### 쉬운 말로 풀면

- 합쳐서 돌려도 “센터 차이”가 없어지지 않습니다.
    
- 즉, PNU-SNU 차이는 단순한 tail artifact만은 아닐 수 있습니다.
    

#### 그래서 중요한 이유

- 이것은 cohort difference가 **진짜 timing difference**, **selection/care pathway 차이**, **follow-up structure 차이**를 포함하고 있음을 시사합니다.
    
- 프로젝트 spec도 site effect를 clean treatment effect로 읽지 말라고 합니다.
    

#### 지금 주장할 수 있는 것

- merged 분석에서 site 조정은 필요합니다.
    
- site effect는 cure model을 넣어도 사라지지 않는 중요한 구조 신호입니다.
    

#### 지금 주장하면 안 되는 것

- site coefficient를 치료 효과라고 번역하면 안 됩니다.
    

---

### 0-5. Parameter / coefficient table(s)

**Source:** `@_stage7_coefficients.csv`, `@_stage7_bootstrap_model_qc.csv`

#### 0-5A. interaction instability block

```text
[대표적 불안정 예시]

SNU interaction cure_lnorm / cure_llogis:
  incidence 계수 일부가 수백~수천 단위로 치솟음
  신뢰구간도 비정상적으로 광범위

PNU interaction cure models:
  age_s:sex_num 및 incidence terms가 매우 불안정
  bootstrap success rate도 일부 모델에서 0.887 ~ 0.937 수준
```

#### 숫자가 직접 말하는 것

- interaction을 넣은 cure incidence block은 여러 family에서 사실상 경계 근처 해를 잡고 있습니다.
    

#### 쉬운 말로 풀면

- “남녀에 따라 나이 효과가 달라진다”까지는 말할 수 있어도, 지금 모형은 그걸 정밀하게 말할 만큼 데이터가 충분하지 않을 수 있습니다.
    

#### 그래서 중요한 이유

- interaction branch를 과하게 해석하면, **실제 임상적 신호보다 수치적 과적합**을 읽게 됩니다.
    

#### 지금 주장할 수 있는 것

- interaction model은 주로 **민감도 분석**으로 보는 게 안전합니다.
    
- main message는 base model 중심이 더 방어적입니다.
    

#### 지금 주장하면 안 되는 것

- interaction incidence coefficient 자체를 기전적 발견처럼 해석하면 안 됩니다.
    

---

#### 0-5B. merged site coefficient block

```text
[merged site-adjusted latency signal]

여러 no-cure / cure model에서 site_snu 관련 계수가 일관되게 큼
(모수화에 따라 방향 표시는 달라져도, site 차이 자체는 반복됨)
```

#### 숫자가 직접 말하는 것

- site term은 여러 모델에서 반복적으로 살아남습니다.
    

#### 쉬운 말로 풀면

- PNU와 SNU는 단순히 follow-up 길이만 다른 게 아니라, 사건 발생 시점 구조 자체가 다를 가능성이 큽니다.
    

#### 그래서 중요한 이유

- “PNU가 더 빨리 사건이 난다” 또는 “SNU가 더 완만하다”는 해석이 site proxy를 통해 지지됩니다.
    
- 이는 Stage 4 timing-difference separation과 연결되는 결과입니다. 프로젝트는 timing difference, follow-up immaturity, cure evidence를 분리하라고 요구합니다.
    

#### 지금 주장할 수 있는 것

- cohort/site context는 substantive하게 중요합니다.
    

#### 지금 주장하면 안 되는 것

- 약물 효과나 생물학적 본질 차이라고 단정.
    

---

### 0-6. Prediction or time-specific risk table(s)

**Source:** `@_stage7_threshold_metrics.csv`, `@_stage7_risk_summary.csv`

#### 0-6A. PNU, threshold 0.10, base block

```text
[PNU / base / threshold 0.10]

1 year:
  대부분 모델 false_positive_burden_primary ≈ 1.000
  treat_all net benefit ≈ 0.078
  모델들 net benefit도 거의 비슷
  예외적으로 cure_coxlatency, cure_aft_sensitivity만 약간 낮은 FP burden

2 years:
  사실상 모든 모델이 false_positive_burden_primary = 1.000
  net benefit도 treat_all과 동일 수준
```

#### 숫자가 직접 말하는 것

- PNU에서 10% cutoff를 쓰면 1–2년에는 거의 다 high risk로 분류됩니다.
    

#### 쉬운 말로 풀면

- 이 문턱값에서는 모델이 환자를 잘 가려내는 게 아니라, **거의 다 위험하다고 찍어 버립니다**.
    

#### 그래서 중요한 이유

- threshold-based clinical usefulness는 “모델이 있느냐”보다 “어떤 threshold를 쓰느냐”에 더 크게 좌우됩니다.
    
- PNU에서는 저 threshold에서 cure/no-cure 비교가 임상적으로 크게 남지 않습니다.
    

#### 지금 주장할 수 있는 것

- PNU의 early horizon에서 10% threshold는 너무 느슨해서 decision utility가 낮습니다.
    

#### 지금 주장하면 안 되는 것

- 이걸 근거로 cure model이 무가치하다고 일반화하면 안 됩니다. threshold 문제일 수도 있습니다.
    

---

#### 0-6B. SNU, threshold 0.10 and 0.20, base block

```text
[SNU / base]

threshold 0.10, 1 year:
  false_positive_burden_primary가 모델별로 0.000 ~ 0.141 정도로 퍼짐

threshold 0.10, 5 year:
  거의 대부분 1.000에 가까워짐

threshold 0.20, 1-2 year:
  false_positive burden 거의 0에 가까움

threshold 0.20, 5 year:
  모델 간 spread가 남아 있음
  예: cure_coxlatency는 더 낮고, 일부 no-cure/parametric은 1.000에 가까움
```

#### 숫자가 직접 말하는 것

- SNU에서는 early horizon에서는 threshold choice에 따라 모델 차이가 보입니다.
    
- 시간이 길어질수록 많은 모델이 treat-all처럼 됩니다.
    

#### 쉬운 말로 풀면

- 초반 1–2년에는 모델이 불필요한 high-risk labeling을 어느 정도 줄일 여지가 있습니다.
    
- 그러나 긴 horizon으로 갈수록, 특히 낮은 threshold에서는 다들 비슷하게 “많이 위험하다”고 하게 됩니다.
    

#### 그래서 중요한 이유

- 이 프로젝트의 실제 임상 메시지는 **장기 tail 숫자 싸움**보다, **supported horizon에서 불필요한 high-risk labeling을 얼마나 줄이느냐**에 더 잘 맞습니다. master spec도 false-positive burden과 clinical usefulness를 핵심 목표로 둡니다.
    

#### 지금 주장할 수 있는 것

- SNU에서는 1–2년, 적절한 threshold에서 모델 차이가 임상적으로 의미 있을 수 있습니다.
    

#### 지금 주장하면 안 되는 것

- 10년 threshold 결과를 직접 관찰 근거처럼 쓰면 안 됩니다.
    

---

### 0-7. Diagnostics / convergence / assumption-check table(s)

**Source:** `@_stage7_fit_registry.csv`, `@_stage7_bootstrap_model_qc.csv`, `@_stage7_replicate_overview.csv`, `@_stage7_replicate_failure_by_model.csv`

```text
[fit registry 요약]

fit_component_success      91 / 91 TRUE
prediction_component_success 91 / 91 TRUE
overall_status:
  ok                      89
  ok_with_warnings         2

warning/convergence problem rows:
  SNU / base / cure_llogis
  SNU / interaction / cure_lnorm
```

```text
[bootstrap replicate overview]

1000 replicates 중:
  NO_FAILURE = 799
  mean failed model rows per replicate = 0.249
  mean overall_success_rate ≈ 0.997
  max failed rows in a replicate = 3
```

```text
[failure_by_model]

PNU interaction cure_coxlatency          n_replicates = 113
PNU interaction cure_aft_sensitivity     n_replicates = 63
SNU interaction cure_aft_sensitivity     n_replicates = 25
PNU base cure_aft_sensitivity            n_replicates = 16
PNU base cure_coxlatency                 n_replicates = 12
merged interaction_siteadjusted cure_aft n_replicates = 10
merged interaction_sitefree cure_aft     n_replicates = 10
```

```text
[bootstrap_model_qc qc_flag]

stable          84
usable_for_ci    7
```

#### 숫자가 직접 말하는 것

- 전반적 계산은 잘 돌아갔습니다.
    
- 불안정성은 특정 cure families, 특히 `cure_aft_sensitivity`, `cure_coxlatency`, interaction branch에 집중됩니다.
    

#### 쉬운 말로 풀면

- 전체 시스템이 무너진 건 아니고, **몇몇 복잡한 cure model만 자주 비틀거립니다**.
    

#### 그래서 중요한 이유

- 이런 패턴은 “질병이 정말 그렇게 복잡하다”기보다, **데이터가 그 복잡함을 버티기엔 부족하다**는 신호일 때가 많습니다.
    
- 특히 AFT-latency sensitivity와 Cox-latency cure는 tail/identifiability 영향을 많이 받습니다. Cox와 AFT cure가 충분/불충분 follow-up 영역 밖에서 크게 달라질 수 있다는 문헌과도 맞닿아 있습니다.
    

#### 지금 주장할 수 있는 것

- **base model 중심 해석은 가능**합니다.
    
- 그러나 **unstable branch는 결론의 주축이 아니라 sensitivity**로 둬야 합니다.
    

#### 지금 주장하면 안 되는 것

- unstable branch에서 나온 cure fraction, interaction effect, tail risk를 “결정적 증거”로 쓰면 안 됩니다.
    

---

### 0-8. Cohort-specific comparison table(s)

**Source:** `@_stage7_risk_summary.csv`, `@_stage7_threshold_metrics.csv`, `Integrated_Modeling_Master_Specification_English.md`

```text
[cohort pattern summary]

PNU:
  event earlier, remission 존재, follow-up 짧음
  short horizon에서는 cure/no-cure 차이 작음
  long horizon에서는 model tail dependence 큼

SNU:
  remission 없음, follow-up 더 길음
  1-5년 risk는 cure/no-cure 거의 유사
  cure fraction 크기 자체는 family-sensitive

merged:
  site effect 남음
  site adjustment 후에도 cohort 차이 신호 지속
```

#### 숫자가 직접 말하는 것

- PNU와 SNU는 단지 표본 크기 차이가 아니라 **event process와 follow-up structure가 다른 두 코호트**로 보입니다.
    
- merged는 이를 평균내서 지워 주지 않습니다.
    

#### 쉬운 말로 풀면

- PNU는 “빨리 일이 벌어지고, 호전(remission)도 따로 보이는 짧은 코호트”
    
- SNU는 “더 길게 지켜본, 비교적 단순한 transition-vs-censoring 코호트”
    
- merged는 둘을 합쳤지만 차이가 사라지지 않은 셈입니다.
    

#### 그래서 중요한 이유

- 이 패턴은 “어느 모델이 이겼나”보다 훨씬 중요합니다.
    
- 왜냐하면 cure-like signal이 **질병 자연사** 때문인지, **치료/관리 맥락** 때문인지, **추적 구조 차이** 때문인지 가르는 실마리가 되기 때문입니다.
    

#### 지금 주장할 수 있는 것

- cohort structure 자체가 substantive합니다.
    
- merged analysis는 유용하지만 separate-cohort analysis를 대체하지 못합니다.
    

#### 지금 주장하면 안 되는 것

- 하나의 공통 biological process가 PNU와 SNU에 동일하게 작동한다고 단정.
    

---

## A. Scientific question reframed

이 분석의 진짜 질문은  
**“어떤 cure model이 AIC가 제일 좋은가?”**가 아닙니다.

더 정확한 질문은 다음입니다.

> **이 질병/상태 전이 과정에서 실제로 장기 비감수성(cure-like heterogeneity)이 있는지, 그리고 그 가정을 허용했을 때 임상적으로 중요한 시점에서 불필요한 고위험 분류를 줄일 만큼 예측이 달라지는지**를 묻는 것입니다.

그리고 이 질문에 답할 때는 반드시 세 가지를 분리해야 합니다.

1. **timing difference**
    
2. **follow-up immaturity / tail instability**
    
3. **cure-like heterogeneity**
    

프로젝트 master spec도 이 셋을 절대 섞지 말라고 못 박고 있습니다.

---

## B. Bottom-line integrated interpretation

- **큰 계산 붕괴는 없어서 해석 자체는 가능합니다.** 다만 해석의 무게 중심은 Stage 7 숫자 그 자체보다, 이 숫자가 어떤 follow-up 구조와 outcome 정의 위에서 나왔는지에 둬야 합니다.
    
- **PNU에서 보이는 cure-like signal은 특히 취약합니다.** 이유는 follow-up이 짧고, remission이 실제 competing process인데 main analysis에서는 censoring 처리되기 때문입니다. 그래서 PNU 장기 cure 해석은 관찰 사실이라기보다 model-dependent projection에 가깝습니다.
    
- **SNU는 1–5년 절대위험 예측이 더 안정적이지만, cure fraction의 ‘크기’는 매우 불안정합니다.** 즉, “cure-like heterogeneity가 아예 없다”고 하긴 어렵지만, “몇 %가 cured다”는 숫자는 latency family에 크게 흔들립니다.
    
- **merged 분석에서 site effect가 남는 것은 가장 중요한 과학적 힌트 중 하나입니다.** 이것은 PNU–SNU 차이가 단지 꼬리 추정의 잡음이 아니라, timing, selection, care pathway, follow-up structure의 차이를 반영할 수 있음을 시사합니다. master spec은 site를 clean treatment effect로 해석하지 말라고 합니다.
    
- **실제 출판 가치가 높은 메시지는 “cure fraction 추정치의 크기”보다, supported horizon에서 cure model이 false-positive burden을 얼마나 줄이는가, 그리고 그 패턴이 cohort와 follow-up maturity에 따라 어떻게 달라지는가”입니다.**
    

---

## C. Statistical interpretation

### 1) 무엇이 수치적으로 지지되는가

- 구조적 export는 대체로 정합적입니다.
    
- 전체 fit는 대부분 성공했고, bootstrap failure도 일부 특정 cure branches에 집중됩니다.
    
- supported horizon에서 PNU와 SNU의 risk 예측은 해석 가능하며, 특히 SNU 1–2년은 상대적으로 안정적입니다.
    

### 2) 무엇이 수치적으로 약한가

- SNU와 merged의 cure fraction 크기는 family-sensitive합니다.
    
- PNU 장기 risk와 delta risk는 tail extrapolation 의존성이 매우 큽니다.
    
- interaction branch는 다수 모델에서 경계해/과대계수 문제가 있습니다.
    
- Stage 6 carry-forward 부재 때문에, Stage 7 cure 결과를 “screening 통과한 cure evidence”처럼 바로 읽으면 안 됩니다.
    

### 3) p-value보다 더 중요한 점

- 지금 핵심은 “유의/비유의”가 아니라,
    
    - supported horizon에서 차이가 실제로 얼마나 큰지,
        
    - 그 차이가 threshold-based classification을 바꾸는지,
        
    - 그 차이가 cohort와 outcome definition의 구조적 차이와 어떻게 얽히는지입니다.
        

---

## D. Clinical interpretation

여기서 임상 맥락 자체는 업로드 파일에 명시적으로 적혀 있지 않습니다.  
다만 프로젝트 전체 문맥상 `transition`은 불리한 상태 전이, `remission`은 경쟁적 호전 경로로 보입니다. 이 가정 아래에서 해석하면:

### 1) PNU가 말해 주는 것

PNU는 사건이 더 이르게 발생하고, remission도 함께 존재합니다.  
이 말은 “일찍 나빠지는 사람”과 “다른 경로로 빠지는 사람”이 같이 섞여 있다는 뜻입니다.

쉽게 말하면,

- 어떤 사람은 빨리 transition으로 가고
    
- 어떤 사람은 remission으로 빠져서
    
- 뒤에 남은 곡선이 평평해져 보일 수 있습니다.
    

그래서 **PNU의 plateau는 ‘완치된 하위군’의 증거일 수도 있지만, equally plausible하게는 ‘경쟁적 호전 경로를 censoring 처리한 결과’일 수도 있습니다.**

### 2) SNU가 말해 주는 것

SNU는 remission 없이 transition-vs-censoring 구조입니다.  
그래서 outcome definition이 더 단순하고, 1–5년 absolute risk는 cure/no-cure 간 큰 차이가 없습니다.

쉽게 말하면,

- SNU에서는 “장기적으로 아주 특별한 cured subgroup이 있느냐”보다,
    
- **짧고 중간 정도 기간에서 누가 실제로 transition할 가능성이 높은가**가 더 잘 보입니다.
    

### 3) merged가 말해 주는 것

merged에서 site effect가 남는다는 것은,  
이 시스템이 단일한 자연경과 한 줄로 설명되기보다  
**코호트가 만들어진 방식, 추적 구조, 관리 체계, 환자 선별**에 의해 달라질 수 있음을 뜻합니다.

즉,

- 치료 맥락
    
- care pathway
    
- referral structure
    
- baseline severity mix
    
- follow-up intensity
    

중 하나 이상이 PNU와 SNU를 갈라놓고 있을 가능성이 큽니다.

---

## E. Model-assumption-to-clinical-reality mapping

### 1) KM benchmark

- **가정:** 개별화 없이 관찰된 누적위험을 그대로 benchmark로 사용
    
- **임상적 타당성:** 높음. 가장 덜 가정적입니다.
    
- **한계:** 환자별 위험 차이를 반영하지 못합니다.
    
- **관찰 패턴과 연결:** supported horizon에서 cure/no-cure 차이가 작을 때, KM은 “복잡한 모형이 실제로 임상적으로 얼마나 더 필요한가”를 묻는 좋은 기준점입니다.
    

### 2) no-cure models

- **가정:** 결국 모두 event risk를 가진다
    
- **임상적 타당성:** transition risk가 장기적으로 계속 남아 있거나, 진짜 cure-like subgroup가 작을 때는 더 덜 위배될 수 있습니다.
    
- **관찰 패턴과 연결:** SNU 1–5년에서 cure/no-cure risk가 거의 비슷하다는 것은, **관찰 가능한 기간 안에서는 no-cure assumption이 크게 틀리지 않을 수 있음**을 뜻합니다.
    

### 3) frequentist cure models

- **가정:** population이 susceptible subgroup와 non-susceptible subgroup의 mixture다
    
- **임상적 타당성:** 장기 plateau가 진짜 latent heterogeneity를 반영할 때는 plausible합니다.
    
- **관찰 패턴과 연결:**
    
    - PNU에서는 remission censoring 때문에 이 가정이 과도하게 유리해질 수 있습니다.
        
    - SNU에서는 cure fraction “존재 가능성”은 있지만 “크기”는 family-dependent합니다.
        
- **무엇을 바꾸는가:** 장기 risk tail과 threshold-based FP burden을 바꿀 수 있습니다.
    
- **무엇을 아직 못 말하는가:** “실제 cured subgroup 크기”의 정확한 추정.
    

### 4) Cox-latency vs AFT-latency cure

Parsa & Van Keilegom은 follow-up이 불충분하거나 support가 covariate에 따라 달라질 때 Cox mixture cure와 AFT mixture cure가 충분히 다른 tail을 줄 수 있다고 설명합니다.  
현재 결과에서도 `cure_aft_sensitivity`와 `cure_coxlatency`는 bootstrap failure가 상대적으로 많고, PNU interaction에서 특히 불안정합니다.

쉽게 말하면,

- “event timing 구조를 비율로 볼지”
    
- “시간 자체를 늘고 줄고 하는 식으로 볼지”
    

에 따라 cure 해석이 크게 바뀔 수 있습니다.  
지금 데이터는 그 차이를 안정적으로 가려낼 만큼 충분히 강하지 않은 듯합니다.

---

## F. Why the observed pattern matters

가장 publishable한 포인트는 이겁니다.

> **supported horizon에서는 cure/no-cure 모델이 절대위험을 크게 다르게 예측하지 않지만, 장기 tail로 갈수록 cure fraction의 크기와 risk 차이는 모델 family와 cohort structure에 크게 의존한다. 따라서 이 데이터에서 진짜 핵심은 “몇 %가 cured인가”보다, “어떤 horizon/threshold에서 불필요한 고위험 분류를 줄일 수 있는가”이다.**

이게 중요한 이유는,

- 단순 model ranking을 넘어서
    
- 질병의 자연경과가 단일 경로가 아닐 가능성,
    
- follow-up maturity의 제한,
    
- PNU remission 구조의 왜곡 가능성,
    
- site/care-pathway 차이
    

를 한꺼번에 설명할 수 있기 때문입니다.

---

## G. Alternative explanations and threats to interpretation

### 강하게 지지되는 설명

1. **follow-up maturity 차이**
    
    - PNU는 짧고, SNU는 더 깁니다.
        
    - 장기 tail 차이의 상당 부분은 이 설명으로 충분합니다.
        
2. **outcome definition asymmetry**
    
    - PNU에는 remission competing process가 있고, SNU는 없습니다.
        
    - main analysis가 remission을 censoring 처리하므로 PNU cure-like signal이 과장될 수 있습니다.
        
3. **site/context heterogeneity**
    
    - merged에서 site effect가 반복적으로 남습니다.
        
    - 단순 잡음보다 구조 차이일 가능성이 큽니다.
        

### 그럴듯하지만 아직 덜 확실한 설명

4. **진짜 latent cured subgroup 존재**
    
    - 일부 plateau와 일부 cure fits는 이를 시사합니다.
        
    - 하지만 cure fraction 크기가 family-sensitive라 강하게 못 박기는 어렵습니다.
        
5. **sex-by-age interaction에 의한 subgroup-specific cure pattern**
    
    - interaction branch가 이를 암시할 수는 있습니다.
        
    - 그러나 현재 계수 불안정성이 커서 직접 주장하기 어렵습니다.
        

### 더 speculative한 설명

6. **특정 치료 효과**
    
    - 현재 데이터/메타데이터만으로는 site를 treatment effect로 해석할 수 없습니다. master spec도 금지합니다.
        

### 무엇이 이 설명들을 가를까

- Stage 6 screening flag를 Stage 7에 다시 조인
    
- Stage 9 remission sensitivity / competing risk
    
- common-window 1년, 2년 재비교
    
- merged site-adjusted와 site-free의 threshold metrics 정밀 비교
    
- interaction reduced model 또는 shrinkage model
    

---

## H. Manuscript-ready discussion paragraph

The Stage 7 results suggest that allowing cure-like heterogeneity does not substantially alter short-horizon transition-risk estimates in the supported observation range, particularly in SNU and in the merged analyses, where cure and no-cure models yielded largely overlapping 1–5 year risks. In contrast, longer-horizon estimates were highly sensitive to latency-family specification, with cure-fraction estimates ranging from near-zero to large values despite similar supported-horizon risks, indicating that the magnitude of the putative cured subgroup is not well identified by the current data. This instability was most pronounced in PNU, where shorter follow-up and the presence of remission as a competing clinical pathway complicate interpretation, especially because the main backbone treats remission as censoring. The persistence of site effects in the merged models further suggests that between-cohort differences reflect more than simple tail artefact, and are likely to encode differences in timing structure, care pathways, selection, or follow-up design rather than a clean treatment effect. Collectively, these findings argue against framing the primary conclusion as a model-ranking exercise; instead, the more defensible interpretation is that cure-model utility in this setting lies chiefly in its potential to modify false-positive burden and clinical decision-making at supported horizons, whereas long-horizon cure-fraction estimates should be regarded as provisional and strongly model-dependent.

---

## I. Next-step analyses

1. **Stage 6 flag를 Stage 7 테이블에 재조인해서 다시 export**
    
    - 가장 먼저 해야 합니다.
        
    - 지금은 cure appropriateness screening 맥락이 빠진 채 Stage 7을 읽고 있습니다.
        
2. **PNU 중심 Stage 9 remission competing-risk / multi-state sensitivity**
    
    - PNU의 cure-like signal이 remission censoring artifact인지 가르는 데 가장 중요합니다.
        
3. **common-window 1년, 2년에서 PNU-SNU 직접 비교 강화**
    
    - 이 구간이 가장 관찰로 지지됩니다.
        
    - timing difference가 초기부터 있는지 명확히 보여 줄 수 있습니다.
        
4. **interaction cure model 축소 또는 shrinkage**
    
    - 현재 interaction incidence block은 너무 불안정합니다.
        
    - 단순화된 interaction 또는 regularized approach가 더 낫습니다.
        
5. **threshold 0.10 외에 0.15, 0.20 중심의 DCA/FP burden 재강조**
    
    - 낮은 threshold에서는 많은 모델이 treat-all처럼 수렴합니다.
        
    - 임상적으로 더 구별력 있는 threshold를 앞세워야 합니다.
        

---

## J. 핵심 숫자의 쉬운 해석

### 1) SNU cure fraction이 0.037에서 0.689까지 흔들린다

1. **문자 그대로 뜻하는 것**  
    같은 SNU 데이터에서도 cure family를 바꾸면 추정 cure fraction이 약 4%에서 69%까지 달라집니다.
    
2. **쉬운 말로 바꾸면**  
    “장기적으로 사실상 event를 겪지 않을 사람 비율”을 아직 정확히 못 박기 어렵습니다.
    
3. **관찰 데이터가 얼마나 직접 뒷받침하는지**  
    supported horizon의 risk는 직접 뒷받침되지만, cure fraction 자체는 tail extrapolation과 family choice에 더 의존합니다.
    
4. **그래서 어느 정도까지 믿어도 되는지**  
    “존재 가능성” 정도는 논의 가능하지만, “정확한 비율”은 provisional입니다.
    
5. **임상적/과학적으로 요구하는 해석**  
    cure fraction 숫자 자체보다, supported horizon에서 decision metrics가 얼마나 달라지는지에 더 무게를 둬야 합니다.
    

### 2) PNU 1년 cure/no-cure risk 차이가 작다

1. **문자 그대로 뜻하는 것**  
    PNU 1년 risk는 cure model과 no-cure model 사이 차이가 매우 작습니다.
    
2. **쉬운 말로 바꾸면**  
    짧은 기간 예측에서는 cure를 넣어도 환자별 위험 해석이 크게 달라지지 않습니다.
    
3. **관찰 데이터가 얼마나 직접 뒷받침하는지**  
    1년은 PNU에서 primary-supported horizon입니다.
    
4. **그래서 어느 정도까지 믿어도 되는지**  
    이 부분은 비교적 믿을 수 있습니다.
    
5. **임상적/과학적으로 요구하는 해석**  
    PNU에서 cure model의 가치는 “당장 1년 위험을 뒤엎는 것”보다는, 장기 tail 가정을 다루는 데 있습니다.
    

### 3) PNU 10년 delta risk가 -0.058에서 +0.495까지 벌어진다

1. **문자 그대로 뜻하는 것**  
    no-cure minus cure 장기 risk 차이가 모델에 따라 음수도, 큰 양수도 됩니다.
    
2. **쉬운 말로 바꾸면**  
    긴 기간에서는 모델이 꼬리를 어떻게 그리느냐에 따라 답이 뒤집힙니다.
    
3. **관찰 데이터가 얼마나 직접 뒷받침하는지**  
    PNU 10년은 거의 projection입니다. 직접 관찰 뒷받침이 매우 약합니다.
    
4. **그래서 어느 정도까지 믿어도 되는지**  
    숫자 자체는 “모형 기반 시나리오”로만 봐야 합니다.
    
5. **임상적/과학적으로 요구하는 해석**  
    PNU 장기 cure 주장은 hypothesis-generating 수준이 더 적절합니다.
    

### 4) bootstrap에서 1000회 중 799회가 무실패지만, 실패는 특정 cure families에 몰린다

1. **문자 그대로 뜻하는 것**  
    전체 계산은 안정적이지만, 실패가 일부 복잡한 cure models에 집중됩니다.
    
2. **쉬운 말로 바꾸면**  
    시스템 전체가 불안정한 게 아니라, 몇몇 민감한 cure model만 자주 흔들립니다.
    
3. **관찰 데이터가 얼마나 직접 뒷받침하는지**  
    replicate_overview와 failure_by_model에서 직접 확인됩니다.
    
4. **그래서 어느 정도까지 믿어도 되는지**  
    main/base branch는 꽤 믿을 수 있지만, 불안정 모델은 민감도 분석으로 제한해야 합니다.
    
5. **임상적/과학적으로 요구하는 해석**  
    “복잡한 cure mechanism 발견”보다 “정보량 부족/식별성 취약” 해석이 더 방어적입니다.
    

### 5) PNU는 transition 12명, remission 13명이다

1. **문자 그대로 뜻하는 것**  
    PNU에서는 transition event 수와 remission 수가 비슷합니다.
    
2. **쉬운 말로 바꾸면**  
    나빠지는 경로와 다른 방향으로 빠지는 경로가 둘 다 큽니다.
    
3. **관찰 데이터가 얼마나 직접 뒷받침하는지**  
    데이터 구조 자체가 그렇습니다.
    
4. **그래서 어느 정도까지 믿어도 되는지**  
    매우 직접적인 구조적 사실입니다.
    
5. **임상적/과학적으로 요구하는 해석**  
    PNU cure 해석 전에 remission sensitivity는 사실상 필수입니다.
    

---

## K. 이 결과가 독자에게 실제로 요구하는 해석과 주의점

- **PNU는 1년 결과를 중심으로 읽고**, 2년은 민감도, 그 이후는 projection으로 다루는 것이 안전합니다.
    
- **SNU와 merged는 1–2년을 주력 메시지로 쓰고**, 3–5년은 보조, 6년 이상은 장기 추정으로 명확히 낮춰야 합니다.
    
- **장기 cure fraction 숫자를 소수점까지 비교하지 말아야 합니다.**
    
- **PNU의 cure-like signal은 remission handling의 영향을 분리해서 읽어야 합니다.**
    
- **site effect는 치료 효과가 아니라 care-pathway / selection / follow-up proxy로 해석해야 합니다.**
    
- **threshold-based clinical usefulness는 early supported horizon에서 더 중요하게 보고**, 긴 horizon의 DCA 차이는 projection으로 취급해야 합니다.
    
- **interaction cure model 계수는 설명용이 아니라 민감도 분석용**에 가깝습니다.
    

---

## L. 한 줄 요약

**이번 Stage 7 결과는 “진짜 cure가 몇 %다”를 확정해 주는 분석이라기보다, 코호트 구조와 follow-up 차이를 분리해 보면 short horizon에서는 cure/no-cure 차이가 크지 않고, long horizon의 cure-like 신호는 특히 PNU에서 remission 처리와 tail projection에 크게 흔들린다는 점을 보여 줍니다.**