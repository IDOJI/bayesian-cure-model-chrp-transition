## **QC 요약**

먼저 결론부터 말하면, **중대한 계산적 실패나 자료구조 파손은 확인되지 않았습니다**.

Stage 4는 애초에 supported_short_horizon_contrast, late_horizon_tail_contrast, piecewise_site_effect 사이의 **count alignment**와 **reporting flag alignment**가 맞지 않으면 중단되도록 짜여 있고, 실제로 산출물이 생성되어 있습니다. 또한 metadata에는 이 단계의 역할이 **“cure 해석 전 timing difference를 분리하는 단계”**라고 명시되어 있고, site는 **치료효과 그 자체가 아니라 broader treatment context, care pathway, selection, follow-up structure의 proxy**로 해석하라고 되어 있으며, remission은 본 단계의 main analysis에서 **censoring**으로 처리된다고 적혀 있습니다.     

  

다만 **작은 주의점**은 있습니다.

- restricted Cox와 piecewise Cox의 **세부 optimizer convergence code/message**는 업로드된 CSV/PDF에 직접 노출되어 있지 않아, 그 부분은 완전한 외부 감사를 할 수 없습니다.
    
- 이 단계는 의도적으로 **timing difference**만 분리해 보는 단계라서, **follow-up immaturity**, **cure-fraction evidence**, **remission competing process**와는 분리해서 읽어야 합니다. 이 구분은 프로젝트 master specification에서도 강하게 요구됩니다. 
    

  

또 하나 매우 중요합니다.

**질병명, “transition”의 정확한 임상 정의, 치료 배정/노출 정보, 결과 판정 프로토콜**은 이번 Stage 4 파일만으로는 확인되지 않습니다. 따라서 아래의 임상적 해석은 **“transition이 두 코호트에서 동일한 임상 사건으로 비교 가능하다”**는 전제 아래에서만 강하게 성립합니다. 이 전제가 흔들리면 해석도 달라집니다.

---

# **0. Evidence tables**

  

## **0-1. Data structure / cohort summary**

  

**Source: stage4_dataset_summary.csv**

|**Dataset**|**Rows**|**Unique subjects**|**Transition**|**Remission**|**Right censoring**|**Median follow-up, y**|**Max follow-up, y**|**Person-years**|**Primary-supported horizons**|
|---|---|---|---|---|---|---|---|---|---|
|PNU|54|54|12|13|29|1.11|2.67|54.71|1|
|SNU|208|208|38|0|170|3.80|10.11|735.79|2|
|merged|262|262|50|13|199|2.07|10.11|790.50|2|

**숫자가 직접 말하는 것**

PNU는 표본이 작고 추적기간이 짧습니다. SNU는 표본이 크고 최대 추적이 약 10년까지 있습니다. merged는 단순히 PNU와 SNU를 더한 구조로 보이며, row 수와 subject 수가 각 데이터셋에서 정확히 일치합니다.

  

**쉬운 말로 풀면**

PNU는 “빨리 무슨 일이 생겼는지”는 어느 정도 볼 수 있지만, 오래 가는 꼬리는 거의 못 봅니다. SNU는 오래 가는 꼬리까지 볼 수 있지만 뒤로 갈수록 사람 수가 줄어듭니다.

  

**그래서 중요한 이유**

이 차이 때문에 PNU와 SNU를 같은 해상도로 5년, 10년까지 비교하면 공정하지 않습니다. **초기 1–2년 비교**와 **후반 tail 해석**을 분리해야 합니다.

  

**지금 주장할 수 있는 것**

- export layering 오류나 subject-row 중복은 보이지 않습니다.
    
- merged는 의도된 pooled dataset이지, site-specific stacked export가 중복 삽입된 흔적은 없습니다.
    
    실제로 merged의 subject 수 262는 54 + 208과 정확히 같습니다.
    

  

**지금 주장하면 안 되는 것**

- “merged가 PNU와 SNU를 균등하게 대표한다”
    
- “PNU와 SNU의 5–10년 tail을 동등한 질로 비교할 수 있다”
    

---

## **0-2. Event and follow-up support by horizon**

  

**Source: stage4_risk_trajectory.csv**

|**Dataset**|**Horizon, y**|**At risk beyond horizon**|**Observed KM risk**|**Support tier**|**해석 consequence**|
|---|---|---|---|---|---|
|PNU|1|28|17.0%|primary_supported|직접 관찰로 해석 가능|
|PNU|2|1|29.9%|sensitivity|숫자는 있으나 매우 불안정|
|PNU|5|0|관찰 불가|projection|직접 관찰 해석 불가|
|PNU|10|0|관찰 불가|projection|직접 관찰 해석 불가|
|SNU|1|163|7.3%|primary_supported|직접 관찰로 해석 가능|
|SNU|2|134|12.2%|primary_supported|직접 관찰로 해석 가능|
|SNU|5|64|21.1%|secondary|어느 정도 관찰 기반이나 tail caution 필요|
|SNU|10|1|32.6%|projection|사실상 tail 끝점, 매우 불안정|
|merged|1|191|9.2%|primary_supported|직접 관찰로 해석 가능|
|merged|2|135|15.3%|primary_supported|직접 관찰로 해석 가능|
|merged|5|64|23.9%|secondary|겉보기엔 pooled지만 실질적으론 SNU tail 중심|
|merged|10|1|35.1%|projection|사실상 SNU 한 명의 tail 끝점|

**숫자가 직접 말하는 것**

PNU는 2년 시점에 이미 horizon 너머에 남아 있는 사람이 1명뿐입니다. 반면 SNU는 2년까지 134명, 5년에도 64명이 남아 있습니다. merged의 5년, 10년 risk-set도 SNU와 완전히 같은 64명, 1명입니다.

  

**쉬운 말로 풀면**

PNU는 2년만 넘어가도 “앞으로 어떻게 되는지”를 볼 사람이 거의 없습니다. merged의 5년 이후 숫자는 이름만 merged일 뿐, 실제로는 거의 SNU 꼬리입니다.

  

**그래서 중요한 이유**

이것이 바로 **late-tail instability**입니다.

5–10년 차에서 merged 곡선이 보인다고 해서 “두 코호트를 함께 본 근거가 충분하다”는 뜻이 아닙니다. 사실상 **SNU가 거의 전부를 대표**합니다.

  

**지금 주장할 수 있는 것**

- PNU는 1년 해석이 중심이고, 2년은 민감도 수준입니다.
    
- SNU는 1–2년이 가장 믿을 만하고, 5년은 보조적, 10년은 projection에 가깝습니다.
    
- merged의 5년 이후는 “pooled long-term evidence”라기보다 “SNU-dominant tail”입니다.
    

  

**지금 주장하면 안 되는 것**

- “PNU와 SNU의 10년 위험이 진짜로 다르다”
    
- “merged 10년 위험 35.1%는 두 코호트의 공통 장기위험이다”
    

---

## **0-3. Short-horizon observed risk comparison block**

  

### **0-3a. 1-year observed block**

  

**Source: stage4_supported_short_horizon_contrast.csv**

|**Row**|**Estimate**|**95% interval**|**At risk beyond horizon**|**Events by horizon**|**Reporting note**|
|---|---|---|---|---|---|
|PNU 1y observed risk|17.0%|6.0% to 28.2%|28|8|stable|
|SNU 1y observed risk|7.3%|3.6% to 11.8%|163|14|stable|
|merged 1y pooled observed risk|9.2%|5.4% to 12.7%|191|22|stable|
|PNU minus SNU 1y observed difference|9.7%p|-1.3%p to 22.3%p|191|22|stable|

### **0-3b. 2-year observed block**

  

**Source: stage4_supported_short_horizon_contrast.csv**

|**Row**|**Estimate**|**95% interval**|**At risk beyond horizon**|**Events by horizon**|**Reporting note**|
|---|---|---|---|---|---|
|PNU 2y observed risk|29.9%|15.3% to 45.2%|1|12|bootstrap unstable|
|SNU 2y observed risk|12.2%|7.5% to 17.4%|134|22|stable|
|merged 2y pooled observed risk|15.3%|10.4% to 20.0%|135|34|stable|
|PNU minus SNU 2y observed difference|17.7%p|1.9%p to 33.1%p|135|34|bootstrap unstable|

**숫자가 직접 말하는 것**

1년에서는 PNU가 SNU보다 약 9.7%p 높습니다.

2년에서는 격차가 약 17.7%p로 더 커 보입니다.

  

**쉬운 말로 풀면**

1년 안에 보면 PNU는 100명 중 약 17명, SNU는 약 7명이 transition을 겪은 셈입니다. PNU가 더 빨리 사건에 도달합니다.

2년에서는 격차가 더 커지지만, PNU 2년은 실제로 horizon 너머에 남아 있는 사람이 1명뿐이라 숫자가 흔들립니다.

  

**그래서 중요한 이유**

“PNU가 더 나쁘다”는 신호는 **후반 tail에서 갑자기 생긴 게 아니라**, 이미 **공통 1년 창**에서 시작됩니다.

이건 아주 중요합니다. 왜냐하면 나중 꼬리의 불안정성을 핑계로 초기 차이까지 지워버리면 안 되기 때문입니다.

  

**지금 주장할 수 있는 것**

- PNU의 불리한 방향은 최소한 1년 시점부터 보입니다.
    
- 2년 observed KM 숫자는 방향성은 같지만, **PNU raw KM 자체가 불안정**합니다.
    

  

**지금 주장하면 안 되는 것**

- “2년 observed difference 17.7%p가 정확한 크기다”
    
- “이 observed 차이가 전부 순수한 질병 자연경과 차이다”
    

---

## **0-4. Model comparison block: merged standardized short-horizon contrast**

  

**Source: stage4_supported_short_horizon_contrast.csv**

|**Model**|**Horizon, y**|**Standardized PNU risk**|**Standardized SNU risk**|**PNU minus SNU**|**95% interval**|**Bootstrap success**|
|---|---|---|---|---|---|---|
|restricted_site_added|1|17.8%|7.1%|10.7%p|0.3%p to 24.5%p|100%|
|restricted_site_interaction|1|18.2%|7.1%|11.1%p|0.2%p to 26.6%p|100%|
|restricted_site_added|2|27.3%|12.8%|14.4%p|-0.6%p to 30.0%p|100%|
|restricted_site_interaction|2|28.1%|12.6%|15.5%p|essentially 0 to 31.6%p|100%|

**숫자가 직접 말하는 것**

age와 sex 차이를 맞춘 뒤에도 PNU의 short-horizon risk가 더 큽니다.

1년에서는 약 10.7–11.1%p 차이, 2년에서는 약 14.4–15.5%p 차이입니다.

  

**쉬운 말로 풀면**

“PNU가 더 불리해 보이는 게 단순히 나이나 성별 구성 차이 때문인가?”에 대해, 적어도 이 표는 **그것만으로는 설명이 안 된다**고 말합니다.

  

**그래서 중요한 이유**

이 블록은 단순한 raw KM 비교보다 한 단계 더 강합니다.

즉, **초기 timing difference가 baseline age/sex mix를 맞춘 뒤에도 남는다**는 뜻입니다.

  

**지금 주장할 수 있는 것**

- 초기 site/context 차이는 단순 age/sex 구성만의 결과는 아닐 가능성이 큽니다.
    
- 1년 차 신호는 모델을 바꿔도 방향이 일관됩니다.
    

  

**지금 주장하면 안 되는 것**

- “site effect가 곧 치료 효과다”
    
    metadata는 이를 명시적으로 금지합니다. Stage 4의 site는 broader context proxy입니다. 
    
- “2년 standardized difference가 확정적이다”
    
    2년 값은 0을 거의 스치는 수준이라 아직 경계가 필요합니다.
    

---

## **0-5. Parameter / coefficient block: piecewise site-effect**

  

**Source: stage4_piecewise_interval_support.csv and stage4_piecewise_site_effect.csv**

|**Interval**|**PNU subjects in interval**|**SNU subjects in interval**|**PNU events**|**SNU events**|**Unadjusted HR, SNU vs PNU**|**Age+sex adjusted HR**|**Age+sex+interaction HR**|**Estimable**|
|---|---|---|---|---|---|---|---|---|
|0–1y|54|208|8|14|0.39 (0.16 to 0.93)|0.44 (0.18 to 1.08)|0.42 (0.17 to 1.02)|yes|
|1–2y|28|163|4|8|0.27 (0.08 to 0.87)|0.29 (0.09 to 0.93)|0.28 (0.09 to 0.90)|yes|
|>2y|1|134|0|16|suppressed|suppressed|suppressed|no|

**숫자가 직접 말하는 것**

0–1년과 1–2년 모두에서 SNU의 hazard가 PNU보다 낮습니다.

특히 1–2년에서는 SNU vs PNU HR이 약 0.27–0.29 수준입니다.

하지만 2년 이후는 PNU가 사실상 비어 있어 추정 자체를 억제했습니다.

  

**쉬운 말로 풀면**

초기 2년 동안은 PNU에서 사건이 더 빨리, 더 자주 나옵니다.

그런데 2년이 지나면 “누가 더 위험한지”를 비교할 PNU 데이터가 거의 없어집니다.

  

**그래서 중요한 이유**

이 표는 **site effect가 시간에 따라 바뀐다**는 해석과 잘 맞습니다.

즉, 차이가 있다면 **초기**에 있고, **후반**에는 데이터가 부족해서 더 말할 수 없습니다.

후반부에 말을 멈춘 것은 약점이 아니라, 오히려 **좋은 QC**입니다.

  

**지금 주장할 수 있는 것**

- early site/context disadvantage는 piecewise model에서도 반복됩니다.
    
- late site effect는 “없다”가 아니라 “추정할 수 없다”가 더 정확합니다.
    

  

**지금 주장하면 안 되는 것**

- “2년 이후에는 두 코호트 차이가 사라진다”
    
- “2년 이후에도 SNU가 PNU보다 안전하다”
    

---

## **0-6. Hazard-pattern summary**

  

**Source: stage4_hazard_pattern_summary.csv**

|**Band**|**PNU crude hazard per 100 PY**|**SNU crude hazard per 100 PY**|**merged crude hazard per 100 PY**|**핵심 해석**|
|---|---|---|---|---|
|0–1y|21.5|7.8|10.2|PNU early event pressure가 훨씬 큼|
|1–2y|23.6|5.5|7.4|격차가 계속 유지됨|
|2–5y|0.0, but only 22.3% of band observed|3.8|3.8|PNU는 band 자체가 거의 관찰되지 않음|
|5–10y|no observed time|4.3|4.3|merged late band는 사실상 SNU-only|

**숫자가 직접 말하는 것**

PNU의 초반 hazard는 SNU보다 약 3배 안팎으로 높습니다.

반대로 2년 이후 PNU는 관찰 person-time 자체가 거의 없습니다.

  

**쉬운 말로 풀면**

PNU는 “늦게 조금씩 나빠지는 코호트”가 아니라, **처음부터 빨리 사건이 나오는 코호트**처럼 보입니다.

  

**그래서 중요한 이유**

이건 질병 자연경과가 다르거나, entry 시점의 중증도/선별구조가 다르거나, site별 관리/판정 구조가 다를 가능성을 시사합니다.

핵심은 차이가 **초기 전면부(front-loaded)** 라는 점입니다.

  

**지금 주장할 수 있는 것**

- Stage 4가 보여주는 주된 현상은 tail divergence보다 **early timing difference**입니다.
    
- late band는 pooled biology보다 **support architecture**의 영향을 더 많이 받습니다.
    

  

**지금 주장하면 안 되는 것**

- “후반 hazard는 정말 동일하다”
    
- “PNU는 장기적으로 안전하다” 혹은 반대로 “장기적으로 더 위험하다”
    

---

## **0-7. Diagnostics / convergence / assumption-check block**

  

**Source: stage4_export_manifest.csv, stage4_supported_short_horizon_contrast.csv, stage4_late_horizon_tail_contrast.csv, stage4_piecewise_site_effect.csv, stage4_plot_book.pdf**

|**Check item**|**Result**|**판정**|
|---|---|---|
|Export file structure|Manifest 기준 핵심 산출물 11개 구조 확인|통과|
|Row vs subject mismatch|세 데이터셋 모두 n_rows = n_unique_person_id|통과|
|merged arithmetic consistency|merged counts가 PNU + SNU 합과 일치|통과|
|Duplicate key structure|short-horizon, late-tail block에서 중복 key 미확인|통과|
|Piecewise count alignment|support table와 site-effect table count 일치|통과|
|Reporting flag coherence|suppressed / unstable / reported 상태 모순 없음|통과|
|Short-horizon instability|PNU 2y observed risk, PNU-SNU 2y observed difference만 bootstrap unstable|경고|
|Late-tail estimability|PNU 5y, 10y observed risk suppressed; >2y piecewise HR 모두 suppressed|경고|
|Optimizer code/message audit|CSV/PDF에 상세 convergence code 미노출|경미한 한계|

plot book도 이 패턴과 맞습니다.

1쪽은 PNU가 2년 전후까지 가파르게 상승하고 멈추는 반면 SNU와 merged는 더 길게 이어짐을 보여주고, 3쪽은 piecewise HR이 0–1, 1–2년만 표시되며, 4쪽은 >2년 구간에서 PNU event support가 사실상 없음을 시각적으로 보여줍니다. 6쪽에서는 late-tail plot에 PNU가 빠지고 SNU와 merged만 그려져 있어 PNU의 장기 tail 미지원이 그림에서도 일관됩니다. 

  

**숫자가 직접 말하는 것**

문제는 “계산이 틀렸다”가 아니라 “후반부에 말할 데이터가 없다”입니다.

  

**쉬운 말로 풀면**

컴퓨터가 헛돌아서 이상한 숫자를 만든 흔적은 없고, 오히려 데이터가 없는 곳에서는 숫자를 숨기거나 불안정 플래그를 달아둔 구조입니다.

  

**그래서 중요한 이유**

이건 연구 신뢰도에 유리합니다.

late tail을 억지로 예쁘게 그리지 않고, unsupported 구간을 억제했기 때문입니다.

  

**지금 주장할 수 있는 것**

- no major computational failure.
    
- interpretive limitation이 주 문제이지, coding breakdown이 주 문제는 아닙니다.
    

  

**지금 주장하면 안 되는 것**

- “late-tail suppression이 모델 실패다”
    
- “PNU long-tail risk는 0이다”
    

---

# **A. Scientific question reframed**

  

이 Stage 4의 진짜 질문은

**“어느 모델이 더 좋아 보이느냐”**가 아니라,

  

**“PNU와 SNU의 차이가 late tail에서 생긴 착시인지, 아니면 공통 1–2년 창에서 이미 나타나는 timing difference인지”**

를 먼저 분리해 보는 것입니다.

  

즉, 이 단계의 estimand는 **cure fraction**이 아니라,

**같은 시간축 위에서 관찰되는 초기 transition timing 차이와 그 지지 정도**입니다.

master specification도 timing difference, follow-up immaturity, cure evidence를 절대 섞지 말라고 요구합니다. 

---

# **B. Bottom-line integrated interpretation**

- **PNU는 SNU보다 더 이른 시기에 transition이 집중되는 패턴**을 보입니다. 이건 raw KM, merged standardized risk difference, piecewise HR가 모두 같은 방향을 보여줘서 단순 우연 한 줄로 치부하기 어렵습니다.
    
- 그 차이는 **late tail이 아니라 이미 1년 공통 창에서 시작**됩니다. 따라서 “후반부 곡선이 벌어져 보여서 생긴 착시”로만 설명하기 어렵습니다.
    
- 하지만 **2년 이후 비교는 급격히 질이 떨어집니다.** PNU는 2년을 넘는 순간 정보가 거의 사라지고, merged long tail은 사실상 SNU가 대표합니다.
    
- 따라서 이 결과가 가장 강하게 말하는 것은 **cure**가 아니라 **timing difference 또는 site-context difference**입니다.
    
- 추가로, **PNU에만 remission이 존재하고 그것을 censoring으로 처리했기 때문에**, 현재의 PNU excess risk 일부는 질병 자연경과 그 자체보다 **결과 정의 구조의 비대칭성**을 반영했을 가능성이 있습니다. Stage 9가 꼭 필요합니다. 
    

---

# **C. Statistical interpretation**

  

첫째, **1년 차 차이**는 가장 깨끗합니다.

Observed 1-year risk는 PNU 17.0%, SNU 7.3%로 약 9.7%p 차이입니다.

이 값만 보면 신뢰구간이 0을 살짝 가로질러 아주 단정적이지는 않지만, **standardized 1-year contrast**는 10.7–11.1%p로 양수이며 bootstrap 성공률도 100%입니다. 즉, **초기 차이의 방향성과 대략적 크기**는 age/sex adjustment 후에도 유지됩니다.

  

둘째, **2년 차 차이**도 방향은 동일합니다.

Observed difference는 17.7%p로 더 커 보이지만, PNU 2-year KM는 risk-set이 1명뿐이라 bootstrap unstable입니다. 반면 restricted Cox standardization에서는 14.4–15.5%p 정도의 차이가 계속 남습니다. 즉, **2년 시점 차이는 plausible하지만 raw KM 숫자 그대로 고정해서 믿으면 안 됩니다.**

  

셋째, **piecewise site-effect 결과**는 차이가 front-loaded임을 보여줍니다.

SNU vs PNU HR이 0–1년에서 약 0.39–0.44, 1–2년에서 약 0.27–0.29입니다. 이는 SNU의 short-term transition hazard가 PNU보다 대략 56%–73% 낮다는 뜻입니다. 반면 >2년에서는 추정 자체가 억제됩니다.

이 억제는 실패가 아니라, **support가 없어 late site effect를 주장할 수 없다는 정직한 결과**입니다.

  

넷째, **hazard-pattern summary**는 이 모든 것을 다른 각도에서 확인해 줍니다.

PNU의 crude hazard는 0–1년 21.5, 1–2년 23.6 per 100 PY로 높고, SNU는 각각 7.8과 5.5입니다. 즉, **PNU는 초기에 사건이 몰리는 코호트**입니다. PNU의 장기 plateau는 “안정된 좋은 tail”이 아니라, “사람이 남지 않아 더 볼 수 없는 tail”일 가능성이 훨씬 큽니다.

---

# **D. Clinical interpretation**

  

이 파일들만 놓고 보면 질병명과 치료 세부 내용은 없으므로, 임상적 해석은 조건부입니다.

  

그 조건 아래에서 가장 defensible한 해석은 이렇습니다.

  

**1) 질병 자연경과 측면**

이 코호트는 두 가지 서로 다른 패턴을 보여줍니다.

PNU는 초기에 빠르게 사건이 발생하고, SNU는 더 완만하게 사건이 누적됩니다.

이는 “같은 질병의 같은 단계가 단순히 늦게 더 많이 생긴다”기보다, **entry 시점의 위험도 분포나 초기 취약성 구조가 다를 가능성**과 잘 맞습니다.

  

**2) 치료 맥락 측면**

metadata가 명시하듯이 site는 치료효과의 clean surrogate가 아닙니다. 오히려 **broader treatment context, care pathway, selection, follow-up structure의 proxy**입니다. 

따라서 지금 결과를 “PNU에서 치료가 나빴다”로 읽으면 안 됩니다. 더 그럴듯한 해석은,

- 더 중한 환자가 한쪽에 몰렸거나,
    
- 치료 또는 case management pathway가 달랐거나,
    
- transition 판정 문턱이나 추적 밀도가 달랐거나,
    
- referral structure가 달랐을 수 있다는 것입니다.
    

  

**3) 추적 구조 측면**

SNU의 장기 추적은 분명 강점입니다. 그러나 10년에서는 1명만 남아 있어, 후반 risk는 “숫자가 있다”는 것과 “임상적으로 견고하다”는 것이 다릅니다.

PNU는 2년 이후 거의 비어 있으므로, Stage 4의 pooled late tail은 **실제론 SNU-dominant tail**입니다.

  

**4) 결과 정의 측면**

PNU에는 remission이 13건 있고 SNU에는 0건입니다. Stage 4는 remission을 censoring으로 처리합니다. 이런 구조에서는 PNU의 transition risk가 현재 방식에서 **과대추정될 가능성**이 있습니다.

왜냐하면 remission이 실제로 transition을 막는 competing path라면, censoring으로 둘 경우 “그들이 계속 transition 위험군에 남아 있었을 것”처럼 취급하기 때문입니다.

따라서 현재 PNU excess risk의 일부는 **질병 악화 자체**가 아니라 **결과 정의 방식**의 산물일 수 있습니다.

---

# **E. Model-assumption-to-clinical-reality mapping**

  

## **Separate-cohort observed KM**

  

핵심 가정은 최소입니다.

관찰된 event와 censoring만으로 누적위험을 보여줍니다.

- **임상적으로 그럴듯한 점**: 모델 가정이 가장 약해서 “있는 그대로의 데이터”를 보여줍니다.
    
- **임상적으로 취약한 점**: cohort mix, follow-up length, remission censoring의 영향을 그대로 받습니다.
    

  

이 접근은 **초기 차이의 존재**를 보여주는 데는 좋지만,

그 차이가 왜 생겼는지를 설명하지는 못합니다.

  

## **Merged restricted Cox standardization**

  

핵심 가정은 “1년” 또는 “2년” 제한 창 안에서는 site effect를 age, sex와 함께 Cox 구조로 요약할 수 있다는 것입니다.

- **임상적으로 더 그럴듯한 점**: long tail 전체에 비례위험을 강요하지 않고, 공통 짧은 창에서만 비교합니다.
    
- **약한 점**: 여전히 site를 proxy로 쓰므로 치료, selection, ascertainment가 섞여 있습니다.
    

  

하지만 이 Stage 4 목적에는 꽤 잘 맞습니다.

왜냐하면 지금 궁금한 것은 “PNU와 SNU의 차이가 초기에 이미 있나?”이지, “10년짜리 site HR이 정확히 얼마인가?”가 아니기 때문입니다.

  

## **Merged piecewise site-effect Cox**

  

핵심 가정은 site effect가 시간 구간마다 달라질 수 있다는 것입니다.

- **임상적으로 매우 그럴듯한 점**: 현재 결과가 보여주는 패턴이 정확히 이 형태입니다. 차이는 초기에 강하고, 뒤로 갈수록 support가 사라집니다.
    
- **약한 점**: 후반부는 데이터가 거의 없어 모형이 맞는지 틀린지가 아니라 **추정 불능**이 됩니다.
    

  

즉, 이 모델은 “late effect가 없다”를 보여주는 모델이 아니라,

**“late effect를 말할 근거가 없다”**를 보여주는 모델입니다.

---

# **F. Why the observed pattern matters**

  

이 결과에서 가장 publishable한 포인트는 **모형 우열이 아니라 현상의 위치**입니다.

  

즉,

  

> **PNU-SNU 차이는 late tail에서 처음 나타나는 현상이 아니라, 이미 1년 공통 창에서 시작되는 초기 timing 차이이다.**

  

이건 중요한 메시지입니다.

왜냐하면 많은 survival 비교가 tail 모양에 끌려가서 “장기 생존 이질성” 이야기를 먼저 하게 되는데, 이번 Stage 4는 오히려 그 반대를 보여주기 때문입니다.

  

동시에 또 하나 중요한 점은,

  

> **2년 이후 pooled long-term picture는 사실상 SNU가 대표한다.**

  

따라서 late pooled curves를 보고 두 코호트의 장기 biology를 함께 논하면 안 됩니다.

이 억제된 late comparison 자체가 결과의 핵심 일부입니다.

---

# **G. Alternative explanations and threats to interpretation**

  

## **Strongly supported**

  

**1) Follow-up asymmetry**

PNU는 2년 이후 거의 정보가 없고, merged late tail은 거의 전적으로 SNU가 담당합니다.

이건 데이터 구조 자체가 직접 보여줍니다.

  

**2) Outcome-definition asymmetry**

PNU에는 remission 13건, SNU에는 0건입니다.

그리고 Stage 4 main analysis는 remission을 censoring으로 처리합니다.

이 차이는 해석에 직접 영향을 줄 수 있습니다.

  

## **Plausible**

  

**3) Baseline severity / referral / care-pathway difference**

초기 1–2년 차이가 age/sex를 맞춘 뒤에도 남는 점을 보면, site별 patient mix 또는 care context 차이가 있었을 가능성이 큽니다.

  

**4) Ascertainment intensity difference**

한 사이트가 더 촘촘히 추적하거나 더 민감하게 transition을 판정했을 수 있습니다.

  

## **Speculative**

  

**5) True cure-like heterogeneity difference**

Stage 4만으로는 이 결론을 내릴 수 없습니다.

master specification도 timing difference, follow-up immaturity, cure evidence를 분리하라고 명시합니다. 

  

### **무엇이 이 설명들을 구분해 줄까**

- remission을 competing event로 둔 Stage 9
    
- baseline severity, symptom burden, medication exposure, referral source 같은 site-level mix 변수
    
- transition adjudication rule과 visit schedule 정보
    
- Stage 6 screening과 Stage 7/8 cure model 결과
    

---

# **H. Manuscript-ready discussion paragraph**

  

본 timing-difference 분석은 PNU와 SNU 사이의 차이가 장기 추적의 말단부에서만 나타나는 현상이 아니라, 이미 공통 1년 창에서 관찰되는 초기 transition timing 차이라는 점을 보여주었다. Raw Kaplan-Meier 위험, age/sex를 보정한 restricted standardized risk, 그리고 piecewise site-effect Cox 모델은 모두 PNU에서 더 높은 초기 transition burden을 시사하였으며, 이는 단순한 demographic mix 차이만으로는 충분히 설명되지 않았다. 반면 2년 이후에는 PNU의 risk-set이 급격히 소실되어 pooled late-tail pattern이 사실상 SNU에 의해 결정되었고, 따라서 후반부 곡선 차이를 두 코호트의 장기 자연경과 차이로 직접 해석하는 것은 정당화되기 어렵다. 또한 PNU에서만 remission이 관찰되었고 이를 본 단계에서 censoring으로 처리하였으므로, 관찰된 site difference의 일부는 결과 정의 및 competing process의 비대칭성을 반영할 가능성이 있다. 따라서 본 단계의 가장 타당한 해석은 cure-like heterogeneity의 직접 증거라기보다, follow-up structure와 outcome handling의 영향을 동반한 초기 timing/context difference의 확인이라고 보는 것이다.

---

# **I. Next-step analyses**

1. **Stage 9 remission sensitivity를 최우선으로 수행**
    
    특히 PNU의 remission 13건 중 상당수가 초기 2년 안에 몰려 있어, 현재 transition excess가 competing-risk reanalysis에서 얼마나 줄어드는지 확인해야 합니다.
    
2. **Stage 6 screening과 Stage 7/8 cure block을 연결**
    
    현재 Stage 4는 cure를 말하는 단계가 아닙니다. follow-up sufficiency와 cure appropriateness screening을 확인한 뒤에야 cure 해석이 정당해집니다. follow-up이 불충분한 구간에서 cure/no-cure 차이가 커 보이는 것은 쉽게 왜곡될 수 있습니다.   
    
3. **공통 1년, 2년 restricted models에 더 많은 baseline covariate 추가**
    
    현재는 age, sex 중심입니다. severity, medication, symptom profile, referral source가 있으면 초기 site difference의 성격을 훨씬 더 잘 설명할 수 있습니다.
    
4. **Transition ascertainment protocol audit**
    
    각 site의 visit spacing, 판정 기준, 누가 event를 확정했는지 확인해야 합니다. 초기 차이가 biology가 아니라 detection intensity일 수도 있기 때문입니다.
    
5. **merged late tail은 SNU-dominant로 별도 라벨링**
    
    이후 통합 비교표에서는 merged 5년 이후 결과에 “effectively SNU tail” 경고를 명시하는 것이 좋습니다.
    

---

# **J. 핵심 숫자의 쉬운 해석**

  

## **1) PNU 1년 observed risk 17.0% vs SNU 7.3%**

1. **숫자가 문자 그대로 뜻하는 것**
    
    1년까지의 Kaplan-Meier 기준 누적 transition risk가 PNU 17.0%, SNU 7.3%라는 뜻입니다.
    
2. **쉬운 말로 바꾸면 무엇인지**
    
    100명 중 대략 PNU는 17명, SNU는 7명 정도가 1년 안에 transition에 도달한 셈입니다.
    
3. **관찰 데이터가 이 숫자를 얼마나 직접 뒷받침하는지**
    
    꽤 직접적입니다. PNU도 1년 시점에는 28명이 horizon 너머로 남아 있고, SNU는 163명이 남아 있어 둘 다 short horizon으로는 해석이 가능합니다.
    
4. **그래서 이 숫자를 어느 정도까지 믿어도 되는지**
    
    방향성은 비교적 믿을 만합니다. 다만 PNU 표본 자체가 작아 신뢰구간은 넓습니다.
    
5. **이 숫자가 임상적/과학적으로 요구하는 해석이 무엇인지**
    
    코호트 차이는 이미 초기 1년 안에 존재할 가능성이 높고, 이는 late-tail artifact만으로 설명되기 어렵습니다.
    

  

## **2) 1년 standardized PNU minus SNU difference 10.7%p to 11.1%p**

1. **숫자가 문자 그대로 뜻하는 것**
    
    age와 sex 분포를 맞춘 뒤에도 PNU의 1년 risk가 SNU보다 약 10.7–11.1%p 높다는 뜻입니다.
    
2. **쉬운 말로 바꾸면 무엇인지**
    
    “PNU가 더 나빠 보이는 이유가 나이와 성별 때문만은 아니다”는 뜻입니다.
    
3. **관찰 데이터가 이 숫자를 얼마나 직접 뒷받침하는지**
    
    1년 창은 primary-supported이고, bootstrap 성공률도 100%라 상대적으로 안정적입니다.
    
4. **그래서 이 숫자를 어느 정도까지 믿어도 되는지**
    
    크기의 exact decimal보다는, **양의 차이가 반복해서 나온다**는 사실을 믿는 편이 맞습니다.
    
5. **이 숫자가 임상적/과학적으로 요구하는 해석이 무엇인지**
    
    site/context difference는 적어도 일부는 composition beyond age/sex를 반영합니다.
    

  

## **3) PNU 2년 observed risk 29.9%, but n_risk = 1**

1. **숫자가 문자 그대로 뜻하는 것**
    
    2년 KM risk는 29.9%로 계산되지만, 2년 시점 너머에 남아 있는 사람이 1명뿐입니다.
    
2. **쉬운 말로 바꾸면 무엇인지**
    
    숫자는 찍히지만, 거의 마지막 한두 명의 관찰에 기대고 있다는 뜻입니다.
    
3. **관찰 데이터가 이 숫자를 얼마나 직접 뒷받침하는지**
    
    매우 약합니다. 실제로 bootstrap instability flag가 붙어 있습니다.
    
4. **그래서 이 숫자를 어느 정도까지 믿어도 되는지**
    
    방향성 참고용 정도로만 써야 하고, 크기 자체는 과신하면 안 됩니다.
    
5. **이 숫자가 임상적/과학적으로 요구하는 해석이 무엇인지**
    
    2년 PNU raw KM은 “설명용 숫자”이지 “정밀 추정치”가 아닙니다.
    

  

## **4) SNU 5년 observed risk 21.1%, n_risk = 64**

1. **숫자가 문자 그대로 뜻하는 것**
    
    SNU에서는 5년 risk가 21.1%입니다.
    
2. **쉬운 말로 바꾸면 무엇인지**
    
    100명 중 약 21명 정도가 5년 안에 transition을 경험한 수준입니다.
    
3. **관찰 데이터가 이 숫자를 얼마나 직접 뒷받침하는지**
    
    중간 정도입니다. 5년에도 64명이 남아 있어 완전한 tail projection만은 아닙니다.
    
4. **그래서 이 숫자를 어느 정도까지 믿어도 되는지**
    
    short horizon만큼 튼튼하진 않지만, SNU 단독으로는 아직 의미 있는 관찰 근거가 있습니다.
    
5. **이 숫자가 임상적/과학적으로 요구하는 해석이 무엇인지**
    
    SNU는 장기적으로 event-free subgroup이 상당수 존재하는 것처럼 보일 수 있으나, 이것을 바로 cure로 읽어선 안 됩니다.
    

  

## **5) SNU 10년 risk 32.6%, merged 10년 risk 35.1%, 둘 다 n_risk = 1**

1. **숫자가 문자 그대로 뜻하는 것**
    
    10년 risk 숫자는 존재하지만 horizon 너머에 남은 사람이 1명뿐입니다.
    
2. **쉬운 말로 바꾸면 무엇인지**
    
    10년 숫자는 사실상 마지막 한 사람의 꼬리 끝점입니다.
    
3. **관찰 데이터가 이 숫자를 얼마나 직접 뒷받침하는지**
    
    매우 약합니다. 형식상 observed이지만 해석상 projection에 가깝습니다.
    
4. **그래서 이 숫자를 어느 정도까지 믿어도 되는지**
    
    “대충 이 정도까지 올라갈 수 있다”는 수준 이상으로 쓰면 안 됩니다.
    
5. **이 숫자가 임상적/과학적으로 요구하는 해석이 무엇인지**
    
    10년 risk로 biological tail difference나 cure fraction을 논하는 것은 시기상조입니다.
    

  

## **6) merged >2 interval support: PNU 1 subject, 0 events vs SNU 134 subjects, 16 events**

1. **숫자가 문자 그대로 뜻하는 것**
    
    merged의 >2년 구간은 PNU가 사실상 비어 있고 SNU가 거의 전부를 차지합니다.
    
2. **쉬운 말로 바꾸면 무엇인지**
    
    2년 이후 merged는 겉으로만 합쳐진 데이터이고, 실제 내용물은 거의 SNU입니다.
    
3. **관찰 데이터가 이 숫자를 얼마나 직접 뒷받침하는지**
    
    매우 직접적입니다. interval support table 자체가 그렇게 말합니다.
    
4. **그래서 이 숫자를 어느 정도까지 믿어도 되는지**
    
    이건 거의 확정적인 구조적 사실입니다.
    
5. **이 숫자가 임상적/과학적으로 요구하는 해석이 무엇인지**
    
    merged late-tail 결과는 pooled cohort biology가 아니라 **SNU-dominant follow-up architecture**로 읽어야 합니다.
    

---

# **K. 이 결과가 독자에게 실제로 요구하는 해석과 주의점**

  

독자는 이 결과를 다음처럼 사용해야 합니다.

- **강조할 것**: 1년 중심, 그다음 2년까지의 early timing difference
    
- **보류할 것**: 5–10년 decimal-level 차이
    
- **분리할 것**: timing difference, follow-up immaturity, remission handling, cure evidence
    
- **경계할 것**: site effect를 treatment effect로 읽는 것
    
- **라벨링할 것**: merged 5년 이후는 사실상 SNU tail
    
- **우선순위 둘 것**: remission competing-risk sensitivity와 cure appropriateness screening
    

  

즉, 지금 Stage 4가 독자에게 요구하는 태도는

**“초기 차이는 진지하게 보되, 후반 꼬리는 겸손하게 읽어라”** 입니다.

---

# **L. 한 줄 요약**


**이 Stage 4 결과는 PNU가 SNU보다 초기에 더 빨리 transition에 도달하는 패턴을 보여주지만, 2년 이후 비교는 거의 SNU tail과 PNU의 remission 처리 방식에 기대므로, 지금 내릴 결론은 ‘초기 timing 차이’이지 ‘cure 존재’가 아닙니다.**



