
# 0. Evidence tables

먼저 결론부터 말하면, **업로드된 Stage 5 결과물에서는 중대한 계산 실패는 확인되지 않았습니다.** 다만, **해석상 주의가 매우 큰 구간**이 분명합니다. 특히 **PNU의 3년 이후**, 그리고 **SNU·merged의 8–10년**은 관찰자료가 아니라 **모형 꼬리 구간 해석**의 비중이 큽니다. 또한 이 프로젝트의 공통 규칙상 **transition만 사건**, **remission은 main 분석에서 censoring**으로 처리되므로, Stage 5만으로는 “진짜 cure”를 주장할 수 없습니다. 이 프로젝트 문서도 **timing difference**, **follow-up immaturity**, **cure-like heterogeneity**를 절대 섞어 해석하지 말라고 명시합니다. 또 `site`는 치료 그 자체가 아니라 **care pathway, selection, follow-up structure의 대리변수**로 보라고 되어 있습니다.

---

## 0-1. Data structure / cohort summary

**Source:** `stage5_export_manifest.csv`, `stage5_model_registry.csv`, `stage5_subject_horizon_risk_long.csv`, `stage5_qc_summary.csv`

|항목|PNU|SNU|merged|
|---|--:|--:|--:|
|분석 대상자 수 `n_obs`|54|208|262|
|transition 사건 수 `n_event`|12|38|50|
|공식 수|2|2|4|
|모델 family 수|5|5|5|
|model registry 행 수|10|10|20|
|subject-horizon risk 행 수|5400|20800|52400|
|마지막 관찰 사건 시점, 년|1.34|7.57|7.57|
|마지막 추적 종료 시점, 년|2.67|10.11|10.11|

추가 QC 요약:

|QC 항목|값|
|---|--:|
|planned model count|40|
|fitted registry rows|40|
|converged model count|40|
|expected prediction rows|78600|
|actual prediction rows|78600|
|risk range violation|0|
|survival range violation|0|
|monotone risk violation|0|
|unsupported metric population violation|0|
|unstable calibration estimate nonmissing|0|
|calibration boundary rows|11|

### 숫자가 직접 말하는 것

- 40개 모델이 모두 적합되었고, 예측 행 수가 기대값과 정확히 일치했습니다.
    
- 위험도 범위 이탈, 생존확률 범위 이탈, 단조성 위반이 모두 0입니다.
    
- 즉, **Stage 5 코드가 의도한 구조적 QC는 통과**했습니다.
    

### 쉬운 말로 풀면

- 계산이 깨졌거나, 위험도가 1보다 크거나 음수가 나온다거나, 행 수가 꼬였다거나 하는 **큰 사고는 없었다**는 뜻입니다.
    
- 다만 `risk_long` 행 수가 큰 것은 중복 오류가 아니라, **사람 × horizon × 공식 × family로 길게 쌓인 구조**라서 그렇습니다.
    

### 그래서 중요한 이유

- 이 정도 QC가 통과했다면, 뒤의 해석은 “계산이 잘못돼서 생긴 착시”보다는 **자료 구조와 추적 구조, 모형 가정**의 문제로 봐야 합니다.
    

### 지금 주장할 수 있는 것

- **중대한 계산적 실패는 없다.**
    
- `merged`는 `PNU + SNU`를 다시 쌓은 export layer라서, **세 dataset의 사람 수를 합산하면 안 된다.**
    

### 지금 주장하면 안 되는 것

- “QC가 통과했으니 과학적 결론도 튼튼하다.”
    
- 그건 아닙니다. QC 통과는 **계산이 멀쩡하다**는 뜻이지, **late tail이 믿을 만하다**는 뜻은 아닙니다.
    

---

## 0-2. Event and follow-up summary

**Source:** `stage5_model_registry.csv`, `stage5_model_performance_long.csv`

|dataset|last event year|last follow-up year|1년 observed risk|2년 observed risk|5년 observed risk|10년 observed risk|
|---|--:|--:|--:|--:|--:|--:|
|PNU|1.34|2.67|0.170|0.299|0.299|0.299|
|SNU|7.57|10.11|0.073|0.122|0.211|0.326|
|merged|7.57|10.11|0.092|0.153|0.239|0.351|

자연빈도로 바꾸면:

- **PNU**: 1년 약 **100명 중 17명**, 2년 약 **100명 중 30명**
    
- **SNU**: 1년 약 **100명 중 7명**, 2년 약 **100명 중 12명**, 5년 약 **100명 중 21명**
    
- **merged**: 1년 약 **100명 중 9명**, 2년 약 **100명 중 15명**, 5년 약 **100명 중 24명**
    

단, PNU의 5년·10년 수치는 아래 support table에서 보듯 **관찰로 직접 뒷받침되는 값이 아닙니다.**

### 숫자가 직접 말하는 것

- PNU는 사건이 **1.34년 이후 더 이상 관찰되지 않았고**, 전체 추적도 **2.67년**이 최대입니다.
    
- 반면 SNU와 merged는 사건이 **7.57년까지 관찰**됩니다.
    

### 쉬운 말로 풀면

- PNU는 “빨리 일어날 사람은 빨리 일어나고, 그 뒤는 많이 비어 있는” 구조입니다.
    
- SNU와 merged는 그보다 훨씬 늦은 시점까지 실제 사건이 있어서, **중장기 위험을 조금 더 직접 볼 수 있습니다.**
    

### 그래서 중요한 이유

- 똑같은 10년 위험이라도, PNU의 10년 위험은 거의 **상상해서 그린 꼬리**이고, SNU·merged의 10년 위험은 **일부는 관찰, 일부는 꼬리**입니다.
    
- 숫자가 같아 보여도 믿을 수 있는 정도가 다릅니다.
    

### 지금 주장할 수 있는 것

- **조기 위험은 PNU가 더 가파르다.**
    
- **중장기 누적 위험은 SNU와 merged에서 더 직접적으로 평가 가능하다.**
    

### 지금 주장하면 안 되는 것

- PNU의 5년, 10년 plateau를 보고 “질병이 멈췄다” 또는 “cure가 있다.”
    
- 그건 Stage 5만으로는 말할 수 없습니다.
    

---

## 0-3. Follow-up support by horizon

**Source:** `stage5_model_performance_long.csv`, `stage5_model_registry.csv`  
참고로 이 프로젝트는 horizon 지원 강도를 PNU는 1년 primary, 2년 sensitivity, 3년 이상 projection, SNU·merged는 1–2년 primary, 3–5년 secondary, 6년 이상 projection으로 해석하라고 명시합니다.

|dataset|horizon|known_count|binary outcome support|data status|해석 consequence|
|---|--:|--:|--:|---|---|
|PNU|1|36|있음|primary-supported|직접 관찰에 가장 가까움|
|PNU|2|13|있음|sensitivity|아주 조심해서만 사용|
|PNU|5|12|없음|projection beyond last follow-up|사실상 모형 꼬리|
|PNU|10|12|없음|projection beyond last follow-up|사실상 모형 꼬리|
|SNU|1|177|있음|primary-supported|안정적|
|SNU|2|156|있음|primary-supported|안정적|
|SNU|5|98|있음|secondary|여전히 꽤 쓸 만함|
|SNU|10|39|있음|projection|late tail 주의|
|merged|1|213|있음|primary-supported|안정적|
|merged|2|169|있음|primary-supported|안정적|
|merged|5|110|있음|secondary|중기 해석 가능|
|merged|10|51|있음|projection|late tail 주의|

보강 설명:

- **PNU**는 2년 이후 `binary_outcome_support_flag = 0`
    
- **SNU·merged**는 10년까지 형식상 support는 남아 있지만, **마지막 사건이 7.57년**이라 8–10년은 실제로는 **사건이 더 안 쌓인 상태에서 꼬리 모양에 의존**합니다.
    

### 숫자가 직접 말하는 것

- PNU는 2년 이후에 fixed-horizon outcome을 제대로 만들 수 있는 nonevent support가 사라집니다.
    
- SNU·merged는 5년까지는 꽤 버티지만, 10년에서는 known_count가 39, 51까지 줄어듭니다.
    

### 쉬운 말로 풀면

- PNU의 3년 이후 숫자는 “본 숫자”가 아니라 “그려진 숫자”에 가깝습니다.
    
- SNU와 merged의 10년 숫자도 완전히 허공은 아니지만, **상당히 model-dependent**입니다.
    

### 그래서 중요한 이유

- 같은 10년 위험이라도 **PNU 10년**은 거의 전적으로 projection,
    
- **SNU·merged 10년**은 partially observed + projection입니다.
    

### 지금 주장할 수 있는 것

- **PNU는 1년이 핵심, 2년은 민감도 확인 수준, 3년 이상은 투영**입니다.
    
- **SNU·merged는 5년까지가 상대적으로 실질적 해석 구간**입니다.
    

### 지금 주장하면 안 되는 것

- late horizon의 소수점 차이를 실재 biological phenomenon처럼 해석하는 것.
    
- 특히 PNU 5년과 10년의 family 차이는 **과학적 사실**이라기보다 **분포 꼬리 선택의 산물**입니다.
    

---

## 0-4. Model comparison table(s)

### 0-4-1. PNU Base complete comparison block

**Source:** `stage5_model_performance_long.csv`

|family|1y predicted risk|2y predicted risk|5y predicted risk|10y predicted risk|1y AUC|2y AUC|1y Brier|2y Brier|
|---|--:|--:|--:|--:|--:|--:|--:|--:|
|observed reference|0.170|0.299|0.299|0.299|NA|NA|NA|NA|
|coxph|0.171|0.308|0.308|0.308|0.673|0.736|0.138|0.180|
|exponential|0.213|0.373|0.664|0.865|0.673|0.736|0.139|0.173|
|weibull|0.220|0.332|0.531|0.699|0.673|0.736|0.139|0.177|
|lognormal|0.223|0.312|0.447|0.555|0.673|0.736|0.139|0.179|
|loglogistic|0.220|0.320|0.480|0.606|0.673|0.736|0.139|0.176|

해석 포인트:

- 1–2년에서는 family 간 차이가 크지 않습니다.
    
- 5–10년에서는 **0.308 대 0.865**까지 크게 벌어집니다.
    
- 그런데 이 벌어짐은 **PNU 2년 이후가 unsupported**라는 사실 위에서 생깁니다.
    

### 숫자가 직접 말하는 것

- PNU에서 10년 예측위험은 coxph 0.308, exponential 0.865로 엄청 다릅니다.
    
- 반면 1년 AUC와 2년 AUC는 family 모두 거의 같습니다.
    

### 쉬운 말로 풀면

- **초기 구간은 데이터가 말하고**, 뒤 구간은 **모형이 대신 말합니다.**
    
- 그래서 앞부분은 비슷하고, 뒷부분은 family마다 제멋대로 벌어집니다.
    

### 그래서 중요한 이유

- PNU에서 “late tail 모양”으로 과학적 이야기를 하는 건 매우 위험합니다.
    
- 오히려 이 표는 **cure를 말해준다기보다 follow-up immaturity를 말해줍니다.**
    

### 지금 주장할 수 있는 것

- PNU에서는 **1년 위험 추정은 비교적 안정적**입니다.
    
- 2년은 이미 민감도 분석 수준이고, 그 이후는 **distribution-driven projection**입니다.
    

### 지금 주장하면 안 되는 것

- PNU에서 어떤 parametric family가 late risk를 더 높게 줬으니 그게 진실이라고 말하는 것.
    
- 특히 coxph의 평평한 tail은 “질병이 멈췄다”가 아니라, **마지막 사건 이후 baseline hazard를 더 올리지 않는 비교자 규칙**의 영향입니다.
    

---

### 0-4-2. SNU Base complete comparison block

**Source:** `stage5_model_performance_long.csv`

|family|1y predicted risk|2y predicted risk|5y predicted risk|10y predicted risk|1y AUC|2y AUC|5y AUC|10y AUC|1y Brier|2y Brier|5y Brier|10y Brier|
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
|observed reference|0.073|0.122|0.211|0.326|NA|NA|NA|NA|NA|NA|NA|NA|
|coxph|0.073|0.121|0.211|0.326|0.484|0.575|0.474|0.252|0.067|0.106|0.166|0.222|
|exponential|0.051|0.099|0.228|0.404|0.493|0.585|0.479|0.252|0.068|0.107|0.167|0.228|
|weibull|0.077|0.123|0.222|0.338|0.496|0.588|0.478|0.252|0.068|0.106|0.166|0.222|
|lognormal|0.084|0.133|0.222|0.307|0.502|0.592|0.481|0.270|0.068|0.106|0.166|0.221|
|loglogistic|0.078|0.125|0.223|0.328|0.502|0.591|0.481|0.270|0.068|0.106|0.166|0.221|

### 숫자가 직접 말하는 것

- SNU는 1–5년 예측 위험이 observed risk와 꽤 비슷하게 모입니다.
    
- 하지만 10년 AUC는 0.252–0.270까지 떨어집니다.
    

### 쉬운 말로 풀면

- SNU는 **5년까지는 family를 바꿔도 큰 줄기가 비슷**합니다.
    
- 그런데 10년쯤 가면 “누가 더 위험한 사람인지 순위를 잘 맞추는 능력”이 거의 무너집니다.
    

### 그래서 중요한 이유

- 이는 “질병 메커니즘이 10년에 갑자기 뒤집혔다”기보다, **후반부 사건이 적고 risk set이 줄어드는 late-tail 문제**로 보는 게 더 타당합니다. 도표에서도 9–10년 AUC 급락과 Brier 상승이 시각적으로 확인됩니다.
    

### 지금 주장할 수 있는 것

- SNU는 **1–5년 절대위험 추정은 비교적 안정적**입니다.
    
- 8–10년은 **숫자는 있으나 신뢰도는 눈에 띄게 떨어진다**고 말할 수 있습니다.
    

### 지금 주장하면 안 되는 것

- 10년 AUC 하락을 생물학적 역전이나 subgroup reversal로 단정하는 것.
    

---

### 0-4-3. merged formula comparison block

**Source:** `stage5_model_performance_long.csv`

|merged formula|1y predicted risk range|1y AUC range|1y Brier range|2y predicted risk range|2y AUC range|2y Brier range|5y predicted risk range|5y AUC range|5y Brier range|10y predicted risk range|10y AUC range|10y Brier range|
|---|---|---|---|---|---|---|---|---|---|---|---|---|
|Base|0.065–0.109|0.577–0.579|0.083–0.084|0.125–0.162|0.656–0.659|0.125–0.126|0.245–0.282|0.563–0.576|0.178–0.179|0.339–0.478|0.437–0.559|0.218–0.226|
|Interaction|0.065–0.108|0.570–0.577|0.083–0.084|0.125–0.161|0.650–0.654|0.125–0.126|0.245–0.281|0.566–0.581|0.177–0.178|0.338–0.477|약 0.559|0.216–0.224|
|Site-added|0.083–0.112|0.601–0.604|0.081|0.154–0.168|0.723–0.728|0.110–0.116|0.266–0.321|0.617–0.635|0.156–0.166|0.356–0.504|약 0.583|0.205–0.212|
|Site + interaction|0.082–0.112|0.592–0.602|0.081|0.153–0.168|0.706–0.721|0.110–0.115|0.266–0.320|0.621–0.660|0.154–0.165|0.355–0.503|약 0.583|0.202–0.209|

### 숫자가 직접 말하는 것

- merged에서 **site를 넣은 공식**은 1–5년 AUC가 확실히 높고, Brier는 더 낮습니다.
    
- 예를 들어 2년 AUC는 site-free Base 0.656–0.659인데, Site-added는 0.723–0.728입니다.
    

### 쉬운 말로 풀면

- **PNU와 SNU를 한데 섞어서 site를 무시하면**, 누가 더 빨리 transition할지 구분하는 능력이 떨어집니다.
    
- 반대로 **site를 넣으면**, 초기와 중기 예측이 훨씬 정돈됩니다.
    

### 그래서 중요한 이유

- 이건 “어느 분포가 이겼냐”보다 훨씬 중요한 신호입니다.
    
- **조기 위험 구조 자체가 cohort/site 맥락에 따라 다르다**는 뜻이기 때문입니다.
    
- 프로젝트 규칙상 `site`는 순수 치료효과가 아니라 **care pathway, selection, follow-up structure의 대리변수**로 해석해야 하므로, 이 결과는 “치료가 다르다”보다는 **관찰된 전이 timing이 cohort context에 따라 다르다**는 쪽이 더 방어적입니다.
    

### 지금 주장할 수 있는 것

- **site-free merged model의 교환가능성 가정이 덜 그럴듯하다.**
    
- **PNU와 SNU를 같은 초기 위험 구조로 보는 것은 무리다.**
    

### 지금 주장하면 안 되는 것

- site coefficient improvement를 곧바로 특정 약물효과나 센터 치료효과로 번역하는 것.
    

---

## 0-5. Parameter / coefficient table(s)

**Source status:** 업로드된 Stage 5 결과물에는 **계수표 자체가 포함되어 있지 않습니다.**

확인 가능한 것은:

- `stage5_model_registry.csv`의 적합 상태, 표본 수, 사건 수, logLik/AIC/BIC, 마지막 사건 시점, 마지막 추적 시점
    
- 그러나 **각 covariate의 HR, coefficient, standard error, p-value**는 Stage 5 업로드 묶음만으로는 직접 검토할 수 없습니다.
    

### 숫자가 직접 말하는 것

- Stage 5는 성능과 예측 결과 중심 export이고, coefficient export는 빠져 있습니다.
    

### 쉬운 말로 풀면

- “site가 중요해 보인다”는 건 **예측 성능**으로는 보이지만, 지금은 **계수의 방향과 크기**를 테이블로 직접 확인할 수는 없습니다.
    

### 그래서 중요한 이유

- 어떤 site effect가 “발병 시점을 앞당긴다/늦춘다” 같은 구체적 기전 해석은, 계수표 없이 너무 멀리 나가면 안 됩니다.
    

### 지금 주장할 수 있는 것

- **예측 성능상 site 정보가 중요하다.**
    

### 지금 주장하면 안 되는 것

- age, sex, interaction의 정량 effect를 coefficient 수준으로 단정하는 것.
    

---

## 0-6. Prediction or time-specific risk table(s)

### 0-6-1. Cohort-specific observed risk comparison block

**Source:** `stage5_model_performance_long.csv`

|horizon|PNU observed risk|SNU observed risk|merged observed risk|해석|
|---|--:|--:|--:|---|
|1년|0.170|0.073|0.092|PNU가 조기 위험이 가장 높음|
|2년|0.299|0.122|0.153|PNU 조기 누적 전이 가파름|
|5년|0.299|0.211|0.239|PNU는 projection, SNU·merged는 관찰 뒷받침 존재|
|10년|0.299|0.326|0.351|PNU는 projection, SNU·merged는 late tail mixed|

### 숫자가 직접 말하는 것

- 1–2년에는 PNU가 훨씬 높고,
    
- 5년 이후에는 SNU와 merged가 계속 쌓입니다.
    

### 쉬운 말로 풀면

- PNU는 “빨리 전환될 사람은 빨리 전환된다”는 패턴,
    
- SNU는 “초기엔 낮지만 중기 이후까지 사건이 조금씩 이어진다”는 패턴입니다.
    

### 그래서 중요한 이유

- 이것은 **질병 자연경과가 cohort마다 다르게 보인다**는 뜻일 수 있지만,
    
- 동시에 **선별 구조, 추적 구조, care pathway 차이**의 흔적일 수도 있습니다.
    

### 지금 주장할 수 있는 것

- PNU와 SNU는 **같은 timing process로 보기 어렵다.**
    

### 지금 주장하면 안 되는 것

- PNU는 biologically more severe, SNU는 biologically milder라고 단정하는 것.
    

---

### 0-6-2. Clinical usefulness comparison block

**Source:** `stage5_model_performance_long.csv`, 시각적 corroboration `stage5_summary_plots.pdf` page 5.

|dataset and horizon|threshold|model net benefit range|treat-all|treat-none|해석|
|---|--:|--:|--:|--:|---|
|PNU Base 2y|0.15|0.175|0.175|0|모델이 사실상 treat-all과 같음|
|PNU Base 2y|0.20|0.124|0.124|0|모델 선택 이득 거의 없음|
|SNU Base 2y|0.15|0.005–0.017|-0.034|0|moderate threshold에서 treat-all보다 낫다|
|merged Base 2y|0.15|0.028–0.042|0.004|0|모델이 약간의 실용적 이득 제공|
|merged Base 2y|0.20|0.019–0.021|-0.058|0|treat-all보다 분명히 낫다|
|merged Site-added 5y|0.20|0.050–0.065|0.050|0|site-aware 모델의 추가 이득은 존재하지만 크지는 않음|

### 숫자가 직접 말하는 것

- PNU 2년은 모델들이 거의 treat-all과 겹칩니다.
    
- SNU와 merged는 2년 moderate threshold에서 model-based decision이 treat-all보다 낫습니다.
    
- 그러나 그 차이가 아주 거대하진 않습니다.
    

### 쉬운 말로 풀면

- PNU에서는 “누가 더 위험한지 고르는 모델”보다, 그냥 대부분을 고위험으로 보는 전략과 별 차이가 없습니다.
    
- SNU와 merged에서는 그보다 조금 낫지만, **임상적으로 극적인 개선**은 아닙니다.
    

### 그래서 중요한 이유

- 이 결과는 Stage 5 no-cure block의 핵심이 “한 분포의 우승”이 아니라,
    
- **어느 cohort에서는 모델을 써도 분류 이득이 작고, 어느 cohort에서는 moderate threshold에서만 조금 도움이 된다**는 점을 보여줍니다.
    

### 지금 주장할 수 있는 것

- **모델의 임상적 분류 이득은 낮은 threshold에서는 제한적**입니다.
    
- **moderate threshold에서만 modest gain**이 보입니다.
    

### 지금 주장하면 안 되는 것

- Stage 5 no-cure model만으로 임상 의사결정이 크게 개선된다고 말하는 것.
    

---

## 0-7. Diagnostics / convergence / assumption-check table(s)

**Source:** `stage5_qc_summary.csv`, `stage5_calibration_diagnostics_long.csv`

### 0-7-1. Global diagnostics block

|diagnostic item|value|의미|
|---|--:|---|
|converged models|40 of 40|적합 자체는 성공|
|supported calibration rows|320|calibration regression 시도 가능 구간|
|successful calibration rows|309|대부분 성공|
|calibration boundary rows|11|일부 late/sparse 구간 불안정|
|calibration nonconverged rows|2|극소수 불안정|
|unstable estimate nonmissing|0|불안정 수치는 NA로 잘 억제됨|
|max Uno-Harrell abs diff|0.0209|시간범위 concordance 계산은 대체로 일관|

### 0-7-2. Calibration failure complete block

|dataset|formula|family|horizon|failure reason|
|---|---|---|--:|---|
|PNU|Base|lognormal|2|slope_boundary|
|SNU|Base|coxph|9|slope_boundary|
|SNU|Base|coxph|10|slope_boundary|
|SNU|Base|exponential|9|slope_boundary|
|SNU|Base|loglogistic|9|slope_boundary|
|SNU|Base|loglogistic|10|slope_boundary|
|SNU|Base|lognormal|9|slope_boundary|
|SNU|Base|lognormal|10|slope_boundary|
|SNU|Base|weibull|9|offset_model_not_converged; slope_boundary|
|SNU|Base|weibull|10|offset_model_not_converged; slope_boundary|
|merged|Site-added|lognormal|10|slope_boundary|

### 숫자가 직접 말하는 것

- 전반적 적합 실패는 없지만, **late horizon calibration**은 몇 군데서 경계값에 걸립니다.
    

### 쉬운 말로 풀면

- 모델이 “돌아가긴 다 돌았다.”
    
- 하지만 드문 사건과 적은 risk set이 겹치는 뒤쪽에서는 calibration slope를 안정적으로 잡기 어려웠습니다.
    

### 그래서 중요한 이유

- 이건 late tail 숫자를 덜 믿어야 한다는, 아주 직접적인 경고입니다.
    
- 다행히도 불안정한 calibration 계수는 억지로 보고되지 않고 NA로 잘 막혔습니다.
    

### 지금 주장할 수 있는 것

- **수치적 안정성은 전반적으로 양호하나, late horizon calibration은 약하다.**
    

### 지금 주장하면 안 되는 것

- 9–10년 calibration 점 몇 개를 보고 장기 예측이 잘 맞는다고 주장하는 것.
    

---

## 0-8. Cohort-specific comparison table(s)

**Source:** `stage5_model_registry.csv`, `stage5_model_performance_long.csv`, `stage5_summary_plots.pdf`

도표상으로도:

- **PNU는 AUC와 Brier가 1–2년에만 존재**
    
- **SNU·merged는 9–10년에서 AUC가 급락**
    
- **merged에서 site를 넣은 공식들이 site-free보다 전반적으로 더 높은 AUC와 더 낮은 Brier를 보임**
    
- **10년 calibration 패널은 소수 점만 남음**이 확인됩니다.
    

### 숫자가 직접 말하는 것

- 표와 그림이 같은 이야기를 합니다.
    

### 쉬운 말로 풀면

- 이건 숫자 한두 개만 튀는 게 아니라, **전체 패턴이 일관되다**는 뜻입니다.
    

### 그래서 중요한 이유

- 해석을 숫자와 그림이 서로 지지합니다.
    

### 지금 주장할 수 있는 것

- **PNU는 조기-risk cohort**, **SNU는 더 긴 tail cohort**, **merged는 site-aware 접근이 필요**하다는 패턴은 robust합니다.
    

### 지금 주장하면 안 되는 것

- 그림만 보고 정확한 수치를 새로 읽어내는 것.
    

---

# 1. QC summary

## Step 1. 파일 및 결과 QC 요약

업로드된 Stage 5 결과물 기준으로 확인된 사항은 다음과 같습니다.

### 확인된 점

- `manifest`, `registry`, `performance`, `risk_long`, `calibration`, `qc_summary`, `pdf plot` 사이 구조가 서로 맞습니다.
    
- 40개 모델 모두 converged입니다.
    
- `planned_model_count = fitted_model_registry_rows = 40`
    
- `expected_prediction_rows = actual_prediction_rows = 78600`
    
- risk/survival range, monotonicity, support gating 관련 위반 0건입니다.
    
- `subject_horizon_risk_long`의 row inflation은 **중복 오류가 아니라 model × horizon stacking**으로 설명됩니다.
    
- 단, **export layering**이 있으므로 `PNU + SNU + merged`를 하나의 독립적 대상자 집합처럼 합산하면 안 됩니다.
    

### 해석상 주의점

- **PNU 3년 이후는 구조적으로 unsupported**
    
- **SNU·merged 8–10년은 사건이 더 이상 거의 쌓이지 않는 late tail**
    
- calibration regression 일부가 boundary에 걸렸지만, 이건 fatal error가 아니라 **late sparse tail 경고**입니다.
    
- coefficient table이 없어서, 예측 성능 차이를 **계수 수준 기전**으로 직접 해석할 수는 없습니다.
    
- Stage 2 follow-up maturity source table과 Stage 4 timing-difference source table이 이번 업로드에는 없어서, follow-up completeness와 timing difference를 원본 exported table로 재확인하지는 못했습니다. 다만 프로젝트 설계상 이 둘은 Stage 5보다 앞선 해석 층위입니다.
    

---

# 2. Branching

## 2-1. Major computational failure 여부

**No major computational failure was identified.**

즉:

- 지금 문제는 “계산이 틀렸다”가 아니라,
    
- **무엇을 어디까지 말해도 되느냐**의 문제입니다.
    

## 2-2. Fatal problem vs minor caution

### Fatal problem

- 없음
    

### Minor caution

- PNU 장기 꼬리
    
- SNU·merged 8–10년 calibration과 discrimination 약화
    
- remission이 main backbone에서 censoring 처리됨
    
- site를 treatment effect처럼 읽으면 안 됨
    
- Stage 5 alone으로 cure를 말하면 안 됨
    

따라서 **코드를 전면 재작성하거나 Stage 5를 다시 돌려야 할 정도의 계산적 붕괴는 없습니다.**

---

# A. Scientific question reframed

이 결과로 다시 써야 할 과학 질문은 이겁니다.

**과학 질문**  
이 자료에서 보이는 늦은 비전환 plateau가 진짜 cure-like heterogeneity를 시사하는가, 아니면 조기 timing difference와 follow-up 구조 차이 때문에 그렇게 보이는가?

**추정하려는 임상 quantity**  
각 cohort에서 **baseline 이후 1–10년의 누적 transition risk**와, 그 위험을 바탕으로 했을 때 **몇 명을 불필요하게 고위험으로 분류하는가**입니다.

**통계 모델의 역할**  
Stage 5의 no-cure 모델들은 “cure가 맞다/아니다”를 결정하는 모델이 아니라, **everyone remains susceptible라고 가정했을 때 short-to-mid horizon prediction이 얼마나 견디는지 보여주는 comparator**입니다. 이 프로젝트도 Stage 6을 cure appropriateness screening, Stage 7–8을 cure fitting, Stage 10을 최종 통합 비교로 분리해 놓고 있습니다.

---

# B. Bottom-line integrated interpretation

1. **PNU의 늦은 평탄화는 우선 cure보다 follow-up immaturity로 읽는 것이 맞습니다.**  
    2년 이후 fixed-horizon binary outcome support가 사라지고, 마지막 사건이 1.34년에 멈춰 있기 때문에, 5–10년 위험은 관찰이 아니라 꼬리 가정의 산물입니다.
    
2. **SNU와 merged에서는 1–5년 절대위험은 비교적 안정적이지만, 8–10년은 점점 model-dependent해집니다.**  
    late tail에서 AUC가 급락하고 calibration이 boundary에 걸리는 것은, 장기 순위예측과 calibration이 약해졌다는 뜻이지 질병 생물학이 갑자기 뒤집혔다는 뜻은 아닙니다.
    
3. **merged에서 site를 포함하면 예측 성능이 눈에 띄게 좋아집니다.**  
    이건 가장 publishable한 신호입니다. 즉, cohort를 섞어놓고 하나의 공통 timing process라고 가정하는 것이 덜 그럴듯하며, 실제로는 조기 transition timing이 cohort/site 맥락에 따라 달랐다는 뜻입니다.
    
4. **임상적 분류 이득은 생각보다 크지 않습니다.**  
    아주 낮은 threshold에서는 모델이 treat-all과 거의 겹치고, 분명한 practical gain은 주로 moderate threshold에서만 작게 나타납니다.
    
5. **따라서 현재 자료가 말해주는 핵심은 “어느 분포가 이겼다”가 아니라, “조기 위험 구조와 추적 성숙도가 cohort마다 다르다”입니다.**
    

---

# C. Statistical interpretation

## 3-1. 무엇이 수치적으로 뒷받침되는가

### PNU

- 1년 observed risk 0.170, 2년 observed risk 0.299
    
- family 간 1–2년 AUC 거의 동일
    
- 5–10년 risk는 family에 따라 크게 갈라짐
    
- 이는 **앞부분은 data-driven, 뒷부분은 tail-driven**이라는 전형적 패턴입니다.
    

### SNU

- 1년 0.073, 2년 0.122, 5년 0.211까지는 안정적
    
- 10년 0.326은 숫자 자체보다, 그 숫자를 둘러싼 late-tail uncertainty를 함께 봐야 합니다.
    
- AUC가 10년에 0.25–0.27 수준으로 떨어지는 건 late horizon ranking failure를 뜻합니다.
    

### merged

- site-free보다 site-aware formulas에서 1–5년 AUC가 높고 Brier가 낮습니다.
    
- 특히 2년 AUC가 약 0.65대에서 0.72대까지 올라가는 건, **site information이 예측상 실제로 중요하다**는 강한 신호입니다.
    

## 3-2. 무엇이 수치적으로 뒷받침되지 않는가

- PNU 장기 plateau를 cure evidence로 보는 것
    
- 10년 SNU/merged의 decimal-level family 차이를 실제 차이로 보는 것
    
- no-cure model family 간 작은 차이를 biological truth ranking으로 읽는 것
    

## 3-3. p-value보다 중요한 것

여기서는 p-value보다:

- **effect size**
    
- **absolute risk**
    
- **time dependence**
    
- **horizon support**
    
- **cohort heterogeneity**  
    가 더 중요합니다.
    

특히 merged 결과는 “site를 넣느냐 마느냐”에 따라 예측력이 크게 바뀌므로, 이 자료는 애초에 **하나의 균질 코호트가 아니었다**는 메시지를 줍니다.

---

# D. Clinical interpretation

현재 임상 맥락이 완전히 제공되진 않았지만, 업로드된 자료와 프로젝트 규칙만 놓고 가장 방어적으로 해석하면 다음과 같습니다.

1. **이 질환 과정은 cohort마다 같은 속도로 진행되지 않습니다.**  
    PNU는 초기 1–2년에 전환이 더 가파르고, SNU는 더 완만하지만 오래 이어집니다.  
    쉬운 말로 하면, PNU는 “빨리 일어날 사람은 빨리 일어나는” 구조, SNU는 “처음엔 덜 보이지만 중기까지 천천히 이어지는” 구조에 가깝습니다.
    
2. **이 차이를 곧바로 생물학 차이라고 부르면 안 됩니다.**  
    프로젝트 문서상 site는 치료 그 자체가 아니라 care pathway, referral pattern, selection, follow-up maturity의 대리변수입니다. 그러므로 지금 단계의 더 안전한 말은, **관찰된 transition timing이 cohort context에 따라 달랐다**입니다.
    
3. **remission을 censoring으로 처리한 점은 특히 PNU 해석에 중요합니다.**  
    PNU는 remission이 실제로 존재하는 구조이고, SNU는 구조적으로 remission이 없습니다. 따라서 PNU에서 long non-transition이 많아 보인다고 해도, 그 일부는 “영원히 susceptible하지 않음”이라기보다 **다른 임상 경로로 빠진 것**일 수 있습니다. 데이터 사전도 remission이 PNU에만 존재한다고 명시합니다.
    
4. **따라서 Stage 5에서 보이는 late flattening은 질병의 자연경과 자체보다, 추적 구조와 결과 정의의 영향이 섞인 그림**으로 보는 것이 맞습니다.
    

---

# E. Model-assumption-to-clinical-reality mapping

## exponential / weibull / lognormal / loglogistic

핵심 가정은 **모든 사람이 끝까지 susceptible**이고, 장기 tail은 특정 parametric shape를 따른다는 것입니다.

- **무슨 뜻인가**  
    뒤로 갈수록 risk가 어떻게 늘어날지를 분포 모양이 정합니다.
    
- **왜 여기서 문제인가**  
    PNU는 그 tail을 정해줄 실제 데이터가 거의 없습니다.
    
- **무엇이 바뀌는가**  
    5–10년 risk가 family마다 크게 달라집니다.
    
- **지금 할 수 있는 주장**  
    short-to-mid horizon에서는 usable comparator
    
- **지금 하면 안 되는 주장**  
    tail shape 자체를 질병 생물학으로 읽기
    

## coxph no tail extrapolation comparator

핵심 가정은 **마지막 관찰 사건 이후 baseline hazard를 더 증가시키지 않는다**는 비교자 성격입니다.

- **무슨 뜻인가**  
    마지막 사건 이후 risk가 사실상 plateau할 수 있습니다.
    
- **왜 여기서 생겼나**  
    PNU는 마지막 사건이 1.34년, SNU·merged는 7.57년입니다.
    
- **무엇이 바뀌는가**  
    PNU coxph는 2년 이후 risk가 0.308로 그대로 멈춥니다.
    
- **지금 할 수 있는 주장**  
    observed-window comparator로는 유용
    
- **지금 하면 안 되는 주장**  
    cox plateau를 자연치유나 cure로 해석
    

## merged site-free formulas

핵심 가정은 **PNU와 SNU가 같은 기본 timing process를 공유**한다는 것입니다.

- **왜 덜 그럴듯한가**  
    site를 넣으면 AUC와 Brier가 좋아지기 때문입니다.
    
- **지금 할 수 있는 주장**  
    site-free 가정은 conceptually 더 많이 깨져 있다
    

## merged site-aware formulas

핵심 가정은 **site-level context 차이를 허용**한다는 것입니다.

- **왜 더 그럴듯한가**  
    1–5년 예측이 확실히 개선됩니다.
    
- **지금 할 수 있는 주장**  
    이 데이터에서 cohort/context heterogeneity를 무시하면 안 된다
    

---

# F. Why the observed pattern matters

이 패턴이 중요한 이유는 다음입니다.

첫째, **초기 구간에서 모델들이 비슷하다는 것 자체가 정보**입니다.  
이건 1–2년 정도는 모형이 아니라 데이터가 위험을 주도한다는 뜻입니다.

둘째, **후반 구간에서 모델들이 갈라지는 것도 정보**입니다.  
이건 특정 family가 특별해서가 아니라, **data support가 사라졌다는 경고**입니다.

셋째, **merged에서 site를 넣으면 좋아지는 것**은 가장 출판 가치가 높습니다.  
이건 “late cure-like signal이 있다/없다”보다 더 앞선 층위의 질문, 즉  
**PNU와 SNU가 애초에 같은 시간 구조를 공유하느냐**에 답하기 때문입니다.

넷째, **decision-curve 차이가 아주 크지 않다**는 점도 중요합니다.  
이건 statistical fit의 차이가 늘 clinical usefulness의 큰 차이로 이어지는 것은 아니라는 뜻입니다. survival prediction에서는 calibration과 discrimination만으로 충분하지 않고, net benefit까지 봐야 한다는 최근 권고와도 맞습니다.

---

# G. Alternative explanations and threats to interpretation

## strongly supported

1. **Follow-up maturity 차이**
    
    - 특히 PNU의 late plateau는 이 설명이 가장 강합니다.
        
2. **Between-cohort timing/context heterogeneity**
    
    - merged에서 site-aware 공식의 성능 향상이 이를 지지합니다.
        

## plausible

3. **Care pathway 또는 treatment-context 차이**
    
    - site가 proxy라는 프로젝트 규칙상 충분히 가능
        
4. **Remission handling의 왜곡**
    
    - PNU에서 remission이 censoring 처리되므로 apparent long non-transition을 만들 수 있음
        

## speculative

5. **진짜 latent cure fraction**
    
    - Stage 5만으로는 아직 speculative
        
6. **질병 생물학의 근본적 cohort 차이**
    
    - 현재 데이터만으로는 context effects와 분리 어려움
        

## 무엇이 이들을 가를까

- **Stage 2** follow-up maturity 원표
    
- **Stage 4** 1년·2년 공통창 timing difference
    
- **Stage 6** cure appropriateness screening
    
- **Stage 9** remission sensitivity
    
- 가능하면 실제 treatment/process/referral variable
    

---

# H. Manuscript-ready discussion paragraph

본 분석에서 가장 중요한 발견은 특정 no-cure 분포형이 일관되게 우월했다는 점이 아니라, 관찰 가능한 위험 구조 자체가 cohort 맥락에 따라 달랐다는 점이다. PNU에서는 transition 사건이 1.34년 이후 더 이상 관찰되지 않았고 2년 이후 fixed-horizon binary outcome support가 사실상 소실되어, 3년 이후의 위험 차이는 주로 분포 꼬리 가정에 의해 결정되었다. 반면 SNU와 merged에서는 5년까지는 관찰 자료가 절대위험 추정을 비교적 안정적으로 지지하였으나, 8–10년에서는 사건 누적이 거의 멈추고 calibration과 discrimination이 약화되어 장기 추정의 불확실성이 커졌다. 특히 merged 분석에서 site 정보를 포함한 공식이 1–5년 구간의 AUC와 Brier score를 일관되게 개선한 점은, 후기 꼬리의 cure-like 양상보다 오히려 조기 transition timing의 cohort-의존성이 더 강한 신호임을 시사한다. 따라서 현재 자료는 장기 plateau를 곧바로 cure fraction의 증거로 해석하기보다, follow-up maturity, cohort selection, care pathway, 그리고 remission handling이 복합적으로 반영된 결과로 읽는 것이 타당하며, substantive cure interpretation은 후속 cure appropriateness screening과 remission-sensitive 분석을 통해 보강될 필요가 있다.

---

# I. Next-step analyses

1. **Stage 2 원표를 붙여서 horizon support를 숫자로 직접 제시하기**  
    특히 PNU 2년 이후와 SNU·merged 8–10년의 numbers-at-risk, reverse-KM, completeness를 같이 보여야 합니다.
    
2. **Stage 4의 1년·2년 공통창 timing-difference 결과를 Stage 5와 직접 연결하기**  
    그래야 “late cure-like signal”보다 “early timing difference”가 먼저인지 확인됩니다.
    
3. **Stage 6 screening 결과와 Stage 5를 나란히 놓기**  
    PNU/SNU/merged 각각이 cure modeling supportive인지 equivocal인지 unsupportive인지 확인해야 합니다. Stage 6은 이 프로젝트에서 Bayesian cure의 전단계 screening으로 설계돼 있습니다.
    
4. **Stage 9 remission sensitivity를 꼭 수행하기**  
    PNU에서 remission을 censoring으로 본 현재 결과가 얼마나 바뀌는지 확인해야 합니다.
    
5. **merged에서 site-aware 공식의 이득이 어떤 변수 패턴에서 생기는지 coefficient-level로 추적하기**  
    이번 업로드에는 coefficient table이 없으므로, 다음에는 effect table까지 같이 내는 게 좋습니다.
    

---

# J. 핵심 숫자의 쉬운 해석

## 1) PNU 1년 observed risk = 0.170

1. 숫자가 문자 그대로 뜻하는 것  
    1년 안에 transition이 일어날 누적위험이 17.0%입니다.
    
2. 쉬운 말로 바꾸면 무엇인지  
    **100명 중 약 17명**이 1년 안에 전환됩니다.
    
3. 관찰 데이터가 얼마나 직접 뒷받침하는지  
    PNU에서 1년은 primary-supported라 비교적 직접적입니다.
    
4. 그래서 어느 정도까지 믿어도 되는지  
    PNU에서 가장 믿을 만한 숫자 중 하나입니다.
    
5. 임상적/과학적으로 요구하는 해석  
    PNU는 적어도 초기 1년 위험이 가볍지 않습니다.
    

## 2) PNU 10년 predicted risk = 0.308 대 0.865

1. 숫자가 문자 그대로 뜻하는 것  
    coxph와 exponential이 10년 위험을 완전히 다르게 예측합니다.
    
2. 쉬운 말로 바꾸면 무엇인지  
    **100명 중 31명일 수도, 87명일 수도 있다는 말**이 됩니다.
    
3. 관찰 데이터가 얼마나 직접 뒷받침하는지  
    거의 직접 뒷받침되지 않습니다.
    
4. 그래서 어느 정도까지 믿어도 되는지  
    믿기보다 **projection sensitivity**로만 봐야 합니다.
    
5. 임상적/과학적으로 요구하는 해석  
    PNU 장기 꼬리로 cure를 논하는 것은 아직 이릅니다.
    

## 3) SNU 5년 observed risk = 0.211

1. 숫자가 문자 그대로 뜻하는 것  
    5년 누적 transition risk가 21.1%입니다.
    
2. 쉬운 말로 바꾸면 무엇인지  
    **100명 중 약 21명**입니다.
    
3. 관찰 데이터가 얼마나 직접 뒷받침하는지  
    secondary horizon이지만 여전히 known_count 98로 꽤 버팁니다.
    
4. 그래서 어느 정도까지 믿어도 되는지  
    장기보다 훨씬 낫고, 중기 해석용으로 충분히 쓸 수 있습니다.
    
5. 임상적/과학적으로 요구하는 해석  
    SNU는 초기보다 중기 누적이 중요합니다.
    

## 4) merged 2년 AUC, site-free 약 0.656 대 site-added 약 0.723

1. 숫자가 문자 그대로 뜻하는 것  
    site 정보를 넣으면 2년 순위구분 능력이 뚜렷이 좋아집니다.
    
2. 쉬운 말로 바꾸면 무엇인지  
    PNU와 SNU를 그냥 섞으면 누가 더 위험한지 덜 잘 맞추고, site를 알려주면 더 잘 맞춥니다.
    
3. 관찰 데이터가 얼마나 직접 뒷받침하는지  
    2년은 merged primary-supported라 강합니다.
    
4. 그래서 어느 정도까지 믿어도 되는지  
    꽤 믿어도 됩니다.
    
5. 임상적/과학적으로 요구하는 해석  
    cohort context 차이를 무시한 통합모형은 부적절할 수 있습니다.
    

## 5) merged Base 2년 threshold 0.15 net benefit, model 0.028–0.042 대 treat-all 0.004

1. 숫자가 문자 그대로 뜻하는 것  
    moderate threshold에서 model-based selection이 treat-all보다 낫습니다.
    
2. 쉬운 말로 바꾸면 무엇인지  
    무조건 다 고위험으로 보는 것보다, 모델을 쓰면 **불필요한 고위험 분류를 조금 줄일 수 있다**는 뜻입니다.
    
3. 관찰 데이터가 얼마나 직접 뒷받침하는지  
    2년 primary-supported라 비교적 직접적입니다.
    
4. 그래서 어느 정도까지 믿어도 되는지  
    방향성은 믿을 만하지만, 절대 차이는 크지 않습니다.
    
5. 임상적/과학적으로 요구하는 해석  
    모델 사용 이득은 존재하지만, 극적이지는 않습니다.
    

## 6) SNU·merged 마지막 사건 시점 = 7.57년

1. 숫자가 문자 그대로 뜻하는 것  
    그 이후에는 새 사건이 사실상 더 관찰되지 않았습니다.
    
2. 쉬운 말로 바꾸면 무엇인지  
    8–10년은 “그 뒤에 무슨 일이 더 생길지”를 데이터가 거의 새로 말해주지 못합니다.
    
3. 관찰 데이터가 얼마나 직접 뒷받침하는지  
    8년 이후는 약합니다.
    
4. 그래서 어느 정도까지 믿어도 되는지  
    10년 위험은 **숫자는 보되, projection**으로 봐야 합니다.
    
5. 임상적/과학적으로 요구하는 해석  
    SNU·merged 10년 수치는 보고는 하되, 주장 강도는 낮춰야 합니다.
    

---

# K. 이 결과가 독자에게 실제로 요구하는 해석과 주의점

독자는 이 결과를 이렇게 써야 합니다.

- **PNU는 1년 중심으로 읽고, 2년은 민감도, 그 이후는 projection으로 취급하기**
    
- **SNU와 merged는 5년까지를 실질적 해석 구간으로 우선 사용하기**
    
- **8–10년 수치는 “관찰치”보다 “모형 기반 장기 추정”으로 소개하기**
    
- **merged에서 site 효과를 treatment effect처럼 쓰지 말고 context heterogeneity의 신호로 쓰기**
    
- **PNU의 late plateau를 cure로 단정하지 말기**
    
- **Stage 5 결과만으로 cure를 말하지 말고, Stage 6과 Stage 9를 거친 뒤 Stage 10에서 통합 해석하기**
    

---

# L. 한 줄 요약

**이 결과의 핵심은 “어느 no-cure 분포가 이겼나”가 아니라, PNU의 늦은 평탄화는 주로 추적 부족의 흔적이고, SNU·merged에서는 조기 전이 timing이 cohort 맥락에 따라 실제로 다르다는 신호가 더 크다**는 점입니다.