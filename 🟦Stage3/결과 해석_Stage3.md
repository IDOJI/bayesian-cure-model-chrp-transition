# Step 1. QC summary

중대한 계산 실패는 확인되지 않았습니다.  
현재 업로드된 Stage 3 산출물은 **KM benchmark classifier 블록**으로서, 공통 backbone의 사건 정의를 그대로 사용합니다. 즉 **transition만 사건**, **right censoring과 remission은 main analysis에서 censoring**으로 처리됩니다. 또한 이 프로젝트에서 KM은 **기술적 생존곡선**이면서 동시에 **false-positive burden 비교를 위한 formal benchmark classifier** 역할을 하도록 설계되어 있습니다.

다만 해석상 주의는 분명합니다.

- **치명적 문제 없음**
    
    - 파일 구조, row 수, subject 수, benchmark registry, group risk table, classification table 사이의 숫자 일관성은 모두 맞았습니다.
        
    - KM라서 별도의 optimizer convergence code나 fit flag가 없는 것은 **정상**입니다.
        
- **경미하지만 중요한 주의점 있음**
    
    1. `merged`의 `km_site_adjusted`는 구조적으로 맞는 export입니다. 즉 **merged overall**과 **merged site-adjusted**가 함께 있는 것은 중복 오류가 아니라 **의도된 export layering**입니다.
        
    2. 다만 `merged site-adjusted`의 **PNU stratum**은 horizon support label을 `merged` 기준으로 받고 있어, **2년 시점이 primary-supported처럼 보일 수 있습니다.** 그런데 전체 프로젝트의 supported-horizon hierarchy에서는 **PNU 2년은 sensitivity**입니다. 따라서 이건 계산 오류는 아니지만, **해석 강도를 과대평가할 위험이 있는 labeling caution**입니다. 프로젝트 문서도 timing difference, follow-up immaturity, cure-like heterogeneity를 절대 섞지 말라고 분명히 적고 있습니다.
        
    3. Stage 3만으로는 cure 여부를 말하면 안 됩니다. 프로젝트 구조상 **Stage 6이 cure appropriateness screening**, **Stage 10이 최종 통합 해석**입니다.
        

# Step 2. Branching

**No major computational failure was identified.**

따라서 **코드 전체를 갈아엎을 필요는 없습니다.**  
하지만 아래의 **해석상 주의점**은 분리해서 가져가야 합니다.

- `merged site-adjusted`의 PNU 행을 **standalone PNU follow-up maturity 근거처럼 읽으면 안 됩니다.**
    
- PNU의 후반부 plateau는 **cure 신호로 읽으면 안 되고**, 우선 **짧은 follow-up과 remission 처리의 영향**으로 읽어야 합니다.
    
- 이 Stage 3 결과는 어디까지나 **KM benchmark block**이므로, 이후 Stage 4, 6, 9, 10과 연결해서 봐야 합니다.
    

---

# 0. Evidence tables

## 0-1. Data structure / cohort summary

**Source:** `stage3_km_benchmark_registry.csv`, `stage3_km_subject_horizon_predictions.csv`

### 0-1a. Benchmark registry

|dataset|benchmark_id|benchmark_scope|site_branch|assigned groups|n_subjects|
|---|--:|---|---|---|--:|
|PNU|km_overall|primary|site_free|overall|54|
|SNU|km_overall|primary|site_free|overall|208|
|merged|km_overall|primary|site_free|overall|262|
|merged|km_site_adjusted|primary_structural|site_adjusted|PNU, SNU|262|

### 0-1b. Cohort event / follow-up summary

|dataset|unique subjects|transition n|remission n|right censoring n|median follow-up years|max follow-up years|
|---|--:|--:|--:|--:|--:|--:|
|PNU|54|12|13|29|1.11|2.67|
|SNU|208|38|0|170|3.80|10.11|
|merged|262|50|13|199|2.07|10.11|

- **숫자가 직접 말하는 것**  
    PNU는 작고 짧습니다. SNU는 크고 길게 추적되었습니다. merged는 둘을 합친 자료지만, PNU와 SNU가 follow-up 구조가 매우 다릅니다.
    
- **쉬운 말로 풀면**  
    PNU는 “빨리 결과가 나는 대신 오래 본 사람은 거의 없는 코호트”, SNU는 “초반 사건은 적지만 오래 본 사람이 있는 코호트”입니다.
    
- **그래서 중요한 이유**  
    나중 꼬리 구간에서 보이는 차이가 **질병 자연경과 차이**인지, **그냥 누가 더 오래 관찰됐는지** 때문인지 분리해야 합니다.
    
- **지금 주장할 수 있는 것**  
    두 코호트는 follow-up design이 본질적으로 다르며, late tail 비교는 자동으로 공정하지 않습니다.
    
- **지금 주장하면 안 되는 것**  
    PNU의 후반 plateau를 보고 “PNU에 cure fraction이 더 크다”라고 말하면 안 됩니다.
    

## 0-2. Event and follow-up summary / follow-up support by horizon

**Source:** `stage3_km_subject_horizon_predictions.csv`

|dataset|horizon|support label|subjects observed beyond horizon|events by horizon|unknown due to early censoring|interpretation consequence|
|---|--:|---|--:|--:|--:|---|
|PNU|1y|primary-supported|28|8|18|short-horizon only reasonably supported|
|PNU|2y|sensitivity|1|12|41|almost no direct support beyond 2y|
|PNU|5y|projection|0|12|42|essentially unsupported tail|
|PNU|10y|projection|0|12|42|essentially unsupported tail|
|SNU|1y|primary-supported|163|14|31|well supported|
|SNU|2y|primary-supported|134|22|52|well supported|
|SNU|5y|secondary|64|34|110|still informative but thinning|
|SNU|10y|projection|1|38|169|almost fully tail-dependent|
|merged|1y|primary-supported|191|22|49|well supported overall|
|merged|2y|primary-supported|135|34|93|supported overall, but uneven by site|
|merged|5y|secondary|64|46|152|late information mostly from SNU|
|merged|10y|projection|1|50|211|effectively unsupported tail|

PNU는 **3년부터 10년까지 horizon을 뒷받침하는 사람이 0명**입니다. 따라서 3년 이후의 KM 위험은 “관찰된 장기 위험”이라기보다 **마지막 관찰 plateau를 그냥 유지한 값**으로 보는 것이 맞습니다.

- **숫자가 직접 말하는 것**  
    PNU는 2년 시점에도 horizon 너머까지 실제로 관찰된 사람이 1명뿐입니다. SNU와 merged는 1–2년은 괜찮지만, 5년 이후는 빠르게 빈약해집니다.
    
- **쉬운 말로 풀면**  
    “PNU는 1년까진 어느 정도 직접 본 숫자, 2년은 거의 아슬아슬, 그 뒤는 사실상 못 본 숫자”입니다.
    
- **그래서 중요한 이유**  
    Kaplan-Meier는 많은 사람이 아직 관찰 중인 구간에서 가장 안정적입니다. follow-up이 빈약한 tail에서는 숫자가 예뻐 보여도 믿음의 강도는 확 떨어집니다. Betensky도 follow-up은 KM 안정성을 해석하기 위해 반드시 같이 제시되어야 한다고 강조합니다.
    
- **지금 주장할 수 있는 것**  
    PNU는 **1년 결과가 핵심**, 2년은 민감도 수준, 3년 이상은 projection입니다. SNU와 merged는 1–2년이 핵심입니다. 이는 프로젝트의 supported-horizon hierarchy와도 일치합니다.
    
- **지금 주장하면 안 되는 것**  
    PNU의 5년, 10년 위험을 “실제로 100명 중 몇 명”처럼 단정적으로 읽으면 안 됩니다.
    

## 0-3. Prediction / time-specific risk table

**Source:** `stage3_km_group_risk_table.csv`

### 0-3a. Overall KM benchmark risk

|dataset|horizon|km risk|natural frequency|support label|
|---|--:|--:|---|---|
|PNU|1y|0.1699|약 100명 중 17명|primary-supported|
|PNU|2y|0.2988|약 100명 중 30명|sensitivity|
|PNU|5y|0.2988|약 100명 중 30명|projection|
|PNU|10y|0.2988|약 100명 중 30명|projection|
|SNU|1y|0.0728|약 100명 중 7명|primary-supported|
|SNU|2y|0.1215|약 100명 중 12명|primary-supported|
|SNU|5y|0.2108|약 100명 중 21명|secondary|
|SNU|10y|0.3265|약 100명 중 33명|projection|
|merged|1y|0.0915|약 100명 중 9명|primary-supported|
|merged|2y|0.1534|약 100명 중 15명|primary-supported|
|merged|5y|0.2391|약 100명 중 24명|secondary|
|merged|10y|0.3506|약 100명 중 35명|projection|

### 0-3b. Merged site-adjusted KM benchmark risk

|dataset|site stratum|horizon|km risk|natural frequency|support label as exported|
|---|---|--:|--:|---|---|
|merged site-adjusted|PNU|1y|0.1699|약 100명 중 17명|primary-supported|
|merged site-adjusted|PNU|2y|0.2988|약 100명 중 30명|primary-supported|
|merged site-adjusted|PNU|5y|0.2988|약 100명 중 30명|secondary|
|merged site-adjusted|PNU|10y|0.2988|약 100명 중 30명|projection|
|merged site-adjusted|SNU|1y|0.0728|약 100명 중 7명|primary-supported|
|merged site-adjusted|SNU|2y|0.1215|약 100명 중 12명|primary-supported|
|merged site-adjusted|SNU|5y|0.2108|약 100명 중 21명|secondary|
|merged site-adjusted|SNU|10y|0.3265|약 100명 중 33명|projection|

- **숫자가 직접 말하는 것**  
    1년과 2년에는 PNU 위험이 SNU보다 뚜렷하게 높습니다. 1년 절대차는 약 9.7%p, 2년 절대차는 약 17.7%p입니다.
    
- **쉬운 말로 풀면**  
    “PNU는 초반에 빨리 transition이 몰리는 코호트, SNU는 더 천천히 누적되는 코호트”처럼 보입니다.
    
- **그래서 중요한 이유**  
    이 패턴은 이후 cure model을 붙이기 전에 먼저 **timing difference**로 읽어야 합니다. 프로젝트 문서도 timing difference와 cure evidence를 분리하라고 강하게 요구합니다.
    
- **지금 주장할 수 있는 것**  
    초반 1–2년에는 PNU가 SNU보다 더 빠른 사건 누적을 보인다는 관찰은 꽤 견고합니다.
    
- **지금 주장하면 안 되는 것**  
    5년, 10년에서 PNU가 낮거나 평평해 보인다고 해서 “장기적으로 PNU가 더 안전하다”라고 말하면 안 됩니다. 그 부분은 관찰 부족이 너무 큽니다.
    

## 0-4. Threshold-based classification comparison block

**Source:** `stage3_km_benchmark_classification.csv`

### 0-4a. Merged overall KM vs merged site-adjusted KM at 1 year

|horizon|benchmark|threshold|predicted positive n|positive rate|TP|FP|TN|FN|specificity|PPV|unknown due to early censoring|
|---|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
|1y|merged overall|0.05|262|1.000|22|191|0|0|0.000|0.103|49|
|1y|merged overall|0.10|0|0.000|0|0|191|22|1.000|NA|49|
|1y|merged overall|0.15|0|0.000|0|0|191|22|1.000|NA|49|
|1y|merged overall|0.20|0|0.000|0|0|191|22|1.000|NA|49|
|1y|merged site-adjusted|0.05|262|1.000|22|191|0|0|0.000|0.103|49|
|1y|merged site-adjusted|0.10|54|0.206|8|28|163|14|0.853|0.222|49|
|1y|merged site-adjusted|0.15|54|0.206|8|28|163|14|0.853|0.222|49|
|1y|merged site-adjusted|0.20|0|0.000|0|0|191|22|1.000|NA|49|

### 0-4b. Merged overall KM vs merged site-adjusted KM at 2 years

|horizon|benchmark|threshold|predicted positive n|positive rate|TP|FP|TN|FN|specificity|PPV|unknown due to early censoring|
|---|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
|2y|merged overall|0.05|262|1.000|34|135|0|0|0.000|0.201|93|
|2y|merged overall|0.10|262|1.000|34|135|0|0|0.000|0.201|93|
|2y|merged overall|0.15|262|1.000|34|135|0|0|0.000|0.201|93|
|2y|merged overall|0.20|0|0.000|0|0|135|34|1.000|NA|93|
|2y|merged site-adjusted|0.05|262|1.000|34|135|0|0|0.000|0.201|93|
|2y|merged site-adjusted|0.10|262|1.000|34|135|0|0|0.000|0.201|93|
|2y|merged site-adjusted|0.15|54|0.206|12|1|134|22|0.993|0.923|93|
|2y|merged site-adjusted|0.20|54|0.206|12|1|134|22|0.993|0.923|93|

- **숫자가 직접 말하는 것**  
    1년 threshold 0.10 또는 0.15에서는 merged overall KM은 **아무도 고위험으로 분류하지 않는데**, merged site-adjusted KM은 **54명, 즉 전체의 20.6%를 고위험으로 분류**합니다.  
    2년 threshold 0.15 또는 0.20에서도 merged overall은 “전원 고위험” 또는 “전원 저위험”처럼 매우 거칠게 움직이지만, site-adjusted는 다시 **54명만 고위험**으로 분리합니다.
    
- **쉬운 말로 풀면**  
    “site를 섞어서 평균만 보면 중요한 초반 위험 차이가 사라지고, site를 나눠 보면 PNU만 따로 고위험처럼 잡힌다”는 뜻입니다.
    
- **그래서 중요한 이유**  
    이건 질병의 자연경과 차이일 수도 있지만, 프로젝트 문서상 site는 치료 그 자체가 아니라 **care pathway, selection, follow-up 구조의 proxy**로 읽어야 합니다. 즉 평균을 내서 섞어버리면 초기 위험 이질성을 가려버릴 수 있습니다.
    
- **지금 주장할 수 있는 것**  
    merged cohort를 하나의 균질한 집단처럼 취급한 pooled KM benchmark는 초반 분류 의사결정을 왜곡할 수 있습니다.
    
- **지금 주장하면 안 되는 것**  
    2년 threshold 0.15와 0.20에서 site-adjusted PPV가 0.923이라고 해서, “PNU는 2년 고위험 예측이 거의 완벽하다”라고 말하면 안 됩니다. 왜냐하면 고위험으로 분류된 54명 중 **41명은 2년 시점 outcome이 미확정**이기 때문입니다. 즉 예뻐 보이는 숫자 상당 부분이 censoring 구조에 기대고 있습니다.
    

## 0-5. Diagnostics / convergence / assumption-check table

**Source:** uploaded Stage 3 CSV bundle 전체 재계산 QC

|QC item|result|interpretation|
|---|---|---|
|Export manifest uniqueness|pass|8개 파일 구조가 서로 충돌 없이 존재|
|Benchmark registry uniqueness|pass|4개 benchmark row가 중복 없이 정의됨|
|Subject-horizon row construction|pass|각 benchmark에서 subject × 10 horizons 구조가 정확함|
|Group-risk reconstruction|pass|group risk table가 subject predictions와 정확히 일치|
|Classification reconstruction|pass|TP, FP, TN, FN, PPV, specificity 등 전부 재계산 일치|
|Subject-count consistency across files|pass|registry, predictions, group_risk, classification 간 subject 수 일치|
|Convergence / fit flags|not applicable|KM는 iterative optimizer가 아니라 별도 convergence code가 없음|
|Minor structural caution|caution|merged site-adjusted PNU 2y support label은 해석 시 standalone PNU보다 더 조심해야 함|

- **숫자가 직접 말하는 것**  
    계산상 내부 모순은 찾지 못했습니다.
    
- **쉬운 말로 풀면**  
    “숫자가 틀려서 이상하게 나온 Stage 3”은 아닙니다.
    
- **그래서 중요한 이유**  
    지금부터의 핵심은 코드를 고치는 일이 아니라, **이 숫자를 어떻게 읽느냐**입니다.
    
- **지금 주장할 수 있는 것**  
    Stage 3는 benchmark KM block으로서 구조적으로 신뢰할 만합니다.
    
- **지금 주장하면 안 되는 것**  
    내부 일관성이 좋다고 해서 과학적 해석까지 강하다고 말하면 안 됩니다. 해석 강도는 follow-up support에 달려 있습니다.
    

## 0-6. Parameter / coefficient table

**Source:** 해당 없음

Stage 3는 **비모수적 KM benchmark stage**이므로 회귀계수나 cure fraction parameter 표는 없습니다.  
이건 누락이 아니라 **stage의 설계상 정상**입니다. KM은 이 프로젝트에서 benchmark이지 individualized model이 아닙니다.

- **숫자가 직접 말하는 것**  
    계수가 없는 stage입니다.
    
- **쉬운 말로 풀면**  
    “누가 왜 위험한지”를 추정하는 단계가 아니라, “관찰된 코호트 평균 위험을 기준점으로 놓는 단계”입니다.
    
- **그래서 중요한 이유**  
    Stage 3 숫자를 가지고 성별, 나이, site effect의 독립적 인과 해석을 하면 안 됩니다.
    
- **지금 주장할 수 있는 것**  
    Stage 3는 baseline descriptive/benchmark evidence입니다.
    
- **지금 주장하면 안 되는 것**  
    Stage 3만으로 공변량 효과를 결론내리면 안 됩니다.
    

## 0-7. Cohort-specific comparison table

**Source:** `stage3_km_group_risk_table.csv`, `stage3_km_subject_horizon_predictions.csv`

|comparison|1y risk|2y risk|support strength|direct reading|
|---|--:|--:|---|---|
|PNU overall|0.1699|0.2988|1y supported, 2y weak|early event concentration|
|SNU overall|0.0728|0.1215|both supported|slower early accumulation|
|PNU minus SNU|+0.0971|+0.1773|strongest interpretable contrast is early|timing difference likely precedes any cure claim|
|merged overall|0.0915|0.1534|supported overall|pooled average hides site heterogeneity|
|merged site-adjusted PNU|0.1699|0.2988|exported as merged-primary, but interpret PNU 2y cautiously|mirrors standalone PNU risk|
|merged site-adjusted SNU|0.0728|0.1215|supported|mirrors standalone SNU risk|

- **숫자가 직접 말하는 것**  
    site-adjusted merged는 사실상 standalone PNU와 standalone SNU의 early-risk 차이를 그대로 드러냅니다.
    
- **쉬운 말로 풀면**  
    “섞어놓고 평균내면 안 보이던 차이가, site를 나누면 다시 보인다”는 뜻입니다.
    
- **그래서 중요한 이유**  
    이것이 바로 later stages에서 **timing difference를 먼저 떼어놓고**, 그 다음에야 cure-like heterogeneity를 묻도록 설계된 이유입니다.
    
- **지금 주장할 수 있는 것**  
    PNU–SNU 차이는 late tail에서 갑자기 생긴 것이 아니라, **이미 1년과 2년 공통 창(window)에서 보이기 시작한 차이**입니다.
    
- **지금 주장하면 안 되는 것**  
    이 차이를 곧바로 생물학적 차이 또는 치료효과라고 부르면 안 됩니다. site는 현재 치료정보가 없어서 proxy로만 해석해야 합니다.
    

---

# A. Scientific question reframed

이 Stage 3 결과가 진짜 답하려는 질문은  
“**PNU, SNU, merged에서 observed transition risk의 시간 패턴이 어떻게 다른가, 그리고 그 차이가 이후 cure-like 해석 전에 이미 보이는 timing difference인지**”입니다.

즉 지금 단계의 핵심은  
“어느 모델이 더 좋으냐”가 아니라,  
“**관찰된 KM benchmark만 놓고 봐도 초기 위험 차이, follow-up 빈약성, remission 처리 문제를 어떻게 읽어야 하는가**”입니다. 프로젝트 문서도 timing difference, follow-up immaturity, cure-like heterogeneity를 섞지 말라고 못 박고 있습니다.

# B. Bottom-line integrated interpretation

- **PNU는 early-transition cohort처럼 보입니다.** 1년 위험이 약 17%, 2년 위험이 약 30%로 SNU보다 훨씬 높고, 이 차이는 late tail이 아니라 이미 초반 공통 창에서 관찰됩니다.
    
- **SNU는 slower-accumulating cohort처럼 보입니다.** 1–2년 위험은 낮지만 follow-up이 길어서 시간이 지나며 위험이 더 누적됩니다.
    
- **PNU의 late plateau는 cure 증거가 아니라 follow-up collapse의 결과로 읽는 것이 더 타당합니다.** 2년 이후 관찰 지원이 거의 없고 3년 이후는 사실상 직접 관찰이 없습니다.
    
- **merged overall KM는 초반 이질성을 과하게 평균내어 숨깁니다.** 반대로 merged site-adjusted KM는 초반 site 차이를 드러내므로, early benchmark로는 더 개념적으로 맞습니다.
    
- **PNU에서 remission이 censoring으로 처리되는 설계는 late low-risk appearance를 강화할 수 있습니다.** 따라서 PNU의 “안정적인 장기 비전이”처럼 보이는 모습은 Stage 9 경쟁위험 민감도 분석 전까지는 보수적으로 읽어야 합니다.
    

# C. Statistical interpretation

이 결과가 통계적으로 가장 강하게 지지하는 것은 **초기 1–2년 timing difference**입니다.

1년 위험은 PNU 16.99%, SNU 7.28%입니다. 문자 그대로는 **100명 중 약 17명 대 7명**입니다. 이 차이는 직접 관찰 구간에서 나온 값이어서 비교적 견고합니다.  
2년 위험은 PNU 29.88%, SNU 12.15%입니다. 문자 그대로는 **100명 중 약 30명 대 12명**입니다. 하지만 여기서 PNU는 2년 넘게 실제로 관찰된 사람이 1명뿐이므로, 값 자체는 남아 있어도 **지지 강도는 약합니다**.

반면 5년 이후의 PNU risk 29.88% plateau는 숫자만 보면 “위험이 더 이상 안 오른다”입니다. 쉬운 말로는 “2년쯤 지나면 남은 사람들은 거의 transition 안 하는 집단처럼 보인다”는 해석이 가능해 보입니다. 하지만 그 숫자는 **거의 관찰 지원이 없는 tail carry-forward**입니다. 따라서 이 숫자가 요구하는 통계적 태도는 **treat as descriptive only**, 더 정확히는 **do not over-interpret**입니다.

merged overall과 merged site-adjusted의 차이도 중요합니다. 1년 threshold 0.10에서 merged overall은 아무도 고위험으로 분류하지 않지만, merged site-adjusted는 54명, 즉 전체의 20.6%를 고위험으로 분류합니다. 이 숫자가 현실에서 뜻하는 바는, **site를 무시하고 풀링하면 실제 존재하는 초기 위험 이질성이 benchmark에서 사라진다**는 것입니다. 따라서 later individualized models와 비교할 baseline benchmark는 pooled overall 하나만 쓰면 불공정할 수 있습니다. 프로젝트 문서도 incompatible scale 비교를 하면 안 된다고 명시합니다.

# D. Clinical interpretation

임상 맥락이 현재 업로드된 Stage 3 파일만으로 충분히 구체적이지는 않지만, 최소한 다음 정도는 방어적으로 말할 수 있습니다.

첫째, **질병의 자연경과가 한 가지 tempo로만 진행되는 것 같지 않습니다.**  
PNU는 “초반에 event가 몰리는 경로”가 더 강하고, SNU는 “초기 위험은 낮지만 더 길게 event가 쌓이는 경로”가 더 강합니다. 쉬운 말로 하면, 어떤 코호트는 **빨리 나빠지는 사람들**이 더 많고, 어떤 코호트는 **천천히 누적되는 사람들**이 더 많아 보입니다.

둘째, 이 차이를 곧바로 **생물학적 본질 차이**로 읽으면 안 됩니다.  
프로젝트 문서상 site effect는 현재 치료 데이터가 없어서 **treatment context, care pathway, referral structure, selection, follow-up maturity의 proxy**로 해석해야 합니다. 즉 PNU가 더 “위험한 질병”이라기보다, **더 높은 위험 환자가 먼저 모였거나, 더 이른 ascertainment가 되었거나, 추적구조가 달랐을 가능성**이 큽니다.

셋째, **outcome definition 자체가 임상 해석을 바꿉니다.**  
현재 main backbone은 transition만 event로 보고 remission은 censoring으로 둡니다. 그런데 데이터 사전에 따르면 remission은 PNU에만 실제로 존재합니다. 따라서 PNU의 late low-risk appearance는 “정말 transition-capable하지 않은 집단” 때문일 수도 있지만, 동시에 “remission이라는 다른 경로로 risk set에서 빠져나간 사람들” 때문일 수도 있습니다. 쉬운 말로 하면, **좋아져서 빠진 사람을 그냥 추적중단처럼 다루면, 남은 곡선이 인위적으로 평평해 보일 수 있다**는 뜻입니다. 그래서 이 Stage 3 결과는 Stage 9를 반드시 요구합니다.

# E. Model-assumption-to-clinical-reality mapping

## 1) Overall KM benchmark

- **핵심 가정**  
    개별화 없이 코호트 평균 위험을 기술하고, censoring이 event process와 심하게 뒤엉키지 않았다고 보는 것입니다.
    
- **여기서 얼마나 그럴듯한가**  
    PNU 1년, SNU 1–2년, merged 1–2년에서는 비교적 그럴듯합니다.
    
- **무엇이 약화시키는가**  
    PNU의 짧은 follow-up, late sparse risk set, 그리고 remission censoring 처리입니다.
    
- **쉬운 말**  
    “짧게는 믿을 만한 평균, 길게는 조심해야 하는 평균”입니다.
    

## 2) Merged site-adjusted KM benchmark

- **핵심 가정**  
    merged를 한 덩어리로 보지 않고, site별 KM를 subject에게 할당해 structural heterogeneity를 보존합니다.
    
- **여기서 얼마나 그럴듯한가**  
    early merged classification에는 pooled overall보다 훨씬 그럴듯합니다. 실제로 1년과 2년의 PNU–SNU 차이를 되살립니다.
    
- **무엇이 약화시키는가**  
    site-specific maturity를 merged label이 충분히 표현하지 못하는 부분, 특히 PNU 2년 label입니다.
    
- **쉬운 말**  
    “평균내서 흐려진 차이를 다시 보이게 하는 benchmark”입니다.
    

## 3) Thresholded KM classification

- **핵심 가정**  
    같은 group에 속하면 같은 risk를 적용받는다는 매우 단순한 classifier입니다.
    
- **여기서 얼마나 그럴듯한가**  
    individualized truth가 아니라 **benchmark**로서는 좋습니다. 프로젝트 문서도 KM을 benchmark로 위치시킵니다.
    
- **무엇이 약화시키는가**  
    threshold 근처에서 pooled overall은 “전원 high” 또는 “전원 low”처럼 너무 거칠게 움직입니다.
    
- **쉬운 말**  
    “정밀한 진단도구가 아니라, 앞으로 나올 individualized 모델과 비교하기 위한 기준점”입니다.
    

# F. Why the observed pattern matters

이 Stage 3에서 가장 publishable한 포인트는 **“late tail보다 early benchmark가 이미 site 구조의 중요성을 보여준다”**는 점입니다.

많은 survival 연구에서는 꼬리에서 plateau가 보이면 곧바로 cure를 떠올립니다. 그런데 여기서는 그보다 먼저, **1년과 2년이라는 supported horizon에서조차 pooled overall benchmark와 site-adjusted benchmark가 완전히 다른 분류 결정을 내립니다.**  
즉 이 데이터의 핵심 이야기는 “누가 cure인가”보다 먼저, **“누가 초반에 빨리 transition하는가, 그리고 그 차이가 cohort/site 구조를 평균내면 얼마나 가려지는가”**입니다.

이게 중요한 이유는, 나중에 cure model이 false-positive burden을 줄였다고 나와도 그것이 진짜 cure-like heterogeneity를 잡은 것인지, 아니면 그냥 **초기 cohort heterogeneity를 더 잘 반영한 것인지** 구분해야 하기 때문입니다. 프로젝트 문서가 Stage 4와 Stage 6을 Stage 7/8보다 먼저 두는 이유가 바로 이것입니다.

# G. Alternative explanations and threats to interpretation

## Strongly supported

1. **Timing difference**
    
    - PNU가 SNU보다 초반 transition이 빠르다.
        
    - 1년과 2년 공통 창에서 이미 관찰된다.
        
2. **Late-tail instability**
    
    - 특히 PNU는 2년 이후 직접 관찰 지원이 거의 없다.
        
    - 따라서 plateau는 신호라기보다 정보 부족의 산물일 가능성이 크다.
        

## Plausible

3. **Site as care-pathway / selection proxy**
    
    - 치료 강도, referral pattern, baseline severity mix, ascertainment timing 차이가 섞여 있을 수 있다.
        
4. **Remission handling artifact**
    
    - PNU에서 remission이 censoring으로 들어가 late risk를 인위적으로 낮춰 보일 수 있다.
        

## Speculative

5. **True cure-like heterogeneity**
    
    - 일부 비전이성 subgroup이 존재할 수는 있다.
        
    - 하지만 Stage 3 alone으로는 방어되지 않는다.
        
6. **Intrinsic biological site difference**
    
    - 가능성은 있으나, 현재 치료/프로토콜 정보 없이 주장하기 어렵다.
        

이 competing explanation들을 가르는 데 가장 중요한 것은

- Stage 4의 1년/2년 common-window timing analysis,
    
- Stage 6의 cure appropriateness screening,
    
- Stage 9의 remission competing-risk sensitivity,
    
- Stage 10의 false-positive burden 및 net benefit 통합 비교입니다.
    

# H. Manuscript-ready discussion paragraph

In this Kaplan-Meier benchmark analysis, the most robust signal was not late cure-like flattening but an early timing difference between cohorts. PNU showed substantially higher 1-year and 2-year transition risk than SNU, indicating a more front-loaded event pattern, whereas SNU accrued risk more gradually over longer follow-up. Importantly, this difference was already visible within supported early horizons and therefore should be interpreted before invoking cure-like heterogeneity. By contrast, the apparent plateau in PNU beyond 2 years was weakly supported by observed follow-up and is better regarded as a tail-unstable feature than as direct evidence of a non-susceptible subgroup. The pooled merged benchmark also obscured meaningful early site heterogeneity, whereas the site-adjusted benchmark preserved it, suggesting that cohort context materially affects risk communication even before individualized modeling. Because remission occurs structurally in PNU and is treated as censoring in the main backbone, any apparent long-term attenuation of transition risk in PNU should remain provisional until remission is reanalyzed as a competing event. Overall, these findings argue that subsequent cure-model interpretation should be anchored to early supported horizons, explicitly separated from timing differences and follow-up immaturity, and judged by whether cure structures reduce false-positive burden and improve clinical usefulness rather than by tail shape alone.

# I. Next-step analyses

1. **Stage 4 공통 1년, 2년 restricted timing-difference 분석을 바로 연결**  
    초반 차이가 hazard timing 차이인지, 단순 누적위험 차이인지 분리해야 합니다.
    
2. **Stage 9 remission competing-risk sensitivity를 PNU 중심으로 우선 실행**  
    PNU late plateau가 remission censoring artifact인지 확인해야 합니다.
    
3. **Stage 6 cure-appropriateness screening 결과와 반드시 결합**  
    RECeUS-AIC, HSU family set, Xie contradiction check 없이 late cure 해석을 강하게 하면 안 됩니다.
    
4. **Stage 5 no-cure individualized model과 동일 horizon/threshold grid에서 비교**  
    Stage 3 benchmark가 보여준 초반 site heterogeneity를 individualized model이 어떻게 흡수하는지 봐야 합니다.
    
5. **Merged에서 pooled overall과 site-adjusted를 모두 유지하되, early primary claim은 site-adjusted 중심으로 보고**  
    pooled overall은 sensitivity benchmark로 두는 것이 더 안전합니다.
    

# J. 핵심 숫자의 쉬운 해석

## 1) PNU 1년 risk = 0.1699

1. **숫자가 문자 그대로 뜻하는 것**  
    1년 내 transition 누적위험이 16.99%입니다.
    
2. **쉬운 말로 바꾸면 무엇인지**  
    PNU에서는 100명 중 약 17명 정도가 1년 안에 transition한 셈입니다.
    
3. **관찰 데이터가 얼마나 직접 뒷받침하는지**  
    1년 넘게 관찰된 사람이 28명 있어, PNU 안에서는 비교적 직접적으로 지지됩니다.
    
4. **그래서 어느 정도 믿어도 되는지**  
    PNU에서 가장 믿을 만한 숫자 중 하나입니다.
    
5. **임상적/과학적으로 요구하는 해석**  
    PNU는 초반 monitoring intensity를 높여야 하는 cohort일 가능성이 큽니다.
    

## 2) SNU 1년 risk = 0.0728

1. **숫자가 문자 그대로 뜻하는 것**  
    1년 내 transition 누적위험이 7.28%입니다.
    
2. **쉬운 말로 바꾸면 무엇인지**  
    SNU에서는 100명 중 약 7명 정도가 1년 안에 transition합니다.
    
3. **관찰 데이터가 얼마나 직접 뒷받침하는지**  
    1년 넘게 관찰된 사람이 163명이라 매우 잘 지지됩니다.
    
4. **그래서 어느 정도 믿어도 되는지**  
    초반 benchmark로 꽤 신뢰할 수 있습니다.
    
5. **임상적/과학적으로 요구하는 해석**  
    SNU는 PNU보다 초반 급격한 위험은 낮지만, 더 길게 누적되는 패턴을 볼 준비가 필요합니다.
    

## 3) PNU minus SNU 2년 risk difference = +0.1773

1. **숫자가 문자 그대로 뜻하는 것**  
    2년 누적위험 차이가 17.73%p입니다.
    
2. **쉬운 말로 바꾸면 무엇인지**  
    2년 기준으로 보면 PNU가 SNU보다 100명 중 약 18명 더 많이 transition한 셈입니다.
    
3. **관찰 데이터가 얼마나 직접 뒷받침하는지**  
    SNU 쪽은 잘 지지되지만, PNU 쪽은 2년 이후 관찰이 거의 없어 지지 강도는 약해집니다.
    
4. **그래서 어느 정도 믿어도 되는지**  
    방향성은 믿되, 정밀한 소수점 수준 비교는 조심해야 합니다.
    
5. **임상적/과학적으로 요구하는 해석**  
    초반 site 차이는 real signal일 가능성이 크지만, 크기 자체는 follow-up 구조의 영향을 받습니다.
    

## 4) merged overall 1년 threshold 0.10에서 predicted positive = 0명

1. **숫자가 문자 그대로 뜻하는 것**  
    pooled merged KM benchmark는 아무도 고위험으로 분류하지 않습니다.
    
2. **쉬운 말로 바꾸면 무엇인지**  
    전체를 평균내면 “1년 10% 넘는 위험자는 없다”고 말하는 셈입니다.
    
3. **관찰 데이터가 얼마나 직접 뒷받침하는지**  
    계산은 정확하지만, cohort heterogeneity를 평균으로 지워버린 결과입니다.
    
4. **그래서 어느 정도 믿어도 되는지**  
    pooled average로서만 의미가 있고, 실제 분류 기준으로는 매우 거칠 수 있습니다.
    
5. **임상적/과학적으로 요구하는 해석**  
    merged cohort를 균질 집단처럼 보는 benchmark는 early decision support에 부적절할 수 있습니다.
    

## 5) merged site-adjusted 1년 threshold 0.10에서 predicted positive = 54명

1. **숫자가 문자 그대로 뜻하는 것**  
    site-adjusted benchmark는 54명을 고위험으로 분류합니다.
    
2. **쉬운 말로 바꾸면 무엇인지**  
    평균을 내지 않으면 “PNU 쪽 사람들만 따로 보면 1년 10%를 넘는 고위험 집단이 있다”는 뜻입니다.
    
3. **관찰 데이터가 얼마나 직접 뒷받침하는지**  
    1년 시점이라 비교적 직접 지지됩니다.
    
4. **그래서 어느 정도 믿어도 되는지**  
    early structural heterogeneity를 보여주는 benchmark로는 믿을 만합니다.
    
5. **임상적/과학적으로 요구하는 해석**  
    이후 모델 비교에서 “cure model이 나아졌다”는 말은, 먼저 이 early site heterogeneity를 얼마나 잘 반영했는지와 구분해서 봐야 합니다.
    

# K. 이 결과가 독자에게 실제로 요구하는 해석과 주의점

독자는 이 결과를 이렇게 써야 합니다.

- **PNU는 1년 중심으로 읽고, 2년은 보조적, 3년 이상은 projection으로만 취급**해야 합니다.
    
- **SNU와 merged는 1–2년을 핵심 해석 구간**으로 삼아야 합니다.
    
- **late plateau를 cure 신호로 과잉해석하지 말고**, follow-up maturity와 risk-set collapse를 먼저 보여줘야 합니다.
    
- **merged overall benchmark만으로 의사결정 논리를 만들지 말고**, early horizon에서는 **site-adjusted benchmark를 함께 제시**해야 합니다.
    
- **PNU의 remission 구조를 무시한 장기 transition 해석은 금물**입니다.
    
- later stages에서 long-horizon risk는 **“직접 관찰된 빈도”가 아니라 “unsupported tail description 또는 projection”**이라고 명시해야 합니다.
    

# L. 한 줄 요약

**이 Stage 3 결과가 가장 강하게 말해주는 것은 “cure”가 아니라 “PNU는 초반에 더 빨리 transition이 몰리고, PNU의 late plateau는 follow-up 부족과 remission 처리 때문에 생긴 착시일 가능성이 크다”는 점입니다.**