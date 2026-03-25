# 0. QC summary

중대한 계산 실패는 확인되지 않았습니다.  
이번 Stage 2는 원래부터 **“cure claim 전에 follow-up maturity를 기술하는 단계”**로 설계되어 있고, 메타데이터에도 **`model_fitting_allowed = FALSE`**가 명시되어 있습니다. 따라서 사용자가 요구한 QC 항목 중 **convergence code/message, fit flag, fitted-model registry 불일치**는 이 단계에서는 **점검 대상 자체가 아닙니다**. 이 단계의 핵심은 KM, reverse-KM, numbers-at-risk, completeness, site contribution의 구조적 일관성과 해석 가능성을 점검하는 것입니다.

내가 확인한 범위에서는 구조적 QC가 잘 되어 있습니다. 업로드된 `stage2_validation_summary.csv`에는 13개 검사가 모두 pass로 표시되어 있고, `stage2_validation_issue_rows.csv`는 비어 있습니다. 코드 수준에서도 Stage 2는 numbers-at-risk 단조성, composition 합 1, horizon monotonicity, merged site consistency, raw reconstruction 등을 점검하도록 설계돼 있습니다. 따라서 **치명적인 코딩/집계 오류의 흔적은 현재 파일들에서는 보이지 않습니다**.

다만 **사소하지만 분명한 패키지 단위 주의점**은 있습니다. Run log와 manifest상 Stage 2는 19개 파일을 export했고 그 안에 `stage2_followup_bundle.rds`도 포함돼야 하는데, 이번에 내가 직접 열어본 업로드 묶음에는 그 `.rds`가 빠져 있습니다. 이것은 **Stage 2 계산 실패라기보다 “검토용 업로드 패키지의 일부 누락”**에 가깝고, 다행히 핵심 source-of-truth CSV들은 모두 있어 해석 자체는 가능합니다.

또 하나 중요합니다. **행 수가 subject 수와 안 맞는 것 자체는 Stage 2에서 정상**입니다. 왜냐하면 이 산출물은 한 사람 한 행 자료가 아니라,

- overall, sex, site subgroup,
    
- transition_only / remission_sensitive view,
    
- merged의 site_free / site_adjusted layering,
    
- 1–10년 horizon 반복  
    을 쌓아 둔 **layered export**이기 때문입니다. 이 layering은 metadata와 코드에 명시된 설계이며, 중복 집계 버그의 증거가 아닙니다. 특히 merged overall/sex는 site_free와 site_adjusted 두 구조를 모두 보존하고, merged site-stratified는 site_adjusted만 쓰도록 되어 있습니다.
    

---

# Step 2. Branching

**No major computational failure was identified.**

따라서 downstream interpretation은 가능합니다.  
하지만 아래 해석은 **“계산적으로 안전하다”**는 뜻이지, **“과학적으로 확정적이다”**는 뜻은 아닙니다.  
특히 이번 결과는 본질적으로 **follow-up maturity 결과**이므로, 해석상 주의점은 별도로 아주 강하게 남습니다.

---

# 0. Evidence tables

## 0-1. Data structure / cohort summary

**Source:** `stage2_numbers_at_risk.csv` + `stage2_followup_composition.csv` + `stage2_reverse_km_summary.csv`  
**보여주는 것:** 코호트 크기, 10년 시점 상태 분포, reverse-KM 중앙 추적기간

|코호트|대상자수|전이|remission|우측검열|10년 시점 미해결|역KM 중앙추적(전이만 사건, 년)|역KM 중앙추적(remission 민감, 년)|최대 관찰추적(년)|
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
|PNU|54|12|13|29|0|1.37|1.57|2.67|
|SNU|208|38|0|169|1|4.91|4.91|10.11|
|merged|262|50|13|198|1|3.83|4.21|10.11|

**숫자가 직접 말하는 것**  
PNU는 전체 54명 중 10년 시점에 전이 12명, remission 13명, 우측검열 29명이고, reverse-KM 중앙 추적기간은 약 1.4–1.6년입니다. 반면 SNU는 208명, 전이 38명, remission 0명, 우측검열 169명, reverse-KM 중앙 추적기간 약 4.9년입니다. merged는 단순 합처럼 보이지만 추적기간 중앙값은 3.8–4.2년으로 SNU 쪽에 더 끌려 있습니다.

**쉬운 말로 풀면**  
PNU는 “대부분의 일이 처음 1–2년 안에 벌어지고 빨리 관찰이 끝나는 코호트”이고, SNU는 “더 오래 관찰돼서 늦게 생기는 일도 어느 정도 볼 수 있는 코호트”입니다.

**그래서 중요한 이유**  
같은 disease process라도, 한 코호트는 짧게 보고 다른 코호트는 길게 보면 나중 꼬리 모양이 달라 보입니다. 이 차이를 cure, plateau, 장기안정으로 오해하면 안 됩니다.

**지금 주장할 수 있는 것**  
PNU와 SNU는 **follow-up structure가 매우 다르다**. 따라서 long-horizon 비교에서 두 코호트를 같은 무게로 읽으면 안 됩니다.

**지금 주장하면 안 되는 것**  
“PNU는 진짜로 장기 plateau가 있다”, “SNU가 biologically 더 chronic하다”, “merged가 세 코호트의 평균적 자연경과를 잘 보여준다” 같은 말은 아직 할 수 없습니다.

---

## 0-2. Follow-up support by horizon

**Source:** `stage2_followup_horizon_summary.csv`  
**보여주는 것:** overall, transition_only 기준의 horizon support와 completeness  
**주의:** merged는 metadata rule에 따라 **site_adjusted**를 주보고용으로 사용해야 합니다.

|코호트|시점(년)|지원등급|시점 이후 관찰중 n|총 n|시점 이후 관찰중 비율(%)|잠재적으로 완전추적 가능 n|잠재적 완전추적 가능 비율(%)|Percentage(%)|CCI(%)|SPT(%)|플롯상태|
|:--|--:|:--|--:|--:|--:|--:|--:|--:|--:|--:|:--|
|PNU|1|primary_supported|28|54|51.9|47|87.0|76.6|86.2|88.0|stable|
|PNU|2|sensitivity|1|54|1.9|35|64.8|31.4|82.9|86.2|very_sparse|
|PNU|5|projection|0|54|0.0|0|0.0|NA|NA|NA|masked_no_eligible|
|PNU|10|projection|0|54|0.0|0|0.0|NA|NA|NA|masked_no_eligible|
|SNU|1|primary_supported|163|208|78.4|191|91.8|91.1|92.6|92.8|stable|
|SNU|2|primary_supported|134|208|64.4|175|84.1|86.9|89.5|90.1|stable|
|SNU|5|secondary|64|208|30.8|146|70.2|63.7|79.2|81.9|stable|
|SNU|10|projection|1|208|0.5|64|30.8|29.7|54.1|64.3|very_sparse|
|merged|1|primary_supported|191|262|72.9|238|90.8|88.2|91.4|91.9|stable|
|merged|2|primary_supported|135|262|51.5|210|80.2|77.6|88.5|89.5|stable|
|merged|5|secondary|64|262|24.4|146|55.7|63.7|79.2|81.9|stable|
|merged|10|projection|1|262|0.4|64|24.4|29.7|54.1|64.3|very_sparse|

**숫자가 직접 말하는 것**  
PNU는 1년까지만 primary-supported이고, 2년에는 risk set이 1명뿐입니다. SNU와 merged는 1–2년은 비교적 안정적이고, 5년은 secondary, 10년은 사실상 tail 1명 수준입니다. 이 지원등급 체계는 프로젝트 스펙과도 일치합니다.

**쉬운 말로 풀면**  
PNU는 “1년 정도까지만 자신 있게 말할 수 있고”, SNU와 merged는 “2년은 꽤 괜찮고 5년은 조심스럽게, 10년은 거의 참고용”입니다.

**그래서 중요한 이유**  
나중에 Stage 7, 8, 10에서 어떤 cure/no-cure model이든 5년, 10년 수치를 내더라도, 그 숫자가 관찰에 얼마나 기대고 있는지 먼저 알아야 합니다. follow-up이 부족하면 model이 tail을 마음대로 그릴 위험이 커집니다. 이런 문제는 cure model에서 특히 심해질 수 있습니다.

**지금 주장할 수 있는 것**  
PNU의 2년 이후, 그리고 SNU/merged의 10년은 **직접 관찰이 아니라 tail 의존 해석 영역**입니다.

**지금 주장하면 안 되는 것**  
PNU 5년, 10년 plateau를 근거로 장기 저위험군이나 cure fraction을 강하게 말하면 안 됩니다.

---

## 0-3. Outcome-definition sensitivity: transition_only vs remission_sensitive

**Source:** `stage2_followup_curve_data.csv` + `stage2_followup_side_by_side.csv`  
**보여주는 것:** overall risk가 결과 정의에 얼마나 민감한지

|코호트|시점(년)|remission 민감 위험(%)|전이만 사건 위험(%)|두 정의 차이(%p)|
|:--|--:|--:|--:|--:|
|PNU|1|24.1|17.0|7.1|
|PNU|2|67.8|29.9|38.0|
|PNU|5|100.0|29.9|70.1|
|PNU|10|100.0|29.9|70.1|
|SNU|1|7.3|7.3|0.0|
|SNU|2|12.2|12.2|0.0|
|SNU|5|21.1|21.1|0.0|
|SNU|10|32.6|32.6|0.0|
|merged|1|10.5|9.2|1.3|
|merged|2|21.2|15.3|5.9|
|merged|5|29.8|23.9|5.9|
|merged|10|40.1|35.1|5.0|

**숫자가 직접 말하는 것**  
PNU에서는 2년 위험이 결과 정의에 따라 29.9% vs 67.8%로 **무려 38.0%p** 차이 납니다. SNU는 remission이 구조적으로 없어 두 정의가 완전히 같습니다. merged는 PNU가 섞이면서 초반에만 차이가 조금 남습니다.

**쉬운 말로 풀면**  
PNU에서는 remission을 “그냥 관찰 중단처럼 취급”할지, 아니면 “임상적으로 중요한 다른 경로”로 볼지에 따라 완전히 다른 그림이 나옵니다. SNU는 그런 문제가 거의 없습니다.

**그래서 중요한 이유**  
이건 결과 정의가 disease process 해석을 바꿀 수 있다는 뜻입니다. PNU의 plateau 비슷한 모습 일부는 “전이가 안 생겨서”가 아니라, **remission이 많이 발생해서 전이 위험 집합에서 빠져나간 결과**일 수 있습니다.

**지금 주장할 수 있는 것**  
PNU에서는 **transition-only라는 주 결과 정의의 가정이 SNU보다 훨씬 더 취약**합니다.

**지금 주장하면 안 되는 것**  
PNU의 전이만 사건 곡선만 보고 “장기적으로 안전하다”거나 “cure-like subgroup가 보인다”고 말하면 안 됩니다. 먼저 remission competing path를 따로 봐야 합니다.

---

## 0-4. Merged structure: site_free vs site_adjusted, and late-tail site dominance

**Source:** `stage2_followup_horizon_summary.csv` + `stage2_merged_site_contribution.csv`  
**보여주는 것:** merged long-horizon 해석이 왜 site-adjusted를 써야 하는지

### 0-4a. merged overall completeness comparison

|구조|시점(년)|완전추적 잠재 n|완전추적 잠재 비율(%)|시점 이후 관찰중 n|Percentage(%)|CCI(%)|SPT(%)|
|:--|--:|--:|--:|--:|--:|--:|--:|
|site_adjusted|1|238|90.8|191|88.2|91.4|91.9|
|site_free|1|255|97.3|191|83.5|89.4|90.0|
|site_adjusted|2|210|80.2|135|77.6|88.5|89.5|
|site_free|2|243|92.7|135|68.7|83.5|84.9|
|site_adjusted|5|146|55.7|64|63.7|79.2|81.9|
|site_free|5|204|77.9|64|48.0|71.2|74.6|
|site_adjusted|10|64|24.4|1|29.7|54.1|64.3|
|site_free|10|144|55.0|1|22.9|49.3|58.0|

### 0-4b. merged overall late-tail site contribution

|시점(년)|전체 eligible n|기여 site 수|지배 site|지배 site 비율(%)|상태|
|--:|--:|--:|:--|--:|:--|
|1|238|2|SNU|80.3|single_site_dominant|
|2|210|2|SNU|83.3|single_site_dominant|
|5|146|1|SNU|100.0|single_site_only|
|10|64|1|SNU|100.0|single_site_only|

**숫자가 직접 말하는 것**  
merged overall에서 10년 잠재적 완전추적 가능 수는 site_free 144명, site_adjusted 64명으로 큰 차이가 납니다. 그리고 site_adjusted 기준 5년 이후 eligible set은 **100% SNU**입니다.

**쉬운 말로 풀면**  
merged를 그냥 한 덩어리로 오래 따라간 것처럼 보면 follow-up이 훨씬 좋아 보이지만, site별 행정 종료일을 반영하면 그 착시가 크게 줄어듭니다. 그리고 5년 이후 merged 꼬리는 사실상 SNU만 남습니다.

**그래서 중요한 이유**  
late horizon에서 merged를 “PNU와 SNU를 섞은 평균 자연경과”로 읽으면 틀립니다. 실제로는 **초반은 혼합, 후반은 SNU 단독 꼬리**입니다.

**지금 주장할 수 있는 것**  
merged long-horizon 주보고값으로는 **site_adjusted만** 써야 합니다. 이건 metadata에도 명시돼 있습니다.

**지금 주장하면 안 되는 것**  
merged 5년, 10년 결과를 PNU+SNU의 공평한 pooled long-term estimate로 해석하면 안 됩니다.

---

## 0-5. Diagnostics / validation / quality flags

**Source:** `stage2_validation_summary.csv` + `stage2_validation_issue_rows.csv` + `stage2_quality_flags.csv`

### 0-5a. Validation summary

|validation_group|check_name|status|n_problem_rows|
|:--|:--|:--|--:|
|cross_dataset_consistency|merged_site_consistency_horizon|pass|0|
|cross_dataset_consistency|merged_site_consistency_reverse_km|pass|0|
|cross_dataset_consistency|site_support_label_alignment|pass|0|
|internal_consistency|composition_totals|pass|0|
|internal_consistency|horizon_monotonicity|pass|0|
|internal_consistency|numbers_at_risk_monotonic|pass|0|
|internal_consistency|reporting_preference_rule|pass|0|
|internal_consistency|required_stage2_views|pass|0|
|internal_consistency|site_contribution_alignment|pass|0|
|internal_consistency|site_dominance_warning_gate|pass|0|
|raw_reconstruction|group_size_alignment_to_raw|pass|0|
|raw_reconstruction|numbers_at_risk_reconstruction_from_raw|pass|0|
|raw_reconstruction|site_admin_lookup_reconstruction|pass|0|

### 0-5b. Quality flag concentration

- completeness flags: 142개
    
- merged site-contribution flags: 48개
    
- reverse-KM incomplete CI flag: 1개
    
- dataset별로 보면 PNU completeness 54개, SNU completeness 20개, merged completeness 68개, merged site-contribution 48개입니다.
    

**숫자가 직접 말하는 것**  
구조적 검사는 전부 통과했습니다. 문제는 “계산 오류”가 아니라 “late tail이 매우 희박하다”는 해석 경고가 많다는 점입니다.

**쉬운 말로 풀면**  
파일들은 잘 만들어졌습니다. 다만 숫자들이 말해주는 임상적 메시지는 “이 뒤는 조심해서 읽어라”입니다.

**그래서 중요한 이유**  
이 경우 중요한 것은 bug 찾기가 아니라, **어디까지가 관찰이고 어디부터가 사실상 꼬리 해석인지**를 분리하는 것입니다.

**지금 주장할 수 있는 것**  
Stage 2 결과의 핵심 문제는 computational failure가 아니라 **support limitation**입니다.

**지금 주장하면 안 되는 것**  
validation pass를 scientific certainty로 읽으면 안 됩니다. pass는 “표가 맞게 만들어졌다”는 뜻이지, “long-term 결론이 강하다”는 뜻이 아닙니다.

---

# A. Scientific question reframed

진짜 질문은 “어느 모델이 더 좋으냐”가 아닙니다.  
지금 Stage 2가 답하려는 질문은 오히려 이것입니다:

> **이 코호트들에서 1년, 2년, 5년, 10년 위험을 말할 때, 어디까지가 실제 관찰에 기대는 숫자이고, 어디부터가 follow-up 구조와 tail 희소성 때문에 해석이 급격히 불안정해지는가?**

즉, 지금 단계의 핵심은 **모델 우열**이 아니라 **나중에 나올 모델 비교를 과학적으로 어디까지 믿을 수 있는지의 경계 설정**입니다.

---

# B. Bottom-line integrated interpretation

- **PNU는 “초기 분기형” 코호트**입니다. 2년 안에 전이, remission, 우측검열이 거의 다 발생해 버리고, 2년 시점 risk set이 1명뿐입니다. 그래서 PNU의 핵심 임상 질문은 장기 tail보다 **초기 1년 내 분화**입니다.
    
- **SNU는 “장기 추적 가능형” 코호트**입니다. 1–2년은 충분히 해석 가능하고, 5년도 보조적 해석이 가능하지만, 10년은 risk set 1명이라 사실상 projection 성격입니다.
    
- **merged의 후반부는 pooled cohort가 아니라 SNU tail**입니다. 5년 이후 eligible set이 100% SNU이므로, merged long-horizon 결과를 두 코호트의 평균적 자연경과로 읽으면 안 됩니다.
    
- **결과 정의의 취약성은 PNU에 집중**돼 있습니다. remission을 censoring으로 둘 때와 informative terminal outcome처럼 볼 때 2년 위험 차이가 38%p에 달하므로, PNU에서는 outcome definition이 자연경과 해석 자체를 바꿉니다.
    
- 따라서 **나중에 cure/no-cure/Bayesian cure 비교에서 가정이 덜 깨지는 곳은 SNU의 1–2년, 그 다음 merged의 1–2년**이고, **가정이 가장 흔들리는 곳은 PNU의 2년 이후**입니다.
    

---

# C. Statistical interpretation

Stage 2 수치가 가장 강하게 지지하는 것은 세 가지입니다.

첫째, **PNU와 SNU는 follow-up maturity가 전혀 다른 자료**입니다. reverse-KM 중앙 추적기간이 PNU 약 1.4–1.6년, SNU 약 4.9년이라는 점은 단순한 표본 크기 차이가 아니라 시간축 정보량의 질적 차이를 의미합니다. Betensky가 말하듯 follow-up measure는 “현재 Kaplan-Meier가 미래 complete follow-up 하에서 얼마나 안정적일지”를 보는 것이지, 진실한 장기 위험을 직접 보여주는 값이 아닙니다.

둘째, **PNU의 long flat tail은 데이터가 풍부해서 생긴 plateau가 아닙니다.** PNU overall numbers-at-risk는 54 → 28 → 1 → 0 → 0 … 이고, plot도 그 구조를 그대로 보여줍니다. 즉 2년 이후 PNU transition-only 곡선이 평평한 것은 “더 이상 전이가 없었다”는 강한 증거라기보다, “더 볼 사람이 거의 없다”는 뜻에 가깝습니다.

셋째, **merged의 late tail은 구조적으로 SNU에 종속됩니다.** site contribution이 1년부터 이미 SNU 80.3%, 2년 83.3%, 5년과 10년은 100% SNU입니다. 그래서 later-stage에서 merged와 SNU가 비슷한 long-horizon risk를 보인다면, 그것은 independent corroboration이 아니라 **같은 tail을 두 번 본 것**일 가능성이 큽니다.

숫자의 강도 측면에서 보면,

- PNU 1년은 아직 읽을 수 있지만,
    
- PNU 2년은 very sparse,
    
- PNU 3년 이후는 사실상 masking 영역입니다.  
    반면
    
- SNU 1–2년은 robust,
    
- 5년은 secondary,
    
- 10년은 descriptive/projection입니다.
    

그러므로 지금 단계에서 가장 robust한 결론은 **short-horizon 비교**입니다. long-horizon 수치는 “있다”는 사실보다 “얼마나 직접 관찰에 기대는가”를 먼저 말해야 합니다. 이런 구분이 없으면 later cure model에서 apparent cure fraction과 simple administrative immaturity를 구별하기 어려워집니다.

---

# D. Clinical interpretation

임상 맥락 정보가 이번 파일에는 제한적입니다.  
정확한 disease/syndrome 정의, `transition`의 임상적 의미, `remission` 선언 기준, site별 치료/개입 내용은 업로드된 Stage 2 파일만으로는 확정할 수 없습니다. 그래서 아래 해석은 **“transition이 주요 불량 경과, remission이 대안적 비불량 경과”라는 구조적 전제하의 조건부 해석**입니다.

그 전제에서 보면:

PNU는 **초기에 환자 경로가 빨리 갈라지는 코호트**처럼 보입니다. 2년 안에

- 일부는 전이하고,
    
- 일부는 remission으로 빠지고,
    
- 상당수는 우측검열됩니다.  
    따라서 PNU에서 임상적으로 중요한 시간창은 “장기 tail”이 아니라 **첫 1년, 길어야 2년 이내의 경로 분기**입니다.  
    이건 질병 자연경과가 정말 빠를 수도 있지만, 동시에 더 공격적 환자선별, 더 적극적 조기개입, 다른 remission 판정 관행, 혹은 짧은 행정 추적기간 때문일 수도 있습니다.
    

SNU는 **초기 전이율은 낮고, 더 오랜 기간 위험집합이 유지되는 코호트**처럼 보입니다. 즉 “시간이 지나면서 천천히 사건이 누적되는 구조”를 상대적으로 더 잘 관찰할 수 있습니다. 이런 자료는 later-stage cure/no-cure 비교에서 훨씬 유리합니다. 왜냐하면 plateau나 tail deceleration이 보여도, 그게 단순한 관찰 종료 때문인지 실제 장기안정 신호인지 상대적으로 더 분리하기 쉽기 때문입니다.

merged는 임상적으로 “두 site의 평균 환자”를 나타내는 것이 아닙니다.  
초반에는 혼합 cohort지만, 후반에는 사실상 SNU만 남습니다. 따라서 merged long-horizon 해석은 “양쪽을 섞은 보편적 자연경과”가 아니라 **“초반 혼합 + 후반 SNU 중심”**으로 읽어야 합니다.

---

# E. Model-assumption-to-clinical-reality mapping

## 1) Transition-only main view

핵심 가정은 **remission이 전이의 대안 경로가 아니라 censoring처럼 취급 가능하다**는 것입니다.  
이 가정은 **SNU에서는 덜 깨집니다.** remission이 구조적으로 0이라서 정의 민감도가 없습니다.  
반대로 **PNU에서는 더 많이 깨질 가능성**이 큽니다. remission을 포함하느냐 빼느냐에 따라 2년 위험이 38%p 차이 난다는 것은, remission이 단순 censoring이라고 보기 어렵다는 뜻입니다.

## 2) Remission-sensitive view

핵심 가정은 **remission도 임상적으로 “사건성 경로”로 취급해야 한다**는 것입니다.  
PNU에서는 이 view가 disease course의 competing pathway를 더 잘 반영할 수 있습니다.  
하지만 이것도 “주 사건 정의” 자체를 바꾸는 것이므로, later-stage no-cure/cure 비교와는 별도의 sensitivity track으로 유지해야 맞습니다.

## 3) Reverse-KM

핵심 가정은 **follow-up maturity의 기술 도구**라는 것입니다.  
이건 “얼마나 오래 따라봤는가”를 보여주지, “실제 장기 위험이 얼마인가”를 직접 말해주지 않습니다. 따라서 reverse-KM median이 길다고 cure evidence가 되는 것은 아닙니다. Stage 2 metadata도 reverse-KM을 descriptive measure로만 쓰라고 명시합니다.

## 4) Merged site_free

핵심 가정은 **pooled global administrative end date를 써도 merged completeness를 대표할 수 있다**는 것입니다.  
이번 결과는 이 가정이 덜 타당함을 보여줍니다. 특히 late horizon에서 site_free는 PNU와 SNU의 행정 종료 차이를 무시하므로, 실제보다 더 긴 잠재 추적창을 허용합니다.

## 5) Merged site_adjusted

핵심 가정은 **site별 administrative end date를 반영해야 merged maturity를 공정하게 표현한다**는 것입니다.  
이 가정이 이번 자료에서는 훨씬 더 현실적입니다. 그래서 metadata도 site_adjusted를 merged overall/sex의 preferred reporting view로 고정해 둔 것입니다.

---

# F. Why the observed pattern matters

가장 publishable한 포인트는 이것입니다:

> **이번 결과가 보여주는 핵심은 cure evidence가 아니라, “plateau처럼 보이는 것”을 만들어내는 세 요인—short follow-up, remission handling, late SNU dominance—을 명확히 분리했다는 점입니다.**

이게 왜 중요하냐면, 이후 Stage 7/8에서 cure model이 그럴듯한 tail을 만들어도, 그 tail이

- 실제 cured subgroup,
    
- 짧은 추적,
    
- remission competing path,
    
- 혹은 SNU-only late support  
    중 무엇의 반영인지 먼저 따져야 하기 때문입니다.
    

즉, 이번 Stage 2 결과 자체가 later cure modeling의 **해석 안전장치** 역할을 합니다.  
이걸 무시하면 “PNU에서 plateau가 보인다 → cure model 필요” 같은 너무 빠른 결론으로 가기 쉽습니다.

---

# G. Alternative explanations and threats to interpretation

## 강하게 지지되는 설명

1. **PNU의 짧은 follow-up이 long-tail 해석을 크게 제한한다.**
    
2. **PNU에서 remission handling이 결과 해석을 실질적으로 바꾼다.**
    
3. **merged late tail은 SNU가 사실상 독점한다.**
    

## 그럴듯하지만 아직 증명되지 않은 설명

1. **PNU와 SNU 사이에 실제 timing difference가 있다.**  
    즉, PNU는 진짜로 더 빠른 초기 전이 코호트일 수 있습니다.
    
2. **site 차이가 치료 맥락이나 care pathway 차이를 반영한다.**  
    그러나 현재 파일들에는 치료 변수나 중재 강도가 없어 분리할 수 없습니다.
    
3. **referral/severity mix가 다르다.**  
    더 심한 환자가 PNU에 많았다면 초기 사건 집중이 설명될 수 있습니다.
    

## 아직은 투기적인 설명

1. **PNU 또는 merged에 실제 cured subgroup이 존재한다.**  
    지금 단계에서는 follow-up immaturity와 remission pathway를 먼저 제거하지 않고는 강하게 말하기 어렵습니다.
    

## 무엇이 이 설명들을 가를까

- Stage 4 common-window timing analysis
    
- Stage 6 cure appropriateness screening
    
- Stage 9 remission competing-event / multi-state sensitivity
    
- site별 치료/개입/추적 프로토콜 audit
    

---

# H. Manuscript-ready discussion paragraph

In this follow-up maturity analysis, the most important finding was not evidence for cure-like heterogeneity itself but rather the clear separation of three phenomena that are often conflated in long-term time-to-event data: early timing differences between cohorts, limited late-horizon follow-up support, and sensitivity to outcome definition. PNU exhibited a markedly shorter observation window than SNU, with only one subject remaining at risk by 2 years and no observable risk set thereafter, indicating that any apparent late plateau in the PNU Kaplan-Meier curve is poorly supported by observed data. By contrast, SNU provided substantially stronger support through 2 years and usable, though still diminishing, information through 5 years. The merged cohort should not be interpreted as a stable pooled long-term trajectory, because site-adjusted analyses showed that the late eligible set was progressively dominated and then entirely determined by SNU. In addition, remission handling materially altered risk estimates in PNU but not in SNU, implying that the main transition-only estimand is conceptually more fragile in PNU. Together, these results indicate that subsequent cure-model analyses should prioritize short supported horizons, interpret merged late-tail estimates as SNU-driven, and treat apparent long-term stabilization in PNU with substantial caution unless confirmed in remission-sensitive and cure-appropriateness analyses.

---

# I. Next-step analyses

1. **Stage 4 common-window timing-difference analysis를 바로 연결**  
    특히 0–1년, 1–2년 구간에서 PNU와 SNU 차이가 follow-up artifact가 아니라 실제 timing shift인지 확인해야 합니다.
    
2. **Stage 6 cure appropriateness screening을 short supported horizon 해석과 함께 읽기**  
    screening이 positive여도 PNU 2년 이후 tail은 여전히 약합니다.
    
3. **Stage 9에서 remission을 competing event 또는 multi-state path로 재분석**  
    PNU에서는 이 단계가 사실상 필수입니다.
    
4. **site별 치료/개입/추적 프로토콜 정보 추가**  
    지금의 site effect는 치료 효과가 아니라 context bundle일 가능성이 큽니다.
    
5. **merged long-horizon 보고는 site_adjusted만 본문값으로 사용**  
    site_free는 sensitivity appendix로 남기는 것이 맞습니다.
    

---

# J. 핵심 숫자의 쉬운 해석

## 1) PNU reverse-KM 중앙 추적기간 1.37–1.57년 vs SNU 4.91년

1. 숫자가 문자 그대로 뜻하는 것  
    PNU는 대략 1년 반 정도, SNU는 거의 5년 정도의 추적 성숙도를 가졌습니다.
    
2. 쉬운 말로 바꾸면 무엇인지  
    PNU는 “짧게 본 영화”, SNU는 “훨씬 길게 본 영화”입니다.
    
3. 관찰 데이터가 이 숫자를 얼마나 직접 뒷받침하는지  
    직접적인 follow-up summary라서 비교적 직접 뒷받침됩니다.
    
4. 그래서 이 숫자를 어느 정도까지 믿어도 되는지  
    “누가 더 오래 추적됐는가”라는 주장에는 꽤 믿어도 됩니다.
    
5. 이 숫자가 임상적/과학적으로 요구하는 해석이 무엇인지  
    PNU와 SNU의 long-horizon 위험을 같은 해상도로 비교하면 안 됩니다.
    

## 2) PNU 2년 시점 risk set 1/54

1. 숫자가 문자 그대로 뜻하는 것  
    54명 중 2년을 넘겨 실제 관찰 중인 사람이 1명뿐입니다.
    
2. 쉬운 말로 바꾸면 무엇인지  
    2년 이후 PNU는 거의 아무도 남아 있지 않습니다.
    
3. 관찰 데이터가 이 숫자를 얼마나 직접 뒷받침하는지  
    numbers-at-risk 자체라 매우 직접적입니다.
    
4. 그래서 이 숫자를 어느 정도까지 믿어도 되는지  
    이건 강하게 믿어도 됩니다.
    
5. 이 숫자가 임상적/과학적으로 요구하는 해석이 무엇인지  
    PNU의 2년 이후 plateau, 5년 위험, 10년 위험은 **직접 관찰이 아니라 사실상 꼬리 해석**입니다.
    

## 3) PNU 2년 위험: transition-only 29.9% vs remission-sensitive 67.8%

1. 숫자가 문자 그대로 뜻하는 것  
    결과 정의만 바꿔도 2년 누적위험이 38.0%p 달라집니다.
    
2. 쉬운 말로 바꾸면 무엇인지  
    PNU에서는 “무엇을 사건으로 볼지”가 결과를 크게 흔듭니다.
    
3. 관찰 데이터가 이 숫자를 얼마나 직접 뒷받침하는지  
    직접 관찰 곡선에서 나온 값이지만, 2년 시점 risk set이 1이라 tail 안정성은 낮습니다.
    
4. 그래서 이 숫자를 어느 정도까지 믿어도 되는지  
    “결과 정의 민감성이 크다”는 결론은 믿어도 되지만, 정확한 2년 절대위험 수치 자체는 과신하면 안 됩니다.
    
5. 이 숫자가 임상적/과학적으로 요구하는 해석이 무엇인지  
    PNU에서는 remission을 단순 censoring으로 두는 주분석 가정이 취약합니다.
    

## 4) merged 10년 site_adjusted eligible 64 vs site_free eligible 144

1. 숫자가 문자 그대로 뜻하는 것  
    같은 merged 자료라도 구조를 어떻게 잡느냐에 따라 10년 완전추적 잠재 인원이 64명 vs 144명입니다.
    
2. 쉬운 말로 바꾸면 무엇인지  
    pooled admin end date를 쓰면 merged가 훨씬 더 성숙해 보이는 착시가 생깁니다.
    
3. 관찰 데이터가 이 숫자를 얼마나 직접 뒷받침하는지  
    직접 계산된 completeness summary입니다.
    
4. 그래서 이 숫자를 어느 정도까지 믿어도 되는지  
    site-adjusted가 실제 구조에 더 가깝다고 보는 것이 타당합니다.
    
5. 이 숫자가 임상적/과학적으로 요구하는 해석이 무엇인지  
    merged long-horizon 본문 해석은 반드시 site_adjusted 기준이어야 합니다.
    

## 5) merged late tail site dominance: 1년 80.3%, 2년 83.3%, 5년 이후 100% SNU

1. 숫자가 문자 그대로 뜻하는 것  
    merged의 eligible set이 시간이 갈수록 SNU 한쪽으로 쏠리고, 5년 이후는 완전히 SNU만 남습니다.
    
2. 쉬운 말로 바꾸면 무엇인지  
    merged 후반부는 사실상 SNU 이야기입니다.
    
3. 관찰 데이터가 이 숫자를 얼마나 직접 뒷받침하는지  
    site contribution table이 직접 보여줍니다.
    
4. 그래서 이 숫자를 어느 정도까지 믿어도 되는지  
    late-tail dominance라는 구조 해석은 강하게 믿어도 됩니다.
    
5. 이 숫자가 임상적/과학적으로 요구하는 해석이 무엇인지  
    merged 5년, 10년 결과를 PNU+SNU 평균으로 쓰면 안 되고, SNU-driven tail로 써야 합니다.
    

---

# K. 이 결과가 독자에게 실제로 요구하는 해석과 주의점

독자는 이 결과를 이렇게 써야 합니다.

- **1–2년 해석을 중심에 두십시오.**  
    특히 PNU는 1년이 핵심이고, SNU/merged는 1–2년이 핵심입니다.
    
- **PNU long tail은 과해석하지 마십시오.**  
    flat tail은 data-rich plateau가 아니라 지원 없는 연장일 가능성이 큽니다.
    
- **merged late horizon은 pooled biology가 아니라 SNU tail로 읽으십시오.**
    
- **PNU에서는 remission을 무시한 전이-only 해석을 단독 결론으로 쓰지 마십시오.**  
    remission은 단순 주변 현상이 아니라, 결과 정의를 실질적으로 바꾸는 competing path처럼 보입니다.
    
- **reverse-KM과 completeness는 관찰 성숙도의 지표이지 장기 진실의 증거가 아닙니다.**  
    follow-up measure는 KM의 안정성 해석에 도움을 주지만, 장기 위험 그 자체를 보증하지는 않습니다.
    

---

# L. 한 줄 요약

**이번 Stage 2 결과의 핵심은 “PNU는 너무 빨리 끝나고, merged의 늦은 꼬리는 사실상 SNU만 남으며, PNU에서는 remission 처리 방식이 결과를 크게 바꾼다”는 점이지, 아직 “cure가 보인다”는 결론이 아닙니다.**