중대한 계산 실패는 확인되지 않았습니다. 따라서 **Stage 8 결과는 해석 가능합니다.** 다만 해석의 무게는 코호트별로 다릅니다.

특히 **PNU는 긴 구간 cure 해석의 주근거가 아니라 민감도 해석(sensitivity)** 으로 다뤄야 하고, **SNU와 merged도 6년 이후는 본질적으로 projection** 으로 봐야 합니다. 또한 이 Stage 8의 main likelihood는 status_num == 1만 event로 두고 status_num ∈ {0,2}를 censoring으로 처리합니다. remission은 **PNU에만 존재**하고 **SNU에는 구조적으로 없습니다.** 

  

추가로 프로젝트 master spec상 **site 효과는 치료효과 그 자체가 아니라 치료 맥락, care pathway, referral, selection, follow-up 구조의 대리변수**로 해석해야 하며, supported horizon도 **PNU는 1년만 primary-supported, 2년은 sensitivity, 3년 이상은 projection; SNU·merged는 1–2년 primary, 3–5년 secondary, 6년 이상 projection** 으로 읽어야 합니다. 

---

# **0. Evidence tables**

  

## **0-1. Data structure / cohort summary**

  

**출처 파일:** bayes_stage8_dataset_registry.csv, bayes_stage8_model_registry.csv, bayes_stage8_output_audit.csv

```
Dataset-level structure
dataset   n    transition  main censoring  remission  mean age   male proportion
PNU       54   12          42              13         24.87      0.519
SNU       208  38          170             0          21.21      0.683
merged    262  50          212             13         21.96      0.649

Stage 8 model inventory
dataset   models   screening flag
PNU       8        unsupportive
SNU       8        equivocal
merged    32       equivocal

Export-layer consistency
posterior_subject_profile rows   = 10,480
expected from dataset n × models = 54×8 + 208×8 + 262×32 = 10,480

posterior_subject_yearly rows    = 104,800
expected = 10,480 × 10 horizons  = 104,800

Output audit
12 / 12 checks pass
```

- 숫자가 직접 말하는 것
    
    PNU는 표본이 54명으로 작고 transition 12건에 remission 13건이 함께 있습니다. SNU는 208명으로 더 크고 remission은 0입니다. merged는 둘을 합친 262명입니다. Stage 8은 PNU 8개, SNU 8개, merged 32개 모델을 가졌고 output audit 12개 항목이 모두 pass입니다.
    
- 쉬운 말로 풀면
    
    **파일 구조는 깨지지 않았고**, 행 수가 이상하게 불어나거나 subject가 중복 복제된 흔적은 없습니다. 다만 **PNU는 사건 수가 적고 remission 비중이 커서 원래부터 해석이 불안정해질 조건**입니다.
    
- 그래서 중요한 이유
    
    작은 코호트에서 cure fraction은 쉽게 커 보일 수 있습니다. 특히 remission을 censoring으로 두는 설계에서는 **“event가 안 난 사람”이 많아 보이면서 cure처럼 보이는 착시**가 생길 수 있습니다.
    
- 지금 주장할 수 있는 것
    
    구조적 QC 관점에서 **Stage 8 export layering, row multiplication, registry consistency는 정상**입니다.
    
- 지금 주장하면 안 되는 것
    
    구조가 정상이라고 해서 **과학적 해석까지 강하다고 말하면 안 됩니다.** 특히 PNU는 표본 자체가 약합니다.
    

---

## **0-2. Event and follow-up summary**

  

**출처 파일:** bayes_stage8_horizon_metadata.csv, bayes_stage8_ppc_summary.csv

```
Observed KM risk and Stage 8 horizon support
dataset  horizon  observed KM risk  support label in Stage 8 export     classification estimable
PNU      1        0.170             primary-supported                   yes
PNU      2        0.299             secondary                           yes
PNU      5        0.299             projection                          no
PNU      10       0.299             projection                          no

SNU      1        0.073             primary-supported                   yes
SNU      2        0.122             primary-supported                   yes
SNU      5        0.211             secondary                           yes
SNU      10       0.326             projection                          yes, with projection caution

merged   1        0.092             primary-supported                   yes
merged   2        0.153             primary-supported                   yes
merged   5        0.239             secondary                           yes
merged   10       0.351             projection                          yes, with projection caution
```

- 숫자가 직접 말하는 것
    
    PNU는 observed KM risk가 2년 0.299 이후 10년까지 그대로 0.299입니다. 반면 SNU는 1년 0.073에서 10년 0.326까지, merged는 1년 0.092에서 10년 0.351까지 계속 올라갑니다.
    
- 쉬운 말로 풀면
    
    **PNU는 2년 이후 더 이상 실제 관찰된 event 증가가 거의 안 보입니다.** 하지만 그게 “진짜 cure”인지, 아니면 “follow-up이 짧아서 더 못 본 것”인지는 별개입니다. SNU와 merged는 더 오래 관찰되어 뒤쪽 tail 정보가 더 있습니다.
    
- 그래서 중요한 이유
    
    같은 plateau라도 의미가 다릅니다. PNU plateau는 **정보가 없어서 멈춘 것**일 수 있고, SNU plateau는 **좀 더 관찰된 뒤의 완만화**일 수 있습니다.
    
- 지금 주장할 수 있는 것
    
    **PNU의 3년 이후 threshold classification은 계산적으로도 비추정(not estimable)** 입니다.
    
- 지금 주장하면 안 되는 것
    
    PNU의 2년 이후 평평한 KM만 보고 **“장기 cure가 있다”**고 단정하면 안 됩니다.
    

---

## **0-3. Follow-up support by horizon**

  

**출처 파일:** bayes_stage8_horizon_metadata.csv, Integrated_Modeling_Master_Specification_English.md

```
Interpretation support table
dataset   1 year                  2 years                    5 years                       10 years
PNU       directly supported      limited support            mostly extrapolated            mostly extrapolated
SNU       directly supported      directly supported         partly model-dependent         mostly extrapolated
merged    directly supported      directly supported         partly model-dependent         mostly extrapolated
```

여기서 한 가지 **minor caution**이 있습니다.

Stage 8 export에서는 **PNU 2년이 secondary로 라벨링**되어 있는데, 프로젝트 master spec은 **PNU 2년을 sensitivity로 해석하라고 규정**합니다. 이는 **숫자를 바꾸는 계산 오류는 아니고**, **지원수준(label) 표현의 불일치**입니다. 해석 단계에서는 master spec을 우선해 **PNU 2년을 “제한적 지원/민감도 구간”**으로 읽는 것이 맞습니다. 

- 숫자가 직접 말하는 것
    
    PNU는 1년만 진짜 주해석 구간이고, 2년도 약합니다. SNU와 merged는 1–2년은 비교적 안정적이고, 5년부터는 모델 의존성이 커집니다.
    
- 쉬운 말로 풀면
    
    **100명 중 몇 명이 생길지**를 말할 때, PNU는 1년만 꽤 믿을 수 있고 그 뒤는 많이 추정입니다. SNU와 merged는 2년까지는 꽤 믿을 만하지만 10년은 예측 비중이 큽니다.
    
- 그래서 중요한 이유
    
    long-horizon 숫자를 똑같은 무게로 비교하면 안 됩니다. **1년 숫자와 10년 숫자는 증거의 질이 다릅니다.**
    
- 지금 주장할 수 있는 것
    
    SNU·merged의 1–2년 비교는 주해석으로 사용할 수 있습니다.
    
- 지금 주장하면 안 되는 것
    
    PNU 5년, 10년 값을 실제 관찰치처럼 말하면 안 됩니다.
    

---

## **0-4. Model fit, convergence, and structural QC**

  

**출처 파일:** bayes_stage8_model_registry.csv, bayes_stage8_output_audit.csv

```
Dataset-level diagnostics
dataset  models  fit_ok  admissible  reused  max div  max treedepth  max Rhat  min bulk ESS  min tail ESS  LOO warn  WAIC warn  Pareto bad
PNU      8       8       8           8       0        0              1.004     1758          1839          3         7          3
SNU      8       8       8           8       0        0              1.006     1638          1844          0         0          0
merged   32      32      32          32      0        0              1.005     1259           861          3         2          3
```

추가 QC 포인트:

- fit_status == ok는 48/48
    
- admissible_flag == TRUE는 48/48
    
- divergences == 0, treedepth_exceeded == 0는 전 모델
    
- ppc_gross_contradiction_flag == FALSE도 전 모델
    
- 다만 **all models were reused from prior Stage 8 exports** 입니다. 이번 업로드는 재적합이 아니라 **기존 export 재사용** 상태입니다. 이는 Stage 8 코드가 허용하는 경로라서 fatal issue는 아닙니다.
    
- 숫자가 직접 말하는 것
    
    샘플러가 터지거나 발산한 모델은 없습니다. Rhat도 모두 1.01 미만입니다.
    
- 쉬운 말로 풀면
    
    **컴퓨터가 모델을 못 맞춘 건 아닙니다.** 숫자가 튀어서 “수학적으로 실패한 결과”는 아닙니다.
    
- 그래서 중요한 이유
    
    해석을 막는 종류의 계산 실패는 없으므로 downstream interpretation 자체는 가능합니다.
    
- 지금 주장할 수 있는 것
    
    **Fatal computational failure는 없다**고 말할 수 있습니다.
    
- 지금 주장하면 안 되는 것
    
    LOO/WAIC 경고가 일부 없다고 해서 **모델 가정이 참**이라고 말할 수는 없습니다. 그건 가정의 적합성 문제가 아니라 fitting stability 문제와 다른 층위입니다.
    

---

## **0-5. Parameter / shape summary**

  

**출처 파일:** bayes_stage8_posterior_cohort_yearly.csv, bayes_stage8_model_registry.csv

  

아래 표는 raw coefficient가 아니라, **임상적으로 해석 가능한 요약 파라미터**입니다.

cohort_mean_cure_fraction_mean은 cohort 평균 cure-like fraction 요약이고,

hazard_10y / hazard_1y는 10년 hazard가 1년 hazard의 몇 배인지 보여줍니다.

```
Family-level summary
dataset  family        cure fraction range         median hazard(10y)/hazard(1y)
PNU      exponential   0.454 to 0.489              1.000
PNU      loglogistic   0.416 to 0.418              0.183
PNU      lognormal     0.427 to 0.437              0.216
PNU      weibull       0.390 to 0.422              0.677

SNU      exponential   0.383 to 0.401              1.000
SNU      loglogistic   0.226 to 0.230              0.372
SNU      lognormal     0.214 to 0.225              0.289
SNU      weibull       0.237 to 0.238              0.544

merged   exponential   0.402 to 0.562              1.000
merged   loglogistic   0.194 to 0.234              0.323
merged   lognormal     0.192 to 0.237              0.263
merged   weibull       0.208 to 0.286              0.511
```

- 숫자가 직접 말하는 것
    
    exponential family는 모든 dataset에서 hazard ratio가 정확히 1.000입니다. 즉 susceptible hazard를 시간에 따라 **전혀 줄어들지 않는 상수**로 가정합니다. 반면 lognormal/loglogistic는 10년 hazard가 1년 hazard의 약 18%–37% 수준으로 떨어집니다.
    
- 쉬운 말로 풀면
    
    exponential은 **“처음 위험도, 나중 위험도 거의 같다”**는 모델이고, lognormal/loglogistic는 **“초반 위험이 높고 시간이 지나면 위험이 확 줄어든다”**는 모델입니다.
    
- 그래서 중요한 이유
    
    이 데이터의 observed pattern은 특히 PNU에서 **초반 event가 몰리고 뒤가 평평**합니다. 이런 자료에 constant hazard는 개념적으로 덜 맞고, 그 보상으로 cure fraction을 과대하게 밀어 올릴 가능성이 있습니다.
    
- 지금 주장할 수 있는 것
    
    **개념적으로 덜 위배되어 보이는 가정은 exponential보다 declining-hazard family**입니다.
    
- 지금 주장하면 안 되는 것
    
    이 표만 보고 특정 family가 “정답”이라고 단정하면 안 됩니다. 우리는 ranking이 아니라 **가정-현실 매핑**을 보고 있습니다.
    

---

## **0-6. Prediction or time-specific risk table**

  

**출처 파일:** bayes_stage8_posterior_cohort_yearly.csv, bayes_stage8_ppc_summary.csv

  

### **전체 모델 범위 요약**

```
Observed KM risk versus Bayesian risk range across ALL Stage 8 models
dataset  horizon  observed KM   Bayesian risk min  median  max
PNU      1        0.170         0.204              0.208   0.214
PNU      2        0.299         0.274              0.281   0.310
PNU      5        0.299         0.359              0.378   0.441
PNU      10       0.299         0.412              0.442   0.506

SNU      1        0.073         0.060              0.083   0.089
SNU      2        0.122         0.111              0.129   0.137
SNU      5        0.211         0.219              0.222   0.229
SNU      10       0.326         0.296              0.317   0.355

merged   1        0.092         0.087              0.107   0.116
merged   2        0.153         0.151              0.162   0.170
merged   5        0.239         0.249              0.262   0.287
merged   10       0.351         0.323              0.356   0.404
```

### **10년 complete comparison block: PNU**

```
PNU, 10-year, all 8 models
model_id   family        formula        Bayes risk   cure fraction   observed KM   abs diff
PNU-L0-E   exponential   base           0.482        0.489           0.299         0.183
PNU-L0-LL  loglogistic   base           0.421        0.416           0.299         0.122
PNU-L0-LN  lognormal     base           0.412        0.437           0.299         0.114
PNU-L0-W   weibull       base           0.463        0.422           0.299         0.164
PNU-L1-E   exponential   interaction    0.506        0.454           0.299         0.207
PNU-L1-LL  loglogistic   interaction    0.418        0.418           0.299         0.119
PNU-L1-LN  lognormal     interaction    0.414        0.427           0.299         0.115
PNU-L1-W   weibull       interaction    0.478        0.390           0.299         0.179
```

### **10년 complete comparison block: SNU**

```
SNU, 10-year, all 8 models
model_id   family        formula        Bayes risk   cure fraction   observed KM   abs diff
SNU-L0-E   exponential   base           0.353        0.401           0.326         0.027
SNU-L0-LL  loglogistic   base           0.313        0.226           0.326         0.013
SNU-L0-LN  lognormal     base           0.297        0.225           0.326         0.030
SNU-L0-W   weibull       base           0.325        0.238           0.326         0.001
SNU-L1-E   exponential   interaction    0.355        0.383           0.326         0.028
SNU-L1-LL  loglogistic   interaction    0.310        0.230           0.326         0.017
SNU-L1-LN  lognormal     interaction    0.296        0.214           0.326         0.031
SNU-L1-W   weibull       interaction    0.321        0.237           0.326         0.005
```

### **10년 merged branch summary**

```
merged, 10-year
branch                           Bayes risk min  median  max   cure fraction min  median  max
main structural                  0.327           0.355   0.402 0.200              0.231   0.550
incidence-site supplementary     0.323           0.356   0.404 0.192              0.229   0.562
```

- 숫자가 직접 말하는 것
    
    PNU는 2년 이후 observed KM이 0.299에서 멈췄는데, Bayesian cure risk는 10년에 0.412–0.506으로 계속 올라갑니다. 반면 SNU는 10년 observed 0.326에 대해 Bayesian 0.296–0.355로 꽤 가깝습니다. merged도 10년 observed 0.351, Bayesian 0.323–0.404로 상대적으로 가깝습니다.
    
- 쉬운 말로 풀면
    
    **PNU는 모델이 뒤를 많이 상상하고 있고**, **SNU/merged는 모델이 실제 관찰 패턴을 꽤 따라가고 있습니다.**
    
- 그래서 중요한 이유
    
    PNU의 큰 cure-like signal은 “데이터가 분명히 말한 것”이라기보다 **짧은 follow-up 뒤를 모델이 채운 결과**일 가능성이 큽니다.
    
- 지금 주장할 수 있는 것
    
    **SNU와 merged에서는 Bayesian cure block이 observed data와 큰 충돌 없이 작동**합니다.
    
    **PNU에서는 2년 이후 숫자가 projection 중심**입니다.
    
- 지금 주장하면 안 되는 것
    
    PNU 10년 risk 0.41–0.51을 **실제 관찰 기반 장기위험**처럼 말하면 안 됩니다.
    

---

## **0-7. Bayesian cure versus no-cure comparison**

  

**출처 파일:** bayes_stage8_posterior_delta_vs_nocure.csv

  

delta_mean = no-cure value - Bayesian cure value 이므로

양수면 Bayesian cure가 위험 또는 burden을 낮춘 것이고, 음수면 오히려 Bayesian cure가 더 높인 것입니다.

```
Median delta across all Stage 8 models
metric = meanRisk
dataset  1y      2y      5y      10y
PNU      +0.011  +0.043  +0.092  +0.165
SNU      -0.007  -0.006  -0.001  +0.012
merged   -0.005  -0.002  +0.004  +0.022

metric = net benefit
dataset  1y      2y      5y      10y
PNU      +0.001  +0.007    NA      NA/unstable
SNU      +0.002  +0.002  +0.003  +0.001
merged   +0.001  +0.002  +0.001  +0.000
```

- 숫자가 직접 말하는 것
    
    SNU와 merged의 supported horizon에서 meanRisk delta 중앙값은 거의 0에 가깝고, net benefit 개선도 매우 작습니다.
    
- 쉬운 말로 풀면
    
    **초기 1–2년, 심지어 SNU 5년까지도 “Bayesian cure를 썼더니 임상 결정이 확 달라진다”는 신호는 약합니다.**
    
- 그래서 중요한 이유
    
    이 프로젝트의 핵심은 model ranking이 아니라 **false-positive burden과 clinical usefulness가 의미 있게 바뀌는지**인데, Stage 8 단독 결과만 보면 그 변화는 early supported horizons에서 크지 않습니다.
    
- 지금 주장할 수 있는 것
    
    **Bayesian cure block의 주된 차이는 early supported horizon보다는 late tail에서 더 크게 나타납니다.**
    
- 지금 주장하면 안 되는 것
    
    “Bayesian cure가 no-cure보다 임상적으로 확실히 우월하다”는 결론은 이 Stage 8 결과만으로는 못 냅니다.
    

---

## **0-8. Plot/table consistency**

  

**출처 파일:** bayes_stage8_diagnostic_plots.pdf

  

Diagnostic PDF의 요약 그림은 표와 잘 맞습니다.

- **49쪽 risk trajectories**: PNU가 가장 높고 불확실성도 넓으며, SNU는 더 낮고 부드럽고, merged는 중간입니다.
    
- **50쪽 hazard trajectories**: exponential 계열은 수평선처럼 flat hazard를 보이고, 나머지 계열은 시간에 따라 하강합니다.
    
- **51쪽 net benefit**: 특히 SNU와 merged에서는 threshold에 따라 큰 모델 간 벌어짐보다 완만한 차이가 많습니다.
    
- **52쪽 PPC**: PNU는 observed KM marker가 2년 이후 거의 평평한데 posterior line은 더 올라가고, SNU·merged는 observed와 posterior가 상대적으로 더 가깝습니다. 
    
- 숫자가 직접 말하는 것
    
    표에서 본 패턴이 그림에서도 그대로 보입니다.
    
- 쉬운 말로 풀면
    
    **CSV와 PDF가 서로 싸우지 않습니다.** 숫자표와 그림이 같은 얘기를 합니다.
    
- 그래서 중요한 이유
    
    표-그림 불일치가 있으면 export나 plotting bug를 의심해야 하는데, 여기서는 그 징후가 없습니다.
    
- 지금 주장할 수 있는 것
    
    plot/table consistency는 양호합니다.
    
- 지금 주장하면 안 되는 것
    
    그림이 예쁘게 보인다고 long-horizon 숫자의 지지수준이 높아지는 것은 아닙니다.
    

---

# **A. Scientific question reframed**

  

이 결과의 진짜 질문은

**“어떤 family가 AIC나 WAIC에서 이겼는가?”**가 아니라,

**“transition-only outcome을 기준으로 했을 때, 이 데이터가 실제로 cure-like 하위집단을 지지하는지, 아니면 조기 timing difference, follow-up 부족, PNU의 remission 구조가 만들어낸 겉모양인지”** 입니다.

  

조금 더 쉽게 말하면,

**“장기적으로 안 나빠지는 사람이 정말 있는 건지, 아니면 아직 덜 본 건지를 구분할 수 있는가?”**가 핵심입니다.

---

# **B. Bottom-line integrated interpretation**

- **계산 기반은 멀쩡합니다.** 48개 모델 모두 fit status가 ok이고 admissible이며 발산도 없었습니다. 그래서 “컴퓨터가 잘못 맞춘 결과”는 아닙니다.
    
- **하지만 PNU의 장기 cure 해석은 약합니다.** PNU는 54명, transition 12건, remission 13건이고, 2년 이후 observed KM risk가 0.299에서 더 올라가지 않는데 Bayesian 모델은 10년에 0.412–0.506까지 계속 올립니다. 이건 “관찰이 말한 cure”보다 “짧은 follow-up 뒤를 모델이 메운 결과”에 가깝습니다.
    
- **SNU와 merged는 modest한 cure-like heterogeneity는 가능하지만, dramatic한 조기 임상 이득 증거는 약합니다.** SNU와 merged는 1–2년, 그리고 SNU 5년까지 observed KM과 Bayesian risk가 꽤 잘 맞습니다. 다만 no-cure 대비 false-positive burden이나 net benefit의 early supported horizon 개선은 크지 않습니다.
    
- **가장 덜 위배되어 보이는 latency 가정은 declining-hazard 계열입니다.** exponential처럼 hazard를 평평하게 두는 가정은 이 자료의 “초기 집중, 이후 완만화” 패턴과 덜 맞고, 그 보상으로 cure fraction을 더 크게 밀어 올리는 경향이 있습니다.
    
- **따라서 이 Stage 8 결과가 가장 강하게 말하는 것은 ‘공통의 강한 cure 신호’가 아니라 ‘코호트별로 다른 timing and follow-up structure’입니다.** PNU에서 보이는 큰 cure-like signal은 sensitivity 수준이고, SNU/merged에서만 조심스러운 cure-like heterogeneity 가능성을 논하는 편이 더 방어적입니다.
    

---

# **C. Statistical interpretation**

  

첫째, **computationally clean** 합니다. Rhat가 모두 1.01 미만이고 divergences와 treedepth 문제도 없습니다. 이는 posterior summary 자체를 읽어도 된다는 뜻입니다.

  

둘째, **posterior predictive fit은 코호트별로 다릅니다.**

PNU는 1–2년에서는 observed와 아주 심하게 안 맞지는 않지만, 5년과 10년으로 갈수록 absolute difference가 커집니다. 반대로 SNU와 merged는 1–2년뿐 아니라 5년, 심지어 10년에서도 모델 범위가 observed KM과 크게 떨어지지 않습니다.

  

셋째, **family 차이는 주로 tail에서 벌어집니다.**

PNU 10년에서 exponential은 0.482–0.506, lognormal/loglogistic는 0.412–0.421입니다. SNU 10년에서는 exponential 0.353–0.355, Weibull 0.321–0.325, lognormal 0.296–0.297입니다. 즉 **supported early horizon에서는 큰 차이가 아니지만, unsupported tail에서 family choice가 숫자를 많이 바꿉니다.**

  

넷째, **Bayesian cure와 no-cure의 early-horizon 차이는 작습니다.**

SNU와 merged에서 1–2년 meanRisk delta 중앙값이 거의 0이고, net benefit도 +0.001~+0.003 수준입니다.

이건 **cure modeling이 short-horizon clinical decision을 완전히 뒤집는 상황은 아니다**는 뜻입니다.

---

# **D. Clinical interpretation**

  

이 결과를 질병 자연경과 관점에서 보면, **transition 위험은 일정한 속도로 계속 가는 병이라기보다 초반에 몰리고 시간이 지나며 완화되는 구조**로 읽는 편이 더 자연스럽습니다. 즉 “계속 같은 속도로 나빠진다”보다는 “초기 취약 집단이 먼저 사건을 겪고, 남아 있는 집단은 상대적으로 안정적”인 모양입니다.

  

그런데 이걸 곧바로 **“cured subgroup가 분명하다”**고 읽으면 안 됩니다. 왜냐하면:

1. **PNU는 remission이 실제로 존재하고, main Stage 8에서는 remission을 censoring으로 처리**합니다.
    
    그래서 transition을 안 겪은 사람이 많아 보이면, 그 일부는 “진짜 cure”가 아니라 **remission 또는 competing process가 censoring으로 들어간 결과**일 수 있습니다. 
    
2. **PNU follow-up은 짧고 late classification도 not estimable** 입니다.
    
    즉 뒤쪽 숫자는 자연경과를 직접 본 결과가 아니라 projection입니다.
    
3. **site는 treatment 자체가 아니라 broader context의 proxy** 입니다.
    
    따라서 PNU–SNU 차이는 생물학적 차이일 수도 있지만, care pathway, referral, selection, follow-up maturity 차이일 수도 있습니다. 
    

  

결국 임상적으로는 이렇게 읽는 것이 가장 안전합니다.

**“일부 코호트에서 장기적으로 상대적으로 안정적인 하위집단이 있을 가능성은 있다. 그러나 그 신호의 크기와 확실성은 코호트에 따라 다르고, PNU에서는 remission 구조와 추적 제한 때문에 true cure로 보기 어렵다.”**

---

# **E. Model-assumption-to-clinical-reality mapping**

  

## **1) Exponential cure model**

- 핵심 가정: susceptible hazard가 시간에 따라 일정
    
- 여기서의 임상적 타당성: **가장 약함**
    
- 왜냐하면: observed pattern은 초기에 event가 몰리고 뒤가 완만합니다
    
- 결과 패턴: cure fraction을 더 크게 잡는 경향, 특히 PNU·merged에서 두드러짐
    

  

쉽게 말하면,

**“초기 위험 집중”을 설명하지 못하니 그 부족분을 “cure 비율이 높다”로 메우는 모델**처럼 행동합니다.

  

## **2) Weibull cure model**

- 핵심 가정: monotonic hazard
    
- 여기서의 임상적 타당성: **중간**
    
- 왜냐하면: hazard가 내려가는 방향은 허용하지만, 더 복잡한 hump/flattening은 덜 유연
    
- 결과 패턴: exponential보다 현실적이지만 tail family dependence는 남음
    

  

쉽게 말하면,

**“위험이 점점 줄어든다”는 정도는 담지만, 실제 패턴이 더 휘어져 있으면 부족할 수 있습니다.**

  

## **3) Lognormal / loglogistic cure model**

- 핵심 가정: 초기 위험 상승 후 감소 또는 보다 강한 비대칭 tail 허용
    
- 여기서의 임상적 타당성: **가장 덜 위배되어 보임**
    
- 왜냐하면: 초반 사건 집중과 뒤 완만화를 더 잘 수용
    
- 결과 패턴: 특히 PNU와 SNU에서 exponential보다 덜 과격한 long-tail risk
    

  

쉽게 말하면,

**“처음 위험은 높고 나중에 진정된다”는 질병 경과에 더 잘 맞는 모형군**입니다.

  

## **4) merged site-free model**

- 핵심 가정: site 차이를 구조적으로 따로 설명하지 않음
    
- 임상적 타당성: **제한적**
    
- 왜냐하면: PNU와 SNU는 follow-up과 event timing이 다릅니다
    
- 결과 패턴: 서로 다른 코호트를 평균내는 위험
    

  

쉽게 말하면,

**“성격 다른 두 집단을 한 덩어리로 본다”**는 점에서 해석이 흐려질 수 있습니다.

  

## **5) merged site-adjusted model**

- 핵심 가정: site 차이를 incidence 또는 latency에서 일부 구조적으로 설명
    
- 임상적 타당성: **site-free보다 높음**
    
- 왜냐하면: 실제로 cohorts가 다르기 때문
    
- 단, site는 clean treatment effect가 아니므로 과해석 금지
    

  

쉽게 말하면,

**“두 병원을 같은 집단으로 보지 않는다”**는 점에서 더 현실적이지만, 그 차이를 곧 치료효과로 읽으면 안 됩니다.

---

# **F. Why the observed pattern matters**

  

이 결과의 가장 publishable한 점은

**“cure-like signal이 보인다고 해서 모두 같은 의미가 아니다”**는 점입니다.

  

PNU에서는 큰 cure fraction이 나와도 그 신호는:

- 작은 n,
    
- remission 존재,
    
- main analysis에서 remission censoring,
    
- 2년 이후 projection 중심,
    
- Stage 6 unsupportive
    
    라는 조건 위에서 나온 것입니다.
    

  

반면 SNU와 merged에서는:

- 더 긴 follow-up,
    
- observed와 posterior의 더 좋은 정합,
    
- late horizon까지 classification 자체는 가능
    
    이라는 점 때문에 **훨씬 더 방어적인 cure-like heterogeneity 논의가 가능합니다.**
    

  

즉, **같은 “plateau”라도 PNU plateau와 SNU plateau는 과학적 의미가 다릅니다.**

이 차이를 보여주는 것이 이 분석의 핵심 가치입니다.

---

# **G. Alternative explanations and threats to interpretation**

  

## **강하게 지지되는 설명**

1. **PNU의 follow-up immaturity**
    
2. **PNU에서 remission을 censoring으로 둔 구조적 영향**
    
3. **PNU–SNU 간 timing/context 차이**
    

  

이 셋은 현재 데이터가 가장 강하게 지지합니다.

  

## **그럴듯하지만 아직 확정 못 하는 설명**

1. **SNU와 merged에서 modest한 cure-like subgroup 존재**
    
2. **site-specific care pathway or selection 차이가 장기 tail을 바꿈**
    

  

이건 plausible합니다. 하지만 definitive하진 않습니다.

  

## **아직 speculative한 설명**

1. **PNU와 SNU의 진짜 생물학적 cure mechanism 차이**
    
2. **site effect를 실제 treatment effect로 해석하는 것**
    

  

현재 자료로는 이 수준까지는 못 갑니다.

  

무엇이 이 설명들을 가를까?

- remission을 competing risk로 둔 Stage 9 결과
    
- common-window 1년, 2년 비교
    
- site-specific treatment/process variable
    
- 더 긴 PNU follow-up
    

---

# **H. Manuscript-ready discussion paragraph**

  

본 Bayesian mixture-cure analysis는 계산 안정성 측면에서는 양호하였으며, 모든 모델이 수렴하고 posterior predictive contradiction도 관찰되지 않았다. 그러나 이러한 계산적 안정성이 곧바로 강한 생물학적 cure evidence를 의미하지는 않았다. 특히 PNU에서는 표본 수가 작고 remission이 존재함에도 본 단계의 주분석에서 remission을 censoring으로 처리하였으며, 2년 이후는 사실상 projection에 의존하므로 장기 cure fraction 추정은 민감도 수준의 해석에 머물러야 한다. 반면 SNU와 merged에서는 1–2년의 주지지 구간에서 Bayesian risk가 observed KM과 대체로 양호하게 일치하였고, 5년까지도 큰 모순 없이 작동하였다. 다만 no-cure 모델 대비 false-positive burden 및 net benefit의 개선은 supported horizon에서 크지 않아, Bayesian cure block의 주된 차별점은 short-horizon 임상의사결정보다는 late-tail 구조를 어떻게 해석하느냐에 더 가까웠다. 전반적으로 본 결과는 공통된 강한 cure subgroup의 존재보다는, 코호트별 timing difference와 follow-up structure의 이질성이 관찰된 plateau와 장기 위험 추정의 차이를 상당 부분 설명할 가능성을 시사한다.

---

# **I. Next-step analyses**

- **Stage 9 remission sensitivity를 최우선으로** 하십시오.
    
    PNU에서 remission을 censoring으로 둔 현재 backbone이 cure-like signal을 얼마나 부풀렸는지 확인해야 합니다.
    
- **PNU–SNU common-window 1년, 2년 direct comparison을 강화** 하십시오.
    
    장기 tail이 아니라 초반부터 timing difference가 있는지 분리해야 합니다.
    
- **supported horizon 중심의 decision analysis를 재정리** 하십시오.
    
    특히 SNU·merged의 1–2년, 그리고 5년을 중심으로 Bayesian cure vs no-cure의 false-positive burden과 net benefit를 다시 정리하는 것이 좋습니다.
    
- **family별 tail behavior를 본문에서는 shape language로, 부록에서는 수치로** 제시하십시오.
    
    exponential의 flat hazard 가정이 왜 덜 그럴듯한지 독자가 즉시 이해할 수 있게 해야 합니다.
    
- **site-specific process/treatment 정보가 있다면 반드시 추가** 하십시오.
    
    지금은 site가 치료효과인지 care pathway인지 분리할 수 없습니다.
    

---

# **J. 핵심 숫자의 쉬운 해석**

  

## **1) PNU: 54명 중 transition 12건, remission 13건**

1. 숫자가 문자 그대로 뜻하는 것
    
    PNU에서는 transition event 수보다 remission 수가 오히려 더 많습니다.
    
2. 쉬운 말로 바꾸면 무엇인지
    
    **안 나빠진 사람**이 많아 보이는데, 그 안에는 **진짜 장기 안정군**도 있을 수 있고 **remission으로 다른 경로를 탄 사람**도 섞여 있을 수 있습니다.
    
3. 관찰 데이터가 이 숫자를 얼마나 직접 뒷받침하는지
    
    이건 직접 관찰입니다.
    
4. 그래서 이 숫자를 어느 정도까지 믿어도 되는지
    
    숫자 자체는 믿어도 됩니다.
    
5. 이 숫자가 임상적/과학적으로 요구하는 해석이 무엇인지
    
    PNU의 plateau를 곧장 cure로 읽지 말고, **remission 구조를 먼저 분리**해야 합니다.
    

  

## **2) PNU observed KM risk: 2년 0.299, 10년도 0.299**

1. 숫자가 문자 그대로 뜻하는 것
    
    2년까지 약 29.9%가 transition을 겪었고, 그 이후 observed KM 상 추가 상승이 없습니다.
    
2. 쉬운 말로 바꾸면 무엇인지
    
    **100명 중 약 30명은 2년 안에 사건을 겪었고, 그 뒤에는 더 본 사건이 거의 없다**는 뜻입니다.
    
3. 관찰 데이터가 이 숫자를 얼마나 직접 뒷받침하는지
    
    2년까지는 직접적, 그 이후 plateau는 관찰되지만 뒤 risk set 정보는 약합니다.
    
4. 그래서 이 숫자를 어느 정도까지 믿어도 되는지
    
    **2년 값은 제한적으로 믿을 수 있지만**, 2년 이후 plateau의 의미는 조심해야 합니다.
    
5. 이 숫자가 임상적/과학적으로 요구하는 해석이 무엇인지
    
    **“2년 이후 완치”가 아니라 “2년 이후 정보가 빈약하다”**가 우선 해석입니다.
    

  

## **3) PNU Bayesian 10년 risk: 0.412–0.506**

1. 숫자가 문자 그대로 뜻하는 것
    
    모델에 따라 10년 transition risk가 41.2%에서 50.6%로 추정됩니다.
    
2. 쉬운 말로 바꾸면 무엇인지
    
    **100명 중 41명에서 51명 정도**가 장기적으로 transition할 수 있다고 모델이 상상한 것입니다.
    
3. 관찰 데이터가 이 숫자를 얼마나 직접 뒷받침하는지
    
    약합니다. 거의 projection입니다.
    
4. 그래서 이 숫자를 어느 정도까지 믿어도 되는지
    
    **직접 관찰치처럼 믿으면 안 됩니다.**
    
5. 이 숫자가 임상적/과학적으로 요구하는 해석이 무엇인지
    
    PNU의 장기위험은 **데이터가 아니라 모델 가정과 prior의 영향이 크다**는 뜻입니다.
    

  

## **4) SNU observed KM 10년 risk 0.326, Bayesian 0.296–0.355**

1. 숫자가 문자 그대로 뜻하는 것
    
    observed 32.6%, Bayesian range 29.6%–35.5%입니다.
    
2. 쉬운 말로 바꾸면 무엇인지
    
    **100명 중 약 30명대 초중반**이 장기적으로 transition할 것이라는 점에서 observed와 model이 대체로 비슷합니다.
    
3. 관찰 데이터가 이 숫자를 얼마나 직접 뒷받침하는지
    
    PNU보다 훨씬 낫지만, 10년은 여전히 projection caution이 있습니다.
    
4. 그래서 이 숫자를 어느 정도까지 믿어도 되는지
    
    **방향성과 대략적 규모는 믿을 수 있지만 소수점 단위 비교는 금물**입니다.
    
5. 이 숫자가 임상적/과학적으로 요구하는 해석이 무엇인지
    
    SNU에서는 **modest cure-like heterogeneity 가능성**을 방어적으로 논할 수 있습니다.
    

  

## **5) merged Bayesian–no-cure meanRisk delta: 1년 -0.005, 2년 -0.002, 5년 +0.004**

1. 숫자가 문자 그대로 뜻하는 것
    
    no-cure와 Bayesian cure의 위험 차이 중앙값이 early horizon에서 거의 0입니다.
    
2. 쉬운 말로 바꾸면 무엇인지
    
    **Bayesian cure를 써도 초반 위험 예측은 거의 안 바뀝니다.**
    
3. 관찰 데이터가 이 숫자를 얼마나 직접 뒷받침하는지
    
    supported horizon에서는 꽤 직접적입니다.
    
4. 그래서 이 숫자를 어느 정도까지 믿어도 되는지
    
    early horizon에서 “차이가 작다”는 결론은 상당히 믿을 만합니다.
    
5. 이 숫자가 임상적/과학적으로 요구하는 해석이 무엇인지
    
    Stage 8의 실질적 차이는 **초반 임상결정이 아니라 late-tail 해석**에 더 있습니다.
    

  

## **6) 모든 모델에서 prior tail warning이 존재**

1. 숫자가 문자 그대로 뜻하는 것
    
    48개 모델 모두에서 prior_tail_warning_flag_any = TRUE입니다.
    
2. 쉬운 말로 바꾸면 무엇인지
    
    prior가 **매우 긴 susceptible tail**을 허용한다는 뜻입니다.
    
3. 관찰 데이터가 이 숫자를 얼마나 직접 뒷받침하는지
    
    이건 관찰 데이터가 아니라 모델 설정에서 나온 경고입니다.
    
4. 그래서 이 숫자를 어느 정도까지 믿어도 되는지
    
    long-horizon 숫자를 볼 때 항상 같이 기억해야 하는 경고입니다.
    
5. 이 숫자가 임상적/과학적으로 요구하는 해석이 무엇인지
    
    **장기 추정치는 prior-sensitive하다**는 점을 독자에게 분명히 알려야 합니다.
    

---

# **K. 이 결과가 독자에게 실제로 요구하는 해석과 주의점**

- **1–2년을 중심으로 읽으십시오.**
    
    특히 PNU는 1년, SNU·merged는 1–2년이 주해석입니다.
    
- **PNU 5년·10년 숫자를 실측처럼 비교하지 마십시오.**
    
    그 구간은 projection이고 threshold classification도 비추정입니다.
    
- **site 차이를 치료효과로 읽지 마십시오.**
    
    지금 단계에서는 treatment context, care pathway, selection, follow-up의 proxy입니다. 
    
- **PNU plateau를 cure 증거로 과대해석하지 마십시오.**
    
    remission-as-censoring과 짧은 follow-up이 더 직접적인 설명입니다.
    
- **SNU와 merged에서만 조심스러운 cure-like heterogeneity를 논하십시오.**
    
    그것도 “강한 확증”이 아니라 “plausible signal” 수준이 적절합니다.
    
- **family 간 소수점 승부를 피하고, hazard shape의 plausibility를 보십시오.**
    
    exponential의 flat hazard는 이 자료의 초반 집중-후반 완화 패턴과 덜 맞습니다.
    

---

# **L. 한 줄 요약**

  

**이번 Stage 8 결과는 계산적으로는 깨끗하지만, 과학적으로는 “PNU의 큰 cure 신호”보다 “SNU·merged에서는 modest한 cure 가능성, PNU에서는 follow-up과 remission 구조의 영향”을 더 강하게 말해 줍니다.**