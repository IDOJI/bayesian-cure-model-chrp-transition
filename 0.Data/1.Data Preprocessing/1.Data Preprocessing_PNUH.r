# =================== 패키지 ===================
library(dplyr)

# =================== 데이터 로드 ===================
path_data <- '/Volumes/ObsidianVault/🟩Obsidian/☔️Papers_Writing(논문 쓰기)/✅Currently working/⬛A Bayesian Survival Analysis of Transition to Schizophrenia and Remission in Individuals at Clinical High-Risk for Psychosis (CHR-P) (조현병 베이지안 생존분석)/🟩2.Methods/🟨1.데이터 설명/🟪1.생존 데이터 처리/🟧데이터1/attachments/survival.csv'
df <- read.csv(path_data, stringsAsFactors = FALSE)

# =================== time / event_type 생성 ===================
df2 <- df %>%
  mutate(
    time = pmin(dur_transition, dur_remission),
    event_type = as.integer(dplyr::case_when(
      transition == 1 ~ 1L,   # transition
      remission  == 1 ~ 2L,   # remission
      TRUE            ~ 0L    # censored
    ))
  ) %>%
  arrange(time)

# =================== (선택) 데이터 수정/보정 예시 ===================
df2 <- df2 |>
  mutate(
    dur_remission = ifelse(id == 20092201, 730, dur_remission)  # 필요 시만
  )

# =================== factor 라벨 ===================
df2$event_status <- factor(df2$event_type, levels = c(0,1,2),
                           labels = c("censored","transition","remission"))
df2$sex_factor   <- factor(df2$sex, levels = c(0,1), labels = c("Female","Male"))

# =================== 단위 명시 리네임 ===================
df3 <- df2 %>%
  rename(
    time_day           = time,
    dur_transition_day = dur_transition,
    dur_remission_day  = dur_remission
  )

# =================== 관점별 경과일 & 도달 나이(정수, 365 기준) ===================
df4 <- df3 %>%
  mutate(
    t_end_day       = time_day,
    event_tr_def    = as.integer(event_type == 1L),
    event_re_def    = as.integer(event_type == 2L),
    event_any_def   = as.integer(event_type %in% c(1L,2L)),
    
    tr_elapsed_day  = ifelse(event_tr_def == 1L, dur_transition_day, t_end_day),
    re_elapsed_day  = ifelse(event_re_def == 1L, dur_remission_day,  t_end_day),
    any_elapsed_day = t_end_day,
    
    tr_age_add_years  = trunc(tr_elapsed_day  / 365),
    re_age_add_years  = trunc(re_elapsed_day  / 365),
    any_age_add_years = trunc(any_elapsed_day / 365),
    
    tr_age_attained_int  = age + tr_age_add_years,   # transition 기준 (데이터1 고유)
    re_age_attained_int  = age + re_age_add_years,   # remission  기준
    any_age_attained_int = age + any_age_add_years   # 첫 사건/검열 기준
  )

# =================== ▶ 데이터셋2의 `out` 스키마에 맞춘 정렬/추가 ===================
# 데이터1에는 name/hid/birth/enter/date가 없음 → NA로 채움
# day/year는 데이터2 정의에 맞추어 생성: day = time_day, year = day/365
# status/event_status는 "transition 발생 여부" 기준(데이터2와 동일 의미)
out1 <- df4 %>%
  mutate(
    name  = NA_character_,
    hid   = NA_integer_,
    birth = as.Date(NA), enter = as.Date(NA), date = as.Date(NA),
    
    day   = as.integer(time_day),
    year  = day / 365,
    
    # 데이터2의 검증 컬럼 형식 맞춤
    day_calc          = NA_integer_,
    day_mismatch      = NA,                   # 계산 불가 → NA
    year_calc_365     = day / 365,
    year_calc_36525   = day / 365.25,
    year_mismatch_365 = FALSE,               # 정의상 year == year_calc_365 이므로 FALSE로 둠
    
    status       = as.integer(event_type == 1L),                       # 1=transition, 0=censored
    event_status = factor(ifelse(status == 1L, "transition", "censored"),
                          levels = c("censored","transition")),
    
    # 생일 기준 나이(정확)는 계산 불가 → NA
    age_entry_cont_36525    = NA_real_,
    age_attained_cont_36525 = NA_real_,
    
    # 정수 나이: entry=제공된 정수 age, attained=첫 사건/검열 기준(any)
    age_entry_int    = as.integer(age),
    age_attained_int = as.integer(any_age_attained_int),
    
    # 365일 기준 비교용(데이터2와 동일 정의)
    age_add_years_int_365           = trunc(day / 365),
    age_attained_from_entry_int_365 = age_entry_int + age_add_years_int_365,
    attained_int_diff  = age_attained_int - age_attained_from_entry_int_365,
    attained_diff_flag = attained_int_diff != 0,
    
    # 데이터2와 동일 네이밍(transition 기준)
    # (데이터2에서는 tr_age_attained_int == age_attained_int 이지만,
    #  데이터1은 remission 선행 가능 → 다를 수 있음. 그대로 유지)
    t_end_day        = t_end_day,
    event_tr_def     = event_tr_def,
    tr_elapsed_day   = tr_elapsed_day,
    tr_age_add_years = tr_age_add_years,
    tr_age_attained_int = tr_age_attained_int,
    
    # exact 칼럼(placeholder): 현재는 정수 기반 값을 '임시 연속형'으로 복사
    dataset_source = "dataset1",
    age_entry_cont_exact           = as.numeric(age_entry_int),
    age_entry_int_exact            = as.integer(age_entry_int),
    tr_age_attained_cont_exact     = as.numeric(tr_age_attained_int),
    tr_age_attained_int_exact      = as.integer(tr_age_attained_int),
    any_age_attained_cont_exact    = as.numeric(any_age_attained_int),
    any_age_attained_int_exact     = as.integer(any_age_attained_int),
    re_age_attained_cont_exact     = as.numeric(re_age_attained_int),
    re_age_attained_int_exact      = as.integer(re_age_attained_int),
    age_exact_flag_cont = FALSE,
    age_exact_flag_int  = FALSE,
    age_exact_method    = "placeholder: baseline_integer + trunc(elapsed/365)"
  ) %>%
  # 데이터셋2의 out 컬럼 순서와 동일하게 정렬
  select(
    id, name, hid,
    birth, enter, date,
    day, year,
    day_calc, day_mismatch,
    year_calc_365, year_calc_36525, year_mismatch_365,
    status, event_status,
    sex, sex_factor,
    age_entry_cont_36525, age_attained_cont_36525,
    age_entry_int,        age_attained_int,
    age_add_years_int_365, age_attained_from_entry_int_365,
    attained_int_diff, attained_diff_flag,
    t_end_day, event_tr_def, tr_elapsed_day, tr_age_add_years, tr_age_attained_int,
    dataset_source,
    age_entry_cont_exact, age_entry_int_exact,
    tr_age_attained_cont_exact, tr_age_attained_int_exact,
    any_age_attained_cont_exact, any_age_attained_int_exact,
    re_age_attained_cont_exact,  re_age_attained_int_exact,
    age_exact_flag_cont, age_exact_flag_int, age_exact_method
  )




# =================== 저장 ===================
save_path1 <- '/Volumes/ObsidianVault/🟩Obsidian/☔️Papers_Writing(논문 쓰기)/✅Currently working/⬛A Bayesian Survival Analysis of Transition to Schizophrenia and Remission in Individuals at Clinical High-Risk for Psychosis (CHR-P) (조현병 베이지안 생존분석)/🟩2.Methods/🟨1.데이터 설명/🟪1.생존 데이터 처리/🟧데이터1/attachments'
write.csv(out1, file.path(save_path1, "dataset1_processed.csv"), row.names = FALSE)






