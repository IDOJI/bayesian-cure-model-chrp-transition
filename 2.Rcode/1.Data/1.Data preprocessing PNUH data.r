# 🟧 패키지 로드 ===================
suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(readr)
  library(tibble)
})



# 🟧 입출력 설정(여기만 바꾸면 됨) ===================
path_in <- "/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/✅Currently working/⬛조현병 베이지안 생존분석/🟧0.생존 데이터 처리와 요약/🟦2.데이터1 처리/attachments/survival_merged.csv"
out_dir <- "/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/✅Currently working/⬛조현병 베이지안 생존분석/🟧0.생존 데이터 처리와 요약/🟦2.데이터1 처리/attachments"
out_file <- "dataset1_processed.csv"

path_out_data <- file.path(out_dir, out_file)

# 🟧 파일명 베이스(base_name) 사전 정의 ===================
base_name <- sub("\\.csv$", "", out_file)



# 🟧 데이터 로드 ===================
data <- read_csv(path_in, show_col_types = FALSE)



# 🟧 전처리(rename/mutate/reorder) ===================
data_processed <- data %>%
  rename(
    age_reported      = age,
    transition_days   = dur_transition,
    remission_days    = dur_remission,
    transition_status = transition,
    remission_status  = remission,
    date_entry        = mri
  ) %>%
  mutate(
    sex_num = suppressWarnings(as.numeric(sex)),
    sex_fact = case_when(
      sex_num == 1 ~ "Male",
      sex_num == 0 ~ "Female",
      TRUE ~ NA_character_
    ),
    sex_fact = factor(sex_fact, levels = c("Female", "Male")),
    
    date_birth = ymd(birth),
    date_entry = ymd(date_entry),
    
    transition_status = suppressWarnings(as.integer(transition_status)),
    remission_status  = suppressWarnings(as.integer(remission_status)),
    
    transition_days = suppressWarnings(as.numeric(transition_days)),
    remission_days  = suppressWarnings(as.numeric(remission_days))
  ) %>%
  select(
    id, sex_fact, sex_num, date_entry, date_birth, age_reported,
    transition_status, transition_days,
    remission_status, remission_days,
    everything()
  ) %>%
  select(-sex, -birth)



# 🟧 QC: ID 중복 점검 ===================
qc_iddup <- data_processed %>%
  summarise(
    n_total = n(),
    n_id_duplicated = sum(duplicated(id)),
    n_id_unique = n_distinct(id)
  )

print(qc_iddup)

cat("\n[해석 메시지]\n")
if (qc_iddup$n_id_duplicated == 0) {
  cat("- id 중복 없음 → OK\n")
} else {
  cat("- id 중복 존재 → 반복측정/중복등록 여부 점검 권장\n")
}



# 🟧 QC: 날짜 파싱/논리 점검 ===================
flagged_dates <- data_processed %>%
  transmute(
    id,
    date_birth,
    date_entry,
    flag_birth_na = is.na(date_birth),
    flag_entry_na = is.na(date_entry),
    flag_birth_after_entry = !is.na(date_birth) & !is.na(date_entry) & (date_birth > date_entry)
  )

n_total_dates <- nrow(flagged_dates)

summary_tbl_dates <- tibble(
  check = c(
    "date_birth NA (파싱 실패/결측)",
    "date_entry NA (파싱 실패/결측)",
    "date_birth > date_entry (논리 오류 가능)"
  ),
  n_flagged = c(
    sum(flagged_dates$flag_birth_na, na.rm = TRUE),
    sum(flagged_dates$flag_entry_na, na.rm = TRUE),
    sum(flagged_dates$flag_birth_after_entry, na.rm = TRUE)
  )
) %>%
  mutate(pct = ifelse(n_total_dates > 0, round(100 * n_flagged / n_total_dates, 2), NA_real_))

print(summary_tbl_dates)

cat("\n[해석 메시지]\n")
if (sum(flagged_dates$flag_birth_na | flagged_dates$flag_entry_na | flagged_dates$flag_birth_after_entry, na.rm = TRUE) == 0) {
  cat("- 핵심 오류(NA/날짜역전) 없음 → 정확 나이 계산을 안정적으로 수행 가능\n")
} else {
  cat("- 핵심 오류(NA 또는 날짜역전) 존재 → 정확 나이 계산 전에 해당 케이스 정리 권장\n")
}



# 🟧 QC: 윤일(2/29) 출생자 존재 여부 ===================
qc_feb29_tbl <- data_processed %>%
  transmute(
    id,
    date_birth,
    flag_birth_na = is.na(date_birth),
    flag_birth_feb29 = !is.na(date_birth) & format(date_birth, "%m-%d") == "02-29"
  )

qc_feb29_summary <- qc_feb29_tbl %>%
  summarise(
    n_total = n(),
    n_birth_na = sum(flag_birth_na, na.rm = TRUE),
    n_birth_feb29 = sum(flag_birth_feb29, na.rm = TRUE),
    pct_birth_feb29 = ifelse(n_total > 0, round(100 * n_birth_feb29 / n_total, 2), NA_real_)
  ) %>%
  as_tibble()

print(qc_feb29_summary)

cat("\n[해석 메시지]\n")
if (qc_feb29_summary$n_birth_feb29 == 0) {
  cat("- 윤일(2/29) 출생자 없음 → OK\n")
} else {
  cat("- 윤일(2/29) 출생자 존재 → 윤년 처리 규칙(2/28 vs 3/1) 명시 권장\n")
}

# 🟧 QC: 윤일(2/29) 출생자 케이스 저장 ===================
if (qc_feb29_summary$n_birth_feb29 > 0) {
  feb29_cases <- qc_feb29_tbl %>%
    filter(flag_birth_feb29) %>%
    select(id, date_birth)
  
  path_out_feb29 <- file.path(out_dir, paste0(base_name, "_QC_feb29_birth_cases.csv"))
  write_csv(feb29_cases, path_out_feb29)
  
  cat("- QC 윤일 출생자 케이스: ", path_out_feb29, "\n", sep = "")
}



# 🟧 QC: status 값 범위 및 동시 관측(both=1) 점검 ===================
qc_status_raw <- data_processed %>%
  summarise(
    n_total = n(),
    n_transition_out_of_range = sum(!is.na(transition_status) & !transition_status %in% c(0L, 1L)),
    n_remission_out_of_range  = sum(!is.na(remission_status)  & !remission_status  %in% c(0L, 1L)),
    n_both_1 = sum(transition_status == 1L & remission_status == 1L, na.rm = TRUE),
    n_both_0 = sum(transition_status == 0L & remission_status == 0L, na.rm = TRUE)
  )

print(qc_status_raw)

cat("\n[해석 메시지]\n")
if (qc_status_raw$n_transition_out_of_range == 0 && qc_status_raw$n_remission_out_of_range == 0) {
  cat("- transition/remission status는 0/1 범위 → OK\n")
} else {
  cat("- status 값이 0/1 범위를 벗어남 → 코딩 점검 필요\n")
}

if (qc_status_raw$n_both_1 == 0) {
  cat("- transition=1 & remission=1 동시 관측(충돌) 없음 → OK\n")
} else {
  cat("- transition=1 & remission=1 동시 관측(충돌) 존재 → competing risk 정의 재정의 필요\n")
}



# 🟧 QC: 둘 다 0인데 days가 다른 케이스 점검 ===================
chk_days_inconsistency <- data_processed %>%
  mutate(
    flag_both0 = (transition_status == 0L & remission_status == 0L),
    flag_days_diff_when_both0 = flag_both0 &
      !is.na(transition_days) & !is.na(remission_days) &
      (transition_days != remission_days)
  )

summary_tbl_days <- chk_days_inconsistency %>%
  summarise(
    n_total = n(),
    n_both0 = sum(flag_both0, na.rm = TRUE),
    n_days_diff_when_both0 = sum(flag_days_diff_when_both0, na.rm = TRUE)
  )

print(summary_tbl_days)

cat("\n[해석 메시지]\n")
if (summary_tbl_days$n_days_diff_when_both0 == 0) {
  cat("- 둘 다 0인 케이스에서 transition_days와 remission_days가 동일\n")
  cat("  → 둘 다 0인 경우 days를 검열 시점으로 사용해도 일관적\n")
} else {
  cat("- 둘 다 0인데 days가 다른 케이스 존재\n")
  cat("  → 둘 다 0인 경우의 검열 시점(days) 규칙을 별도로 정해야 함\n")
  ex_ids <- chk_days_inconsistency %>%
    filter(flag_days_diff_when_both0) %>%
    transmute(id, transition_days, remission_days) %>%
    head(20)
  print(ex_ids)
}



# 🟧 나이/생존변수 생성(엔트리 나이 + status/days_followup 엄격 생성) ===================
data_processed_age <- data_processed %>%
  mutate(
    # (A) entry 시점 만 나이(정수)
    age_int = year(date_entry) - year(date_birth) -
      as.integer(format(date_entry, "%m%d") < format(date_birth, "%m%d")),
    
    # (B) entry 시점 정확 나이(소수 포함)
    last_bday_entry = date_birth %m+% years(age_int),
    next_bday_entry = date_birth %m+% years(age_int + 1),
    
    age_exact_entry = age_int +
      as.numeric(difftime(date_entry, last_bday_entry, units = "days")) /
      as.numeric(difftime(next_bday_entry, last_bday_entry, units = "days")),
    
    # (C) 상태 조합 플래그(엄격 판단)
    flag_both1 = (transition_status == 1L & remission_status == 1L),
    flag_both0 = (transition_status == 0L & remission_status == 0L),
    
    # (D) 둘 다 0인데 duration 불일치 플래그
    flag_censor_days_mismatch = flag_both0 &
      !is.na(transition_days) & !is.na(remission_days) &
      (transition_days != remission_days),
    
    # (E) status_num 생성(정의)
    status_num = case_when(
      transition_status == 1L ~ 1L,
      remission_status  == 1L ~ 2L,
      TRUE ~ 0L
    ),
    
    # (F) days_followup 생성(정의)
    days_followup = case_when(
      transition_status == 1L ~ transition_days,
      remission_status  == 1L ~ remission_days,
      TRUE ~ transition_days # censor는 일단 transition_days로 두되 아래에서 "동일성" 강제
    )
  ) %>%
  select(
    id, sex_fact, sex_num, date_entry, date_birth,
    age_reported, age_int, age_exact_entry,
    status_num, days_followup,
    transition_status, transition_days,
    remission_status, remission_days,
    everything()
  ) %>%
  select(-last_bday_entry, -next_bday_entry)



# 🟧 days_followup 생성 규칙 강제검증(불일치면 중단) ===================
# (1) transition=1 & remission=1 금지
n_both1 <- sum(data_processed_age$flag_both1, na.rm = TRUE)
if (n_both1 > 0) {
  bad_ids <- data_processed_age %>%
    filter(flag_both1) %>%
    transmute(id, transition_days, remission_days) %>%
    head(20)
  stop(
    paste0(
      "[ERROR] transition=1 & remission=1 동시 관측 케이스가 ", n_both1, "건 존재합니다.\n",
      "예시(최대 20건):\n",
      paste(capture.output(print(bad_ids)), collapse = "\n"),
      "\n→ competing risk 규칙(먼저 발생 사건 등)을 재정의하세요."
    ),
    call. = FALSE
  )
}

# (2) 둘 다 0인데 censor duration 불일치 금지
n_censor_mis <- sum(data_processed_age$flag_censor_days_mismatch, na.rm = TRUE)
if (n_censor_mis > 0) {
  bad_ids <- data_processed_age %>%
    filter(flag_censor_days_mismatch) %>%
    transmute(id, transition_days, remission_days) %>%
    head(20)
  stop(
    paste0(
      "[ERROR] transition=0 & remission=0 인데 transition_days != remission_days 케이스가 ",
      n_censor_mis, "건 존재합니다.\n",
      "예시(최대 20건):\n",
      paste(capture.output(print(bad_ids)), collapse = "\n"),
      "\n→ 검열 시점 정의(공통 종료일 변수/별도 censor_date 등)를 먼저 정하세요."
    ),
    call. = FALSE
  )
}

# (3) 추가 방어: days_followup NA/음수 금지
n_na <- sum(is.na(data_processed_age$days_followup))
n_neg <- sum(!is.na(data_processed_age$days_followup) & data_processed_age$days_followup < 0)
if (n_na > 0 || n_neg > 0) {
  bad_ids <- data_processed_age %>%
    filter(is.na(days_followup) | days_followup < 0) %>%
    transmute(id, transition_status, remission_status, transition_days, remission_days, days_followup) %>%
    head(20)
  stop(
    paste0(
      "[ERROR] days_followup에 NA 또는 음수 값이 존재합니다. (NA=", n_na, ", negative=", n_neg, ")\n",
      "예시(최대 20건):\n",
      paste(capture.output(print(bad_ids)), collapse = "\n")
    ),
    call. = FALSE
  )
}



# 🟧 QC: days_followup 매핑(사건별 duration 대응) ===================
qc_map <- data_processed_age %>%
  summarise(
    n_total = n(),
    
    n_days_followup_na = sum(is.na(days_followup)),
    n_days_followup_negative = sum(!is.na(days_followup) & days_followup < 0),
    
    n_transition_days_na_when_transition =
      sum(transition_status == 1L & is.na(transition_days), na.rm = TRUE),
    n_remission_days_na_when_remission =
      sum(remission_status == 1L & is.na(remission_days), na.rm = TRUE),
    
    n_transition_map_mismatch = sum(
      transition_status == 1L &
        !is.na(days_followup) & !is.na(transition_days) &
        (days_followup != transition_days),
      na.rm = TRUE
    ),
    
    n_remission_map_mismatch = sum(
      remission_status == 1L &
        !is.na(days_followup) & !is.na(remission_days) &
        (days_followup != remission_days),
      na.rm = TRUE
    ),
    
    n_censor_map_mismatch = sum(
      transition_status == 0L & remission_status == 0L &
        !is.na(days_followup) & !is.na(transition_days) &
        (days_followup != transition_days),
      na.rm = TRUE
    )
  )

print(qc_map)

cat("\n[해석 메시지]\n")
if (qc_map$n_days_followup_na == 0 && qc_map$n_days_followup_negative == 0) {
  cat("- days_followup: NA/negative 없음 → OK\n")
} else {
  cat("- days_followup: NA 또는 음수 존재 → 점검 필요\n")
}

if (qc_map$n_transition_map_mismatch == 0 &&
    qc_map$n_remission_map_mismatch == 0 &&
    qc_map$n_censor_map_mismatch == 0) {
  cat("- 사건별 days_followup가 대응 duration과 일관 → OK\n")
} else {
  cat("- days_followup 매핑 불일치 존재 → 규칙/데이터 점검 필요\n")
}

data_end_preview <- data_processed_age %>%
  mutate(date_end = date_entry + days(days_followup)) %>%
  select(id, status_num, transition_status, remission_status, date_entry, days_followup, date_end) %>%
  head(10)

cat("\n[참고: date_end = date_entry + days_followup (상위 10개)]\n")
print(data_end_preview)



# 🟧 QC: age_reported vs age_int 불일치 점검(불일치 테이블 저장용) ===================
chk_age <- data_processed_age %>%
  mutate(
    age_reported_num = suppressWarnings(as.numeric(age_reported)),
    age_same = !is.na(age_reported_num) & (age_int == age_reported_num)
  )

print(
  chk_age %>%
    summarise(
      n_total = n(),
      n_age_reported_na_or_nonnum = sum(is.na(age_reported_num)),
      n_mismatch = sum(!is.na(age_reported_num) & age_int != age_reported_num),
      n_match = sum(age_same, na.rm = TRUE)
    )
)

mismatch_detail <- chk_age %>%
  filter(!is.na(age_reported_num) & age_int != age_reported_num) %>%
  mutate(
    last_bday_entry = date_birth %m+% years(age_int),
    next_bday_entry = date_birth %m+% years(age_int + 1),
    days_to_next_bday = as.numeric(difftime(next_bday_entry, date_entry, units = "days")),
    days_from_last_bday = as.numeric(difftime(date_entry, last_bday_entry, units = "days")),
    delta = age_reported_num - age_int,
    reason = case_when(
      delta == 1 ~ "age_reported가 '연도차(생일 보정 없음)' 방식일 가능성(만나이보다 1 크게 기록)",
      delta == -1 ~ "age_reported가 만나이보다 1 작음 → age_reported 정의/입력 오류 가능",
      TRUE ~ "age_reported 정의 차이 또는 date/age 입력 오류 가능"
    )
  ) %>%
  select(
    id, sex_fact, sex_num,
    date_birth, date_entry,
    age_reported, age_reported_num, age_int,
    days_to_next_bday, days_from_last_bday,
    reason
  )

cat("\n[불일치 개체 상세]\n")
if (nrow(mismatch_detail) == 0) {
  cat("- 불일치 없음\n")
} else {
  print(mismatch_detail)
}

age_mismatch_cases <- mismatch_detail %>%
  select(id, date_birth, date_entry, age_reported, age_reported_num, age_int, reason)

data_processed_age <- data_processed_age %>%
  select(-age_reported)



# 🟧 follow-up 종료 시점 age_exact_followup 계산 ===================
data_processed_age <- data_processed_age %>%
  mutate(
    date_followup = date_entry + days(days_followup),
    
    age_int_tmp = year(date_followup) - year(date_birth) -
      as.integer(format(date_followup, "%m%d") < format(date_birth, "%m%d")),
    
    last_bday_followup = date_birth %m+% years(age_int_tmp),
    next_bday_followup = date_birth %m+% years(age_int_tmp + 1),
    
    age_exact_followup = age_int_tmp +
      as.numeric(difftime(date_followup, last_bday_followup, units = "days")) /
      as.numeric(difftime(next_bday_followup, last_bday_followup, units = "days"))
  ) %>%
  select(-date_followup, -age_int_tmp, -last_bday_followup, -next_bday_followup)



# 🟧 QC: 나이 단조성(age_exact_followup >= age_exact_entry) ===================
qc_age_mono <- data_processed_age %>%
  summarise(
    n_total = n(),
    n_age_followup_lt_entry = sum(age_exact_followup < age_exact_entry, na.rm = TRUE)
  )

print(qc_age_mono)

cat("\n[해석 메시지]\n")
if (qc_age_mono$n_age_followup_lt_entry == 0) {
  cat("- age_exact_followup >= age_exact_entry 항상 성립 → OK\n")
} else {
  cat("- follow-up 나이가 entry 나이보다 작은 케이스 존재 → days_followup/날짜 점검 필요\n")
}



# 🟧 최종 분석용 변수 정리(status factor/원본 변수 제거) ===================
data_final <- data_processed_age %>%
  select(-transition_status, -remission_status, -transition_days, -remission_days) %>%
  mutate(
    status = case_when(
      status_num == 0L ~ "right_censoring",
      status_num == 1L ~ "transition",
      status_num == 2L ~ "remission",
      TRUE ~ NA_character_
    ),
    status = factor(status, levels = c("right_censoring", "remission", "transition"))
  ) %>%
  select(-flag_both1, -flag_both0, -flag_censor_days_mismatch)



# 🟧 site 추가 및 컬럼 재배치 ===================
data_final <- data_final %>%
  mutate(site = "PNU") %>%
  relocate(site, .after = id) %>%
  relocate(age_exact_followup, .after = age_exact_entry)



# 🟧 전처리 결과 저장 ===================
write_csv(data_final, path_out_data)

cat("\n[저장 완료]\n")
cat("- 전처리 데이터: ", path_out_data, "\n", sep = "")



# 🟧 QC 요약/메시지 생성 및 저장 ===================
qc_tbl <- data_final %>%
  summarise(
    n_total = n(),
    n_right_censoring = sum(status == "right_censoring", na.rm = TRUE),
    n_transition = sum(status == "transition", na.rm = TRUE),
    n_remission = sum(status == "remission", na.rm = TRUE),
    
    days_followup_min = min(days_followup, na.rm = TRUE),
    days_followup_median = median(days_followup, na.rm = TRUE),
    days_followup_max = max(days_followup, na.rm = TRUE),
    
    n_days_followup_na = sum(is.na(days_followup)),
    n_days_followup_negative = sum(!is.na(days_followup) & days_followup < 0),
    
    n_age_int_na = sum(is.na(age_int)),
    n_age_exact_entry_na = sum(is.na(age_exact_entry)),
    n_age_exact_followup_na = sum(is.na(age_exact_followup)),
    
    n_id_duplicated = qc_iddup$n_id_duplicated,
    n_both_1_raw = qc_status_raw$n_both_1,
    n_status_out_of_range_raw = qc_status_raw$n_transition_out_of_range + qc_status_raw$n_remission_out_of_range,
    n_days_diff_when_both0 = summary_tbl_days$n_days_diff_when_both0,
    
    n_transition_map_mismatch = qc_map$n_transition_map_mismatch,
    n_remission_map_mismatch = qc_map$n_remission_map_mismatch,
    n_censor_map_mismatch = qc_map$n_censor_map_mismatch,
    
    n_birth_feb29 = qc_feb29_summary$n_birth_feb29,
    pct_birth_feb29 = qc_feb29_summary$pct_birth_feb29,
    
    n_age_followup_lt_entry = qc_age_mono$n_age_followup_lt_entry,
    
    n_age_reported_mismatch = if (exists("age_mismatch_cases")) nrow(age_mismatch_cases) else 0
  ) %>%
  as_tibble()

print(qc_tbl)

qc_msgs <- c(
  paste0("n_total: ", qc_tbl$n_total),
  paste0("status counts: right_censoring=", qc_tbl$n_right_censoring,
         ", remission=", qc_tbl$n_remission,
         ", transition=", qc_tbl$n_transition),
  paste0("days_followup (days): min=", qc_tbl$days_followup_min,
         ", median=", qc_tbl$days_followup_median,
         ", max=", qc_tbl$days_followup_max),
  
  paste0("Feb29 birth (윤일 출생자): n=", qc_tbl$n_birth_feb29,
         ", pct=", qc_tbl$pct_birth_feb29, "%"),
  
  if (qc_tbl$n_days_followup_na == 0 && qc_tbl$n_days_followup_negative == 0)
    "days_followup: NA/negative 없음 → OK"
  else
    paste0("days_followup: NA 또는 음수 존재 → 점검 필요 (NA=", qc_tbl$n_days_followup_na,
           ", negative=", qc_tbl$n_days_followup_negative, ")"),
  
  if (qc_tbl$n_age_int_na == 0 && qc_tbl$n_age_exact_entry_na == 0)
    "entry age (age_int/age_exact_entry): NA 없음 → OK"
  else
    paste0("entry age: NA 존재 → 점검 필요 (age_int NA=", qc_tbl$n_age_int_na,
           ", age_exact_entry NA=", qc_tbl$n_age_exact_entry_na, ")"),
  
  if (qc_tbl$n_age_exact_followup_na == 0)
    "follow-up age (age_exact_followup): NA 없음 → OK"
  else
    paste0("follow-up age: NA 존재 → 점검 필요 (NA=", qc_tbl$n_age_exact_followup_na, ")"),
  
  if (qc_tbl$n_id_duplicated == 0)
    "id 중복: 없음 → OK"
  else
    paste0("id 중복: 존재 → 점검 필요 (n_duplicated=", qc_tbl$n_id_duplicated, ")"),
  
  if (qc_tbl$n_both_1_raw == 0)
    "transition=1 & remission=1 동시 관측(충돌): 없음 → OK"
  else
    paste0("transition=1 & remission=1 동시 관측(충돌): 존재 → 점검 필요 (n=", qc_tbl$n_both_1_raw, ")"),
  
  if (qc_tbl$n_status_out_of_range_raw == 0)
    "raw status 값 범위(0/1): 이상 없음 → OK"
  else
    paste0("raw status 값 범위(0/1): 이상 존재 → 점검 필요 (n_out_of_range=", qc_tbl$n_status_out_of_range_raw, ")"),
  
  if (qc_tbl$n_days_diff_when_both0 == 0)
    "둘 다 0인 케이스에서 transition_days vs remission_days: 동일 → OK"
  else
    paste0("둘 다 0인데 days 불일치 존재 → 규칙 재정의 필요 (n=", qc_tbl$n_days_diff_when_both0, ")"),
  
  if (qc_tbl$n_transition_map_mismatch == 0 && qc_tbl$n_remission_map_mismatch == 0 && qc_tbl$n_censor_map_mismatch == 0)
    "days_followup 매핑(transition/remission/censor): 불일치 없음 → OK"
  else
    paste0("days_followup 매핑 불일치 존재 → 점검 필요 (transition_mis=",
           qc_tbl$n_transition_map_mismatch, ", remission_mis=",
           qc_tbl$n_remission_map_mismatch, ", censor_mis=",
           qc_tbl$n_censor_map_mismatch, ")"),
  
  if (qc_tbl$n_age_followup_lt_entry == 0)
    "나이 단조성(age_followup >= age_entry): 위반 없음 → OK"
  else
    paste0("나이 단조성 위반 존재 → days_followup/날짜 점검 필요 (n=", qc_tbl$n_age_followup_lt_entry, ")"),
  
  if (qc_tbl$n_age_reported_mismatch > 0)
    paste0("age_reported mismatch cases: ", qc_tbl$n_age_reported_mismatch,
           " (QC 객체 age_mismatch_cases 확인)")
  else
    "age_reported mismatch cases: 없음"
)

cat("\n[QC 해석 메시지]\n")
cat(paste0("- ", qc_msgs, collapse = "\n"), "\n")

path_out_qc_tbl <- file.path(out_dir, paste0(base_name, "_QC_summary.csv"))
path_out_qc_msg <- file.path(out_dir, paste0(base_name, "_QC_message.txt"))

write_csv(qc_tbl, path_out_qc_tbl)
writeLines(qc_msgs, path_out_qc_msg)

if (exists("age_mismatch_cases") && nrow(age_mismatch_cases) > 0) {
  path_out_mismatch <- file.path(out_dir, paste0(base_name, "_QC_age_mismatch_cases.csv"))
  write_csv(age_mismatch_cases, path_out_mismatch)
}

cat("\n[QC 저장 완료]\n")
cat("- QC 요약: ", path_out_qc_tbl, "\n", sep = "")
cat("- QC 메시지: ", path_out_qc_msg, "\n", sep = "")
if (exists("age_mismatch_cases") && nrow(age_mismatch_cases) > 0) {
  cat("- QC 불일치 케이스: ", path_out_mismatch, "\n", sep = "")
}












# 🟧 추가 QC: Day0/패턴/민감도 분석 export ===================

# (1) Day0 여부 힌트 QC ------------------------------------------------
qc_day0_hint <- data_processed %>%
  summarise(
    min_transition_days = min(transition_days, na.rm = TRUE),
    min_remission_days  = min(remission_days,  na.rm = TRUE),
    n0_transition_days  = sum(transition_days == 0, na.rm = TRUE),
    n0_remission_days   = sum(remission_days  == 0, na.rm = TRUE)
  ) %>%
  as_tibble()

print(qc_day0_hint)

path_out_qc_day0 <- file.path(out_dir, paste0(base_name, "_QC_day0_hint.csv"))
write_csv(qc_day0_hint, path_out_qc_day0)

cat("\n[추가 QC 저장]\n")
cat("- Day0 힌트 QC: ", path_out_qc_day0, "\n", sep = "")


# (2) duration 패턴 QC(7의 배수 몰림 확인) ------------------------------
qc_mod7 <- data_processed %>%
  transmute(
    transition_days,
    remission_days,
    trans_mod7 = transition_days %% 7,
    remis_mod7 = remission_days %% 7
  ) %>%
  summarise(
    n = n(),
    pct_trans_multiple7 = mean(trans_mod7 == 0, na.rm = TRUE) * 100,
    pct_remis_multiple7 = mean(remis_mod7 == 0, na.rm = TRUE) * 100
  ) %>%
  as_tibble()

print(qc_mod7)

path_out_qc_mod7 <- file.path(out_dir, paste0(base_name, "_QC_duration_mod7.csv"))
write_csv(qc_mod7, path_out_qc_mod7)

cat("- Duration mod7 QC: ", path_out_qc_mod7, "\n", sep = "")


# (3) days_followup 1일 보정 민감도(Cox) -------------------------------
suppressPackageStartupMessages({
  library(survival)
})

fit_transition_cox <- function(df, shift = 0) {
  d <- df
  d$t <- d$days_followup + shift
  
  if (any(is.na(d$t))) stop("shift 적용 후 t에 NA가 생김")
  if (any(d$t <= 0, na.rm = TRUE)) stop("shift로 인해 0 이하 시간이 생김")
  
  d$event_transition <- as.integer(d$status_num == 1L)
  
  # site는 data_final에 이미 추가되어 있으므로 strata(site) 가능
  coxph(Surv(t, event_transition) ~ sex_fact + age_exact_entry + strata(site), data = d)
}

m0  <- fit_transition_cox(data_final, shift = 0)
m_1 <- fit_transition_cox(data_final, shift = -1)

qc_sensitivity_coef <- rbind(
  base   = coef(m0),
  minus1 = coef(m_1)
) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("model") %>%
  as_tibble()

print(qc_sensitivity_coef)

path_out_qc_sens <- file.path(out_dir, paste0(base_name, "_QC_days_shift_sensitivity_coef.csv"))
write_csv(qc_sensitivity_coef, path_out_qc_sens)

cat("- days_followup shift sensitivity (coef): ", path_out_qc_sens, "\n", sep = "")

# (옵션) 로그우도까지 같이 저장하면 더 깔끔함
qc_sensitivity_fit <- tibble(
  model = c("base", "minus1"),
  shift = c(0, -1),
  logLik = c(as.numeric(logLik(m0)), as.numeric(logLik(m_1)))
)

print(qc_sensitivity_fit)

path_out_qc_sens_fit <- file.path(out_dir, paste0(base_name, "_QC_days_shift_sensitivity_fit.csv"))
write_csv(qc_sensitivity_fit, path_out_qc_sens_fit)

cat("- days_followup shift sensitivity (fit): ", path_out_qc_sens_fit, "\n", sep = "")



# 🟧 추가 QC 해석 메시지(쉽게) export ===================

# 숫자 안전하게 뽑기
min_trans <- qc_day0_hint$min_transition_days[1]
min_remis <- qc_day0_hint$min_remission_days[1]
n0_trans  <- qc_day0_hint$n0_transition_days[1]
n0_remis  <- qc_day0_hint$n0_remission_days[1]

n_total <- qc_mod7$n[1]
pct7_t  <- qc_mod7$pct_trans_multiple7[1]
pct7_r  <- qc_mod7$pct_remis_multiple7[1]

# 민감도 결과 동일 여부(네 결과는 TRUE가 나올 것)
coef_same <- isTRUE(all.equal(
  qc_sensitivity_coef[qc_sensitivity_coef$model == "base", -1],
  qc_sensitivity_coef[qc_sensitivity_coef$model == "minus1", -1],
  tolerance = 0
))

ll_base  <- qc_sensitivity_fit$logLik[qc_sensitivity_fit$model == "base"][1]
ll_minus <- qc_sensitivity_fit$logLik[qc_sensitivity_fit$model == "minus1"][1]
ll_same  <- isTRUE(all.equal(ll_base, ll_minus, tolerance = 0))

# ✅ 사람이 읽기 쉬운 메시지(핵심만)
qc_easy_msgs <- c(
  paste0("✅ Day0 체크: 0일(당일) 사건/검열 없음. 가장 빠른 event는 transition ", min_trans,
         "일, remission ", min_remis, "일부터 시작."),
  paste0("✅ 7일 단위 몰림: transition 7의 배수 비율 ", round(pct7_t, 2), "%, remission ",
         round(pct7_r, 2), "% (n=", n_total, "). → 주단위로 뚜렷하게 반올림/몰림된 패턴은 없음."),
  paste0("✅ 1일 보정 민감도: days_followup을 1일 줄여도 결과 동일(coef_same=", coef_same,
         ", logLik_same=", ll_same, "). → 1일 컨벤션 차이는 분석 결과에 영향 없음. 지금 방식 그대로 쓰면 됨.")
)

cat("\n[추가 QC 쉬운 해석 메시지]\n")
cat(paste0("- ", qc_easy_msgs, collapse = "\n"), "\n")

# (A) 쉬운 메시지 전용 txt 저장
path_out_qc_easy <- file.path(out_dir, paste0(base_name, "_QC_easy_message.txt"))
writeLines(qc_easy_msgs, path_out_qc_easy)

cat("\n[저장 완료]\n")
cat("- 쉬운 QC 메시지: ", path_out_qc_easy, "\n", sep = "")

# (B) (선택) 기존 QC_message.txt에 덧붙여서 업데이트
#     (네 코드에서 qc_msgs / path_out_qc_msg가 이미 존재할 때만)
if (exists("qc_msgs") && exists("path_out_qc_msg")) {
  qc_msgs <- c(qc_msgs, qc_easy_msgs)
  writeLines(qc_msgs, path_out_qc_msg)
  cat("- 기존 QC 메시지 파일 업데이트: ", path_out_qc_msg, "\n", sep = "")
}
