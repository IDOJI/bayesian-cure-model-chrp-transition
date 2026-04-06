# *🗂️ Data3 (Merged Data1+Data2) Data Dictionary 📘🧩*

> ***Target file (assumed)**: `MERGED_dataset3_pnu_snu.csv`*  
> ***Composition**: preprocessed Data1 (PNU) + preprocessed Data2 (SNU), merged using `bind_rows()`*  
> ***Unit of observation**: 1 subject = 1 row*  
> ***Time definition (Day 0)**: `date_entry` is defined as day 0, and `days_followup` is calculated as `date_followup - date_entry` (elapsed days: 0/1/2/...) ✅⏱️*

*---*

# *🟧 Core Analysis Schema*

## *🟩 `id` (character)*

- ***Meaning**: Subject ID*
    
- ***Source***
    
    - *PNU: raw `id`*
        
    - *SNU: raw `id`, converted to character*
        
- ***Note**: After merging, **`id` alone may not be a unique key** → using `site + id` as the key is recommended ⚠️*
    

*---*

## *🟩 `site` (character)*

- ***Meaning**: Data source label (center/cohort)*
    
- ***Example values**: `"PNU"`, `"SNU"` (or `"CHR-P"`, etc.)*
    
- ***How created***
    
    - *PNU: constant value `"PNU"`*
        
    - *SNU: `site_label` assigned during preprocessing*
        

*---*

# *🟧 Sex Variables*

## *🟩 `sex_num` (numeric; 0/1)*

- ***Meaning**: Numeric coding for sex*
    
- ***Coding***
    
    - *`0` = Female*
        
    - *`1` = Male*
        
- ***Note**: Missing or other codes are stored as `NA`*
    

## *🟩 `sex_fact` (factor; levels = Female, Male)*

- ***Meaning**: Labeled factor for sex*
    
- ***Levels (fixed)**: `c("Female","Male")`*
    
- ***Created by**: mapping `sex_num` to labels, then converting to factor*
    

*---*

# *🟧 Date Variables*

## *🟩 `date_birth` (Date)*

- ***Meaning**: Date of birth*
    
- ***Source***
    
    - *PNU: raw `birth`*
        
    - *SNU: raw `birth`*
        
- ***Transformation**: parsed with `ymd()` and converted to `Date`*
    

## *🟩 `date_entry` (Date)*

- ***Meaning**: Observation start date (entry/baseline)*
    
- ***Source***
    
    - *PNU: raw `mri` → `date_entry`*
        
    - *SNU: raw `enter` → `date_entry`*
        
- ***Transformation**: parsed with `ymd()` and converted to `Date`*
    

*---*

# *🟧 Age Variables*

## *🟩 `age_int` (numeric/integer; years)*

- ***Meaning**: Integer age at entry*
    
- ***Calculation**:*  
    *`year(date_entry) - year(date_birth) - I(MMDD(entry) < MMDD(birth))`*
    

## *🟩 `age_exact_entry` (numeric; years)*

- ***Meaning**: Exact age at entry (including fractional year)*
    
- ***Calculation**: computed as the day fraction occupied by the entry date within the interval from one birthday to the next*
    

## *🟩 `age_exact_followup` (numeric; years)*

- ***Meaning**: Exact age at the end of follow-up*
    
- ***Derived from***
    
    - *`date_followup = date_entry + days(days_followup)`*
        
- ***Calculation**: computed at `date_followup` using the same method as `age_exact_entry`*
    

*---*

# *🟧 Survival Analysis Variables*

## *🟩 `days_followup` (numeric; days)*

- ***Meaning**: Follow-up duration in days — time-to-event or time-to-censoring*
    
- ***Day 0 definition***
    
    - *the same day as `date_entry` = 0*
        
    - *the day after entry = 1*
        
    - *that is, `date_followup - date_entry`*
        
- ***Source / creation***
    
    - ***PNU**:*
        
        - *if transition: `dur_transition`*
            
        - *if remission (`status_num == 2`, analyzed as right censoring): `dur_remission`*
            
        - *if censored (both = 0): use `dur_transition` (preprocessing QC enforces `dur_transition == dur_remission`)*
            
    - ***SNU**:*
        
        - *`days_followup = as.numeric(date_end - date_entry)` (QC confirmed consistency with raw `day`)*
            
- ***Recommended QC rules***
    
    - *`days_followup >= 0`*
        
    - *`NA` not allowed (must be blocked before analysis)*
        

*---*

## *🟩 `status_num` (integer; status code)*

- ***Meaning**: Numeric coding for raw event status*
    
- ***Coding***
    
    - *`0` = right censoring*
        
    - *`1` = transition*
        
    - *`2` = remission*

- ***Analysis convention**: In the integrated/main transition-only analysis, `status_num == 2` is treated as right censoring*
        
- ***Site-specific difference***
    
    - *PNU: raw `0/1/2` can all appear; in analysis, `2` is collapsed into right censoring*
        
    - *SNU: only `0/1` exist (transition vs censoring); `2` is structurally absent*
        

*---*

## *🟩 `status` (factor; raw levels = right_censoring, remission, transition)*

- ***Meaning**: Labeled factor for raw event status*
    
- ***Labels***
    
    - *`0` → `"right_censoring"`*
        
    - *`1` → `"transition"`*
        
    - *`2` → `"remission"`*
        
- ***Levels (fixed order)**: `c("right_censoring","remission","transition")`*
    
- ***Note**: The `"remission"` label is retained for raw-data traceability, but in the integrated/main transition-only analysis it is treated the same as `"right_censoring"`*

- ***Note**: In SNU, the `"remission"` level exists but has no observed cases (empty level) — it is still included for merge stability ✅*
    

*---*

# *🟧 Derived / Reference Variable (if needed)*

## *🟨 `date_followup` (derived; Date)*

- ***Definition**: `date_entry + days(days_followup)`*
    
- ***Not stored as a column** (calculated only when needed)*
    

*---*

# *🟧 Interpretation Notes for Analysis Design ⚠️*

## *1) Raw remission code exists only in PNU*

- *PNU raw data can contain `status_num == 2` / `"remission"`, but the integrated/main analysis treats this as right censoring rather than as a separate event*
    
- *SNU has only `0/1`, so the integrated dataset is aligned as a transition-only analysis with censoring*
    

*➡️ In integrated modeling for this project, the analysis is standardized to a **transition-only analysis**, and remission-coded observations are handled as right censoring.*
    

## *2) `site + id` is recommended as the identification key*

- *The same `id` may appear in different sites*
    

*---*

# *✅ Final One-Line Definition 🧠*

*This merged dataset is a **standardized schema for center-integrated survival analysis**, using*  
***`days_followup` (elapsed days since Day 0) + `status/status_num` (with `1` as transition and `0/2` treated as censoring in the main analysis)** as the core survival variables,*  
*with `age_*` variables harmonized as covariates derived from `date_birth` and `date_entry`. ✅📌*
