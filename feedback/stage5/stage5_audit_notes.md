
# Stage 5 audit notes

## Overall judgment

Stage 5 is **mostly behaving as intended**. The core structural checks pass:

- 40 planned models were registered and all 40 converged.
- The subject × model × horizon prediction grid is complete: 78,600 expected rows and 78,600 observed rows.
- Risk and survival predictions stay within [0, 1], satisfy `risk + survival = 1`, and are horizon-monotone.
- Unsupported AUC/Brier/threshold outcome metrics are suppressed as designed.
- Horizon governance fields (`support_tier`, `horizon_evidence_class`, `claim_restriction_flag`, `risk_scale`) are appended in the long outputs.

## Findings that matter

### 1) One calibration regression failure remains
There is **one** supported calibration row with failed regression due to a **slope boundary**:

- dataset = SNU
- formula = base
- family = lognormal
- horizon = 9 years

This does **not** look like a global code failure. It looks like a sparse-tail instability / boundary problem confined to one late-horizon model-horizon cell.

### 2) `observed_ipcw_risk` can become misleading when binary support is absent
For **PNU years 3-10**, the exported horizon-reference rows carry:

- `observed_km_risk ≈ 0.299`
- `observed_ipcw_risk = 1.0`
- `nonevent_count_ipcw = 0`

This follows the current implementation, but if these rows are later surfaced in tables or interpretation, the `observed_ipcw_risk = 1.0` value can be misunderstood. A safer rule would be to set `observed_ipcw_risk = NA` whenever binary outcome support is unavailable.

## Recommendation level

- **No major rewrite needed now.**
- **One small code refinement is worth making before downstream interpretation**:
  set `observed_ipcw_risk` to `NA` when `binary_outcome_support_flag == 0` (or when the IPCW denominator has no case/control support).
- The isolated SNU year-9 calibration slope-boundary row should be **logged and tolerated**, not treated as proof that the whole Stage 5 core is broken.
