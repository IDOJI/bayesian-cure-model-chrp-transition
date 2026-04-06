# Latency instability priority diagnostics

## Leave-one-out / leave-last-k-out

- PNU / late_event / rank 4: delta weighted latency median = 34.52 years; delta cure fraction = -0.422; removed = PNU_21070304 | PNU_21072201 | PNU_20101301 | PNU_21020102
- SNU / late_event / rank 5: delta weighted latency median = -33.74 years; delta cure fraction = 0.650; removed = SNU_21093 | SNU_21048 | SNU_21038 | SNU_21026 | SNU_21125
- SNU / late_event / rank 4: delta weighted latency median = -31.19 years; delta cure fraction = 0.558; removed = SNU_21093 | SNU_21048 | SNU_21038 | SNU_21026
- SNU / late_event / rank 1: delta weighted latency median = -27.98 years; delta cure fraction = 0.448; removed = SNU_21093
- SNU / late_event / rank 1: delta weighted latency median = -27.98 years; delta cure fraction = 0.448; removed = SNU_21093
- Merged / late_event / rank 5: delta weighted latency median = -26.67 years; delta cure fraction = 0.572; removed = SNU_21093 | SNU_21048 | SNU_21038 | SNU_21026 | SNU_21125
- Merged / late_event / rank 4: delta weighted latency median = -25.87 years; delta cure fraction = 0.522; removed = SNU_21093 | SNU_21048 | SNU_21038 | SNU_21026
- Merged / late_event / rank 3: delta weighted latency median = -24.92 years; delta cure fraction = 0.477; removed = SNU_21093 | SNU_21048 | SNU_21038
- Merged / late_event / rank 2: delta weighted latency median = -22.22 years; delta cure fraction = 0.373; removed = SNU_21093 | SNU_21048
- Merged (site-adjusted) / late_event / rank 5: delta weighted latency median = -18.39 years; delta cure fraction = 0.487; removed = SNU_21093 | SNU_21048 | SNU_21038 | SNU_21026 | SNU_21125
- Merged / late_event / rank 1: delta weighted latency median = -16.66 years; delta cure fraction = 0.235; removed = SNU_21093
- Merged / late_event / rank 1: delta weighted latency median = -16.66 years; delta cure fraction = 0.235; removed = SNU_21093

## Latency family sensitivity

- PNU / exponential: AIC=66.92, cure fraction=0.630, weighted latency median=0.70 years
- PNU / lognormal: AIC=67.36, cure fraction=0.590, weighted latency median=0.92 years
- PNU / loglogistic: AIC=68.29, cure fraction=0.602, weighted latency median=0.81 years
- PNU / weibull: AIC=68.78, cure fraction=0.615, weighted latency median=0.75 years
- SNU / loglogistic: AIC=306.13, cure fraction=0.034, weighted latency median=24.97 years
- SNU / lognormal: AIC=306.37, cure fraction=0.024, weighted latency median=36.74 years
- SNU / exponential: AIC=306.37, cure fraction=0.600, weighted latency median=4.73 years
- SNU / weibull: AIC=307.01, cure fraction=0.413, weighted latency median=10.44 years
- Merged / lognormal: AIC=370.65, cure fraction=0.084, weighted latency median=28.90 years
- Merged / weibull: AIC=370.96, cure fraction=0.529, weighted latency median=4.70 years
- Merged / loglogistic: AIC=371.05, cure fraction=0.142, weighted latency median=20.01 years
- Merged / exponential: AIC=372.49, cure fraction=0.667, weighted latency median=2.19 years
- Merged (site-adjusted) / lognormal: AIC=366.93, cure fraction=0.173, weighted latency median=20.75 years
- Merged (site-adjusted) / weibull: AIC=367.53, cure fraction=0.217, weighted latency median=12.86 years
- Merged (site-adjusted) / loglogistic: AIC=367.57, cure fraction=0.184, weighted latency median=15.98 years
- Merged (site-adjusted) / exponential: AIC=370.34, cure fraction=0.632, weighted latency median=3.03 years

## Bayesian joint uncertainty

- PNU / anchor_informed: cor(cure, sigma)=-0.550; cor(cure, cohort uncured median)=-0.528; posterior mean uncured median=2.02 y; q975=8.18 y
- PNU / neutral_no_external_info: cor(cure, sigma)=-0.462; cor(cure, cohort uncured median)=-0.516; posterior mean uncured median=4.23 y; q975=15.46 y
- SNU / anchor_informed: cor(cure, sigma)=-0.462; cor(cure, cohort uncured median)=-0.668; posterior mean uncured median=20.34 y; q975=54.95 y
- SNU / neutral_no_external_info: cor(cure, sigma)=-0.398; cor(cure, cohort uncured median)=-0.575; posterior mean uncured median=27.45 y; q975=66.75 y
- Merged / anchor_informed: cor(cure, sigma)=-0.564; cor(cure, cohort uncured median)=-0.719; posterior mean uncured median=15.53 y; q975=41.04 y
- Merged / neutral_no_external_info: cor(cure, sigma)=-0.407; cor(cure, cohort uncured median)=-0.596; posterior mean uncured median=23.27 y; q975=58.26 y
- Merged (site-adjusted) / anchor_informed: cor(cure, sigma)=-0.462; cor(cure, cohort uncured median)=-0.684; posterior mean uncured median=17.18 y; q975=39.89 y
- Merged (site-adjusted) / neutral_no_external_info: cor(cure, sigma)=-0.380; cor(cure, cohort uncured median)=-0.628; posterior mean uncured median=21.15 y; q975=46.96 y


## Supported-horizon overlay

- PNU: frequentist lognormal weighted median = 0.92 y and remains within primary-supported horizon 1 y; Bayesian anchor posterior mean uncured median = 2.02 y (95% interval 0.24 to 8.18 y) and exceeds latest stable horizon 1 y by 1.02 y.
- SNU: frequentist lognormal weighted median = 36.74 y and exceeds latest stable horizon 7 y by 29.74 y; Bayesian anchor posterior mean uncured median = 20.34 y (95% interval 3.57 to 54.95 y) and exceeds latest stable horizon 7 y by 13.34 y.
- Merged: frequentist lognormal weighted median = 28.90 y and exceeds latest stable horizon 6 y by 22.90 y; Bayesian anchor posterior mean uncured median = 15.53 y (95% interval 3.03 to 41.04 y) and exceeds latest stable horizon 6 y by 9.53 y.
- Merged (site-adjusted): frequentist lognormal weighted median = 20.75 y and exceeds latest stable horizon 6 y by 14.75 y; Bayesian anchor posterior mean uncured median = 17.18 y (95% interval 4.24 to 39.89 y) and exceeds latest stable horizon 6 y by 11.18 y.
