# ASCVD Risk Calculator

10-year and lifetime Atherosclerotic Cardiovascular Disease (ASCVD) risk
estimation using the **2013 ACC/AHA Pooled Cohort Equations**.

Pure Python, stdlib only, no dependencies.

## What This Actually Does

Implements the race/sex-specific Cox proportional hazards model published in:

> Goff DC Jr, Lloyd-Jones DM, Bennett G, et al. 2013 ACC/AHA Guideline on
> the Assessment of Cardiovascular Risk. *Circulation*. 2014;129(25 Suppl
> 2):S49-S73.

The calculator produces:

- **10-year ASCVD risk** — probability of a first heart attack or stroke
  within 10 years, for adults aged 40-79
- **Risk category** — low (<5%), borderline (5-7.4%), intermediate (7.5-19.9%),
  or high (≥20%)
- **Statin benefit estimate** — absolute risk reduction with moderate- or
  high-intensity statin therapy (based on CTT meta-analysis, Lancet 2010)
- **Lifetime risk estimate** — simplified model from Lloyd-Jones et al. 2006

## Limitations

This is a **clinical decision support tool**, not a diagnostic device.

- **Valid age range: 40-79.** The PCE were not developed or validated outside
  this range. Results for ages outside 40-79 are not meaningful.
- **Race categories limited to White and African American.** The PCE were
  developed and validated only in these two populations. Applying them to
  other racial/ethnic groups produces unvalidated estimates.
- **10-year horizon only.** The PCE estimate risk over exactly 10 years.
  They do not predict lifetime risk (the lifetime module is a separate,
  simplified model).
- **Does not account for:** family history, coronary artery calcium score,
  inflammatory markers, kidney function, or other risk-enhancing factors
  described in the 2018 AHA/ACC cholesterol guideline.
- **Not a substitute for clinical judgment.** Treatment decisions should
  incorporate shared decision-making, risk-enhancing factors, and patient
  preferences.

## Quick Start

```bash
# Single calculation
python cli.py calculate \
    --age 55 --sex male --race white \
    --tc 213 --hdl 50 --sbp 120

# With statin benefit estimate
python cli.py calculate \
    --age 62 --sex male --race white \
    --tc 230 --hdl 42 --sbp 148 --bp-meds --smoker \
    --ldl 155 --statin-intensity high

# Include lifetime risk
python cli.py calculate \
    --age 55 --sex female --race african_american \
    --tc 213 --hdl 50 --sbp 120 --lifetime

# Batch CSV processing
python cli.py batch -i patients.csv -o results.csv
```

### CSV Format for Batch Mode

```csv
age,sex,race,total_cholesterol,hdl,sbp,bp_meds,smoker,diabetes
55,M,white,213,50,120,false,false,false
62,F,aa,230,42,148,true,true,false
```

## Python API

```python
from ascvd import ten_year_ascvd, lifetime_risk, statin_reduction

# 10-year risk
result = ten_year_ascvd(
    age=55, sex="male", race="white",
    total_cholesterol=213, hdl=50, sbp=120,
    on_bp_meds=False, smoker=False, diabetes=False,
)
print(result.risk_pct)   # 5.3
print(result.category)   # "borderline"

# Statin benefit
sr = statin_reduction(result.risk_pct, ldl_mg_dl=140, intensity="moderate")
print(sr["reduced_pct"])              # risk on statin
print(sr["absolute_reduction_pct"])   # absolute benefit

# Lifetime risk
lr = lifetime_risk("M", 213, 50, 120, False, False, False)
print(lr["lifetime_risk_pct"])
```

## Input Ranges

| Parameter          | Range         | Unit  |
|:-------------------|:--------------|:------|
| Age                | 40 – 79      | years |
| Total cholesterol  | 130 – 320    | mg/dL |
| HDL cholesterol    | 20 – 100     | mg/dL |
| Systolic BP        | 90 – 200     | mmHg  |

Inputs outside these ranges raise `ValueError`.

## Validation

The implementation is validated against the four canonical reference cases
from the PCE publication (55-year-old, TC 213, HDL 50, untreated SBP 120,
non-smoker, no diabetes):

| Group              | Expected | Calculated |
|:-------------------|:---------|:-----------|
| White male         | 5.3%     | 5.3%       |
| White female       | 2.1%     | 2.1%       |
| African American male | 6.1%  | 6.1%       |
| African American female | 3.0% | 3.0%      |

Run the test suite:

```bash
python -m pytest test_ascvd.py -v
```

## Risk Categories

| 10-Year Risk | Category     | Clinical Consideration |
|:-------------|:-------------|:-----------------------|
| <5%          | Low          | Lifestyle counseling   |
| 5 – 7.4%     | Borderline   | Consider risk-enhancing factors, CAC scoring |
| 7.5 – 19.9%  | Intermediate | Statin discussion recommended; consider CAC |
| ≥20%         | High         | High-intensity statin recommended |

## References

1. Goff DC Jr, et al. 2013 ACC/AHA Guideline on the Assessment of
   Cardiovascular Risk. *Circulation*. 2014;129(25 Suppl 2):S49-S73.
2. Lloyd-Jones DM, et al. Lifetime risk of developing coronary heart
   disease. *Lancet*. 2006;367:59-66.
3. Cholesterol Treatment Trialists' (CTT) Collaboration. Efficacy and
   safety of more intensive lowering of LDL cholesterol. *Lancet*.
   2010;376:1670-1681.
4. Grundy SM, et al. 2018 AHA/ACC/AACVPR/AAPA/ABC/ACPM/ADA/AGS/APhA/ASPC/
   NLA/PCNA Guideline on the Management of Blood Cholesterol.
   *Circulation*. 2019;139:e1082-e1143.

## License

MIT License. See [LICENSE](LICENSE).
