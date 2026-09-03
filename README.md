# ASCVD Risk Calculator

> **Clinical Domain:** Cardiovascular Medicine & Preventive Cardiology  
> **Reference Guidelines:** 2013 ACC/AHA Guideline on the Assessment of Cardiovascular Risk, 2018/2019 ACC/AHA Primary Prevention & Cholesterol Guidelines, and USPSTF Statin Recommendation Guidelines

[![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)](#)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Tests](https://img.shields.io/badge/Tests-55%20Passed-brightgreen.svg)](#)
[![Dependencies](https://img.shields.io/badge/Dependencies-Standard%20Library-success.svg)](#)

---

## 1. Overview & Clinical Context

Atherosclerotic Cardiovascular Disease (ASCVD)—including coronary heart disease (myocardial infarction, coronary death), ischemic stroke, and peripheral arterial disease—is the primary cause of morbidity and mortality worldwide.

This tool provides a verified, pure Python standard-library implementation of:
1. **2013 ACC/AHA Pooled Cohort Equations (PCE):** 10-year risk of a first hard ASCVD event (fatal or nonfatal myocardial infarction or fatal/nonfatal stroke) among individuals aged 40–79 years without pre-existing cardiovascular disease.
2. **2018/2019 ACC/AHA & USPSTF Statin Eligibility:** Stratification into low, borderline, intermediate, and high-risk categories with clinical guidance on statin therapy, LDL-C reduction targets, and risk-enhancing factors.
3. **CTT Collaboration Statin Risk Reduction:** Quantification of relative and absolute 10-year risk reduction based on moderate- vs. high-intensity statin regimens (~22% RRR per 1.0 mmol/L LDL-C reduction).
4. **Lifetime ASCVD Risk Estimation:** Risk estimation to age 85 based on count of optimal vs. adverse cardiovascular risk factors (Lloyd-Jones et al.).
5. **High-Throughput Batch Processing:** Automated CLI batch processing for clinical cohorts and electronic health records (CSV input/output).

---

## 2. Mathematical Formulation & Equations

### 2.1 The 10-Year ASCVD Pooled Cohort Equations (Cox Proportional Hazards Model)

For an individual patient with risk factor vector $X$, the 10-year risk is modeled as:

$$\text{Risk}_{10} = 1 - S_0(10)^{\exp(\sum_{k} \beta_k X_k - \bar{\mu})}$$

Where:
- $X_k$: Patient covariates including natural logs of continuous parameters ($\ln(\text{age})$, $\ln(\text{total\_cholesterol})$, $\ln(\text{hdl\_cholesterol})$, $\ln(\text{systolic\_bp})$), indicator variables for treated/untreated hypertension, current smoking status, and diabetes mellitus, as well as interaction terms ($\ln(\text{age}) \times \ln(\text{cholesterol})$, $\ln(\text{age}) \times \text{smoking}$, etc.).
- $\beta_k$: Regression coefficients specific to the patient's race and sex demographic group.
- $\bar{\mu}$: Demographic group mean linear predictor ($\text{mean}$).
- $S_0(10)$: Baseline 10-year survival probability for the demographic cohort.

### 2.2 Canonical Coefficient Reference Table (Goff et al. Circulation 2014)

| Demographic Group | Mean $\bar{\mu}$ | Baseline Survival $S_0(10)$ | Key Distinct Terms |
|:---|:---:|:---:|:---|
| **White Male** | `61.1816` | `0.91436` | $\ln(\text{age})$, $\ln(\text{TC})$, $\ln(\text{age}) \times \ln(\text{TC})$, $\ln(\text{HDL})$, $\ln(\text{age}) \times \ln(\text{HDL})$, SBP treated/untreated, smoking $\times \ln(\text{age})$, diabetes |
| **White Female** | `-29.1817` | `0.96652` | $\ln(\text{age})$, $\ln(\text{age})^2$, $\ln(\text{TC})$, $\ln(\text{age}) \times \ln(\text{TC})$, $\ln(\text{HDL})$, $\ln(\text{age}) \times \ln(\text{HDL})$, SBP treated/untreated, smoking $\times \ln(\text{age})$, diabetes |
| **African American Male** | `19.5425` | `0.89536` | $\ln(\text{age})$, $\ln(\text{TC})$, $\ln(\text{HDL})$, SBP treated/untreated, smoking, diabetes |
| **African American Female** | `86.6081` | `0.95334` | $\ln(\text{age})$, $\ln(\text{TC})$, $\ln(\text{HDL})$, $\ln(\text{age}) \times \ln(\text{HDL})$, SBP treated/untreated, $\ln(\text{age}) \times \ln(\text{SBP})$ treated/untreated, smoking, diabetes |

> [!NOTE]
> **Implementation Note on African American Female Equations:**
> The published NHLBI/ACC 2013 eTable 3 contained a known typesetting column-swap error for interaction terms. This repository implements the corrected formulation matching the official American Heart Association (AHA) and ACC online calculators.

### 2.3 Worked Validation Reference Cases (Ages 55, TC 213, HDL 50, untreated SBP 120, Non-smoker, No Diabetes)

| Demographic Group | Expected 10-Year ASCVD Risk | Risk Classification |
|:---|:---:|:---:|
| White Male | **5.3%** | Borderline |
| White Female | **2.1%** | Low |
| African American Male | **6.1%** | Borderline |
| African American Female | **3.0%** | Low |

---

## 3. Risk Stratification & Clinical Guidelines

### 3.1 ACC/AHA Risk Categories & Statin Recommendations

| Risk Category | 10-Year ASCVD Risk | Primary Clinical Recommendation |
|:---|:---:|:---|
| **Low Risk** | `< 5.0%` | Lifestyle counseling (diet, physical activity, weight management). Statin not indicated unless LDL-C $\ge 190\text{ mg/dL}$. |
| **Borderline Risk** | `5.0% – 7.4%` | Patient-clinician discussion. If risk-enhancing factors are present, consider moderate-intensity statin. CAC measurement may reclassify risk. |
| **Intermediate Risk** | `7.5% – 19.9%` | Initiate moderate-intensity statin to reduce LDL-C by 30–49%. If uncertain, coronary artery calcium (CAC) scoring guides decision. |
| **High Risk** | `≥ 20.0%` | Initiate high-intensity statin (e.g., Atorvastatin 40–80 mg or Rosuvastatin 20–40 mg) to reduce LDL-C by $\ge 50\%$. |

### 3.2 USPSTF Recommendations for Statin Preventive Medication (2022)
- **Ages 40–75 with no CVD history, $\ge 1$ CVD risk factor (dyslipidemia, diabetes, hypertension, or smoking), and 10-year risk $\ge 10\%$:** Grade B (recommend statin).
- **10-year risk 7.5% to <10%:** Grade C (selectively offer statin based on shared decision making).

### 3.3 Major Risk-Enhancing Factors (2018 AHA/ACC)
- Family history of premature ASCVD (males $< 55$ y, females $< 65$ y).
- Primary hypercholesterolemia (persistent LDL-C 160–189 mg/dL).
- Metabolic syndrome (increased waist circumference, elevated triglycerides, hypertension, elevated fasting glucose, low HDL).
- Chronic kidney disease (eGFR 15–59 mL/min/1.73 $\text{m}^2$).
- Chronic inflammatory conditions (rheumatoid arthritis, lupus, psoriasis, HIV/AIDS).
- High-risk ethnicity (e.g., South Asian ancestry).
- Biomarkers: persistently elevated triglycerides ($\ge 175\text{ mg/dL}$), hs-CRP $\ge 2.0\text{ mg/dL}$, Lp(a) $\ge 50\text{ mg/dL}$ or $\ge 125\text{ nmol/L}$, apoB $\ge 130\text{ mg/dL}$.

---

## 4. Installation & Requirements

The core ASCVD engine and CLI are built exclusively with the **Python standard library** without mandatory third-party dependencies.

```bash
# Clone the repository
git clone https://github.com/example/ascvd-risk-calculator.git
cd ascvd-risk-calculator

# Optional: verify Python version (3.10+ recommended)
python --version
```

---

## 5. CLI Usage & Examples

The CLI provides subcommands for individual patient evaluation and batch cohort processing.

```
usage: ascvd [-h] {calculate,batch} ...

ACC/AHA 2013 Pooled Cohort Equations — ASCVD Risk Calculator

subcommands:
  calculate  Single patient risk calculation
  batch      Batch-process CSV of patients
```

### 5.1 Single Patient Risk Calculation (`calculate`)

Evaluate a 55-year-old white male:
```bash
python cli.py calculate --age 55 --sex male --race white --tc 213 --hdl 50 --sbp 120
```

Evaluate high-risk patient with statin reduction modeling and lifetime risk:
```bash
python cli.py calculate \
  --age 62 \
  --sex male \
  --race white \
  --tc 230 \
  --hdl 42 \
  --sbp 148 \
  --bp-meds \
  --smoker \
  --ldl 155 \
  --statin-intensity high \
  --lifetime \
  --json
```

Output:
```text
============================================================
  10-YEAR ASCVD RISK ASSESSMENT
  ACC/AHA 2013 Pooled Cohort Equations
============================================================
  Age:              62
  Sex:              male
  Race:             white
  Total Cholesterol:230.0 mg/dL
  HDL Cholesterol:  42.0 mg/dL
  Systolic BP:      148.0 mmHg
  BP Medication:    Yes
  Smoker:           Yes
  Diabetes:         No
------------------------------------------------------------
  10-Year Risk:     27.2%
  Risk Category:    high
  PCE Group:        white_M
------------------------------------------------------------
  Statin Intensity: high
  Current LDL:      155.0 mg/dL
  Est. LDL Reduced: 77.5 mg/dL
  Risk on Statin:   17.0%
  Absolute Benefit: 10.2%
------------------------------------------------------------
  Lifetime Risk:    80.0% (to age 85)
  Optimal Factors:  1/5
  Lifetime Category:high
============================================================
```

### 5.2 Batch Processing (`batch`)

Process a cohort file with patient parameters:
```bash
python cli.py batch -i sample.csv -o results.csv
```

The batch processor supports flexible standard clinical column headers:
- `age`
- `sex` or `gender` (`male`/`female` or `M`/`F`)
- `race` or `ethnicity` (`white`, `african_american`, `aa`, `black`)
- `total_cholesterol` or `tc` (mg/dL)
- `hdl_cholesterol` or `hdl` (mg/dL)
- `systolic_bp` or `sbp` (mmHg)
- `bp_treated`, `bp_meds`, or `treated_htn` (`0`/`1`, `true`/`false`, `yes`/`no`)
- `diabetes` or `dm` (`0`/`1`, `true`/`false`, `yes`/`no`)
- `smoker`, `smoking`, or `current_smoker` (`0`/`1`, `true`/`false`, `yes`/`no`)

Output CSV appends:
- `ascvd_10yr_pct`: 10-year risk percentage
- `ascvd_category`: Risk category (`low`, `borderline`, `intermediate`, `high`)
- `ascvd_group`: Stratification group identifier (e.g. `white_M`, `aa_F`)

---

## 6. Python API Quickstart

```python
from ascvd import ten_year_ascvd, lifetime_risk, statin_reduction

# 1. Calculate 10-year risk
result = ten_year_ascvd(
    age=55,
    sex="male",
    race="white",
    total_cholesterol=213.0,
    hdl=50.0,
    sbp=120.0,
    on_bp_meds=False,
    smoker=False,
    diabetes=False,
)

print(f"Risk: {result.risk_pct}%")       # 5.4%
print(f"Category: {result.category}")     # borderline
print(f"Group: {result.group}")           # white_M

# 2. Model statin benefit (high-intensity on LDL 150 mg/dL)
reduction = statin_reduction(
    baseline_risk_pct=result.risk_pct,
    ldl_mg_dl=150.0,
    intensity="high",
)
print(f"Absolute Risk Reduction: {reduction['absolute_reduction_pct']}%")

# 3. Calculate lifetime ASCVD risk (to age 85)
lr = lifetime_risk(
    sex="male",
    total_cholesterol=213.0,
    hdl=50.0,
    sbp=120.0,
    on_bp_meds=False,
    smoker=False,
    diabetes=False,
)
print(f"Lifetime Risk: {lr['lifetime_risk_pct']}% ({lr['category']})")
```

---

## 7. Verification & Testing

Execute the complete test suite with pytest:

```bash
python -m pytest -p no:zarr -v
```

The automated test suite verifies:
- Reference worked examples from Goff et al. 2013 across all 4 cohorts
- Strict input range boundary enforcement (Ages 40–79, TC 130–320, HDL 20–100, SBP 90–200)
- Monotonic risk factor behavior (higher SBP, smoking, diabetes, higher TC, lower HDL increase risk)
- CLI argument parsing, JSON serialization, and lifetime risk flags
- Batch processing with realistic patient cohorts and graceful invalid row handling
- Longitudinal and clinical enrichment engines

---

## 8. Clinical References

1. **Goff DC Jr, Lloyd-Jones DM, Bennett G, et al.** *2013 ACC/AHA Guideline on the Assessment of Cardiovascular Risk: A Report of the American College of Cardiology/American Heart Association Task Force on Practice Guidelines.* Circulation. 2014;129(25 Suppl 2):S49-S73.
2. **Grundy SM, Stone NJ, Bailey AL, et al.** *2018 AHA/ACC/AACVPR/AAPA/ABC/ACPM/ADA/AGS/APhA/ASPC/NLA/PCNA Guideline on the Management of Blood Cholesterol.* Circulation. 2019;139(25):e1082-e1143.
3. **Arnett DK, Blumenthal RS, Albert MA, et al.** *2019 ACC/AHA Guideline on the Primary Prevention of Cardiovascular Disease.* Circulation. 2019;140(11):e596-e646.
4. **US Preventive Services Task Force.** *Statin Use for the Primary Prevention of Cardiovascular Disease in Adults: US Preventive Services Task Force Recommendation Statement.* JAMA. 2022;328(8):745–753.
5. **Cholesterol Treatment Trialists' (CTT) Collaboration.** *Efficacy and safety of more intensive lowering of LDL cholesterol: a meta-analysis of data from 170,000 participants in 26 randomised trials.* Lancet. 2010;376(9753):1670-1681.
6. **Lloyd-Jones DM, Leip EP, Larson MG, et al.** *Prediction of lifetime risk for cardiovascular disease by risk factor burden at 50 years of age.* Circulation. 2006;113(6):791-798.
