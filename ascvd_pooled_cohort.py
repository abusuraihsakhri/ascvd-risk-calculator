#!/usr/bin/env python3
"""
ACC/AHA 2013 Pooled Cohort Equations - verified implementation.

Race/sex-specific coefficients from Goff et al. Circulation 2013 / NHLBI
risk-assessment report, validated against the published worked examples
(55 y white male 5.3%, white female 2.1%, African-American female 3.0%,
African-American male 6.1% at TC 213, HDL 50, untreated SBP 120).
Plus: statin eligibility per 2018 AHA/ACC cholesterol guideline and
risk-enhancing factor review.
Stdlib only.
"""

import math
from dataclasses import dataclass, field


COEFFICIENTS = {
    ("white", "M"): {
        "ln_age": 12.344, "ln_age_sq": 0.0, "ln_tc": 11.853,
        "ln_age_ln_tc": -2.664, "ln_hdl": -7.990, "ln_age_ln_hdl": 1.769,
        "ln_sbp_treated": 1.797, "ln_sbp_untreated": 1.764,
        "smoker": 7.837, "ln_age_smoker": -1.795, "diabetes": 0.658,
        "mean": 61.18, "s10": 0.9144,
    },
    ("white", "F"): {
        "ln_age": -29.799, "ln_age_sq": 4.884, "ln_tc": 13.540,
        "ln_age_ln_tc": -3.114, "ln_hdl": -13.578, "ln_age_ln_hdl": 3.149,
        "ln_sbp_treated": 2.019, "ln_sbp_untreated": 1.957,
        "smoker": 7.574, "ln_age_smoker": -1.665, "diabetes": 0.661,
        "mean": -29.18, "s10": 0.9665,
    },
    ("aa", "M"): {
        "ln_age": 2.469, "ln_age_sq": 0.0, "ln_tc": 0.302,
        "ln_age_ln_tc": 0.0, "ln_hdl": -0.307, "ln_age_ln_hdl": 0.0,
        "ln_sbp_treated": 1.916, "ln_sbp_untreated": 1.809,
        "smoker": 0.549, "ln_age_smoker": 0.0, "diabetes": 0.645,
        "mean": 19.54, "s10": 0.8954,
    },
    ("aa", "F"): {
        "ln_age": 17.114, "ln_age_sq": 0.0, "ln_tc": 0.940,
        "ln_age_ln_tc": 0.0, "ln_hdl": -18.920, "ln_age_ln_hdl": 4.475,
        "ln_sbp_treated": 29.291, "ln_age_ln_sbp_treated": -6.432,
        "ln_sbp_untreated": 27.820, "ln_age_ln_sbp_untreated": -6.087,
        "smoker": 0.691, "ln_age_smoker": 0.0, "diabetes": 0.874,
        "mean": 86.61, "s10": 0.9533,
    },
}


@dataclass
class AscvdResult:
    ten_year_risk_pct: float
    group: str
    individual_sum: float
    exponent: float
    category: str


def ten_year_ascvd(age: int, sex: str, race: str, total_cholesterol: float,
                   hdl: float, sbp: float, treated_htn: bool,
                   smoker: bool, diabetes: bool) -> AscvdResult:
    if not (40 <= age <= 79):
        raise ValueError("PCE valid for ages 40-79")
    key = ("aa" if str(race).lower().startswith(("afric", "black", "aa")) else "white",
           "F" if sex.upper().startswith("F") else "M")
    c = COEFFICIENTS[key]
    ln_a = math.log(age)
    ln_tc = math.log(total_cholesterol)
    ln_h = math.log(hdl)
    ln_p = math.log(sbp)
    s = 1 if smoker else 0
    d = 1 if diabetes else 0

    terms = (
        c["ln_age"] * ln_a
        + c["ln_age_sq"] * ln_a * ln_a
        + c["ln_tc"] * ln_tc
        + c.get("ln_age_ln_tc", 0.0) * ln_a * ln_tc
        + c["ln_hdl"] * ln_h
        + c["ln_age_ln_hdl"] * ln_a * ln_h
        + (c["ln_sbp_treated"] * ln_p if treated_htn else c["ln_sbp_untreated"] * ln_p)
    )
    if key == ("aa", "F"):
        terms += ((c["ln_age_ln_sbp_treated"] if treated_htn else
                   c["ln_age_ln_sbp_untreated"]) * ln_a * ln_p)
    terms += c["smoker"] * s + c["ln_age_smoker"] * s * ln_a + c["diabetes"] * d

    exponent = terms - c["mean"]
    risk = 1.0 - math.pow(c["s10"], math.exp(exponent))
    risk_pct = round(100 * min(max(risk, 0.0), 1.0), 1)
    if risk_pct >= 20:
        cat = "high (>=20%)"
    elif risk_pct >= 7.5:
        cat = "intermediate-high (7.5-19.9%)"
    elif risk_pct >= 5:
        cat = "borderline (5-7.4%)"
    else:
        cat = "low (<5%)"
    return AscvdResult(risk_pct, f"{key[0]}_{key[1]}", round(terms, 2),
                       round(exponent, 3), cat)


RISK_ENHANCERS = [
    "family_history_premature_ascvd", "persistent_LDL_C_ge_160",
    "CKD_eGFR_15_59", "metabolic_syndrome",
    "inflammatory_disease_RA_SLE_PSO", "premature_menopause",
    "high_risk_ethnicity_South_Asian", "TG_persistent_ge_175",
]


def statin_recommendation(result: AscvdResult, diabetes: bool,
                          ldl_mg_dl: float, enhancers: list) -> dict:
    n_enhancers = len(enhancers)
    if ldl_mg_dl >= 190:
        return {"statin": "high-intensity", "driver": "LDL-C >=190 mg/dL"}
    if diabetes:
        intensity = "high" if (result.ten_year_risk_pct >= 20 or n_enhancers >= 1) else "moderate"
        return {"statin": f"{intensity}-intensity", "driver": "diabetes 40-79 y"}
    r = result.ten_year_risk_pct
    if r >= 20:
        return {"statin": "high-intensity", "driver": "10-y risk >=20%"}
    if r >= 7.5:
        note = "consider CAC score if decision uncertain"
        if n_enhancers >= 1:
            return {"statin": "moderate-to-high-intensity",
                    "driver": "risk >=7.5% plus risk-enhancing factors"}
        return {"statin": "moderate-intensity", "driver": "risk >=7.5%", "note": note}
    if r >= 5:
        if n_enhancers >= 1:
            return {"statin": "moderate-intensity reasonable",
                    "driver": "borderline risk with enhancers; consider CAC"}
        return {"statin": "lifestyle emphasis; consider CAC for borderline risk",
                "driver": "borderline risk 5-7.4%"}
    return {"statin": "lifestyle only", "driver": "low risk <5%"}


if __name__ == "__main__":
    checks = [
        (55, "M", "white", 213, 50, 120, False, False, False, 5.3),
        (55, "F", "white", 213, 50, 120, False, False, False, 2.1),
        (55, "F", "aa", 213, 50, 120, False, False, False, 3.0),
        (55, "M", "aa", 213, 50, 120, False, False, False, 6.1),
    ]
    for args in checks:
        got = ten_year_ascvd(*args[:9])
        status = "PASS" if abs(got.ten_year_risk_pct - args[9]) < 0.15 else "FAIL"
        print(f"{status} expected {args[9]}% got {got.ten_year_risk_pct}% [{got.group}]")
    print()
    high = ten_year_ascvd(62, "M", "white", 230, 42, 148, True, True, False)
    print(high)
    print(statin_recommendation(high, False, 155,
                                ["family_history_premature_ascvd"]))
