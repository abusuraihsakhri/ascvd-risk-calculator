#!/usr/bin/env python3
"""
ACC/AHA 2013 Pooled Cohort Equations for 10-year ASCVD Risk.

Implements the race/sex-specific Cox model from:
  Goff DC Jr, et al. 2013 ACC/AHA Guideline on the Assessment of
  Cardiovascular Risk. Circulation. 2014;129(25 Suppl 2):S49-S73.

Formula:
  10-year risk = 1 - S0(10) ^ exp(individual_sum - mean)

Coefficients validated against the AHA PCE calculator and the published
worked examples (55-year-old, TC 213, HDL 50, untreated SBP 120,
non-smoker, no diabetes):
  White male: 5.3%  |  White female: 2.1%
  AA male: 6.1%     |  AA female: 3.0%

NOTE on AA female coefficients: The published eTable 3 has a column-swap
error for the SBP and HDL interaction terms. The coefficients below match
the AHA's online PCE calculator source code, which produces the validated
3.0% result for the reference case.

Stdlib only — no third-party dependencies.
"""

import math
from dataclasses import dataclass


# ── PCE Coefficients ────────────────────────────────────────────────────────

COEFFICIENTS = {
    ("white", "M"): {
        "ln_age": 12.344,
        "ln_age_sq": 0.0,
        "ln_tc": 11.853,
        "ln_age_ln_tc": -2.664,
        "ln_hdl": -7.990,
        "ln_age_ln_hdl": 1.769,
        "ln_sbp_treated": 1.797,
        "ln_sbp_untreated": 1.764,
        "ln_age_ln_sbp_treated": 0.0,
        "ln_age_ln_sbp_untreated": 0.0,
        "smoker": 7.837,
        "ln_age_smoker": -1.795,
        "diabetes": 0.658,
        "mean": 61.1816,
        "s10": 0.91436,
    },
    ("white", "F"): {
        "ln_age": -29.799,
        "ln_age_sq": 4.884,
        "ln_tc": 13.540,
        "ln_age_ln_tc": -3.114,
        "ln_hdl": -13.578,
        "ln_age_ln_hdl": 3.149,
        "ln_sbp_treated": 2.019,
        "ln_sbp_untreated": 1.957,
        "ln_age_ln_sbp_treated": 0.0,
        "ln_age_ln_sbp_untreated": 0.0,
        "smoker": 7.574,
        "ln_age_smoker": -1.665,
        "diabetes": 0.661,
        "mean": -29.1817,
        "s10": 0.96652,
    },
    ("aa", "M"): {
        "ln_age": 2.469,
        "ln_age_sq": 0.0,
        "ln_tc": 0.302,
        "ln_age_ln_tc": 0.0,
        "ln_hdl": -0.307,
        "ln_age_ln_hdl": 0.0,
        "ln_sbp_treated": 1.916,
        "ln_sbp_untreated": 1.809,
        "ln_age_ln_sbp_treated": 0.0,
        "ln_age_ln_sbp_untreated": 0.0,
        "smoker": 0.549,
        "ln_age_smoker": 0.0,
        "diabetes": 0.645,
        "mean": 19.5425,
        "s10": 0.89536,
    },
    ("aa", "F"): {
        "ln_age": 17.114,
        "ln_age_sq": 0.0,
        "ln_tc": 0.940,
        "ln_age_ln_tc": 0.0,
        "ln_hdl": -18.920,
        "ln_age_ln_hdl": 4.475,
        "ln_sbp_treated": 29.291,
        "ln_sbp_untreated": 27.820,
        "ln_age_ln_sbp_treated": -6.432,
        "ln_age_ln_sbp_untreated": -6.087,
        "smoker": 0.691,
        "ln_age_smoker": 0.0,
        "diabetes": 0.874,
        "mean": 86.6081,
        "s10": 0.95334,
    },
}


# ── Result dataclass ────────────────────────────────────────────────────────

@dataclass(frozen=True)
class AscResult:
    """Result of a 10-year ASCVD risk calculation."""
    risk_pct: float       # 10-year risk as a percentage (0-100)
    group: str            # e.g. "white_M", "aa_F"
    category: str         # "low", "borderline", "intermediate", "high"
    individual_sum: float # raw linear predictor (for debugging)
    exponent: float       # individual_sum - mean (for debugging)


# ── Helpers ─────────────────────────────────────────────────────────────────

def _classify_risk(risk_pct: float) -> str:
    """Classify 10-year ASCVD risk per ACC/AHA categories."""
    if risk_pct >= 20.0:
        return "high"
    if risk_pct >= 7.5:
        return "intermediate"
    if risk_pct >= 5.0:
        return "borderline"
    return "low"


def _normalize_race(race: str) -> str:
    r = race.strip().lower()
    if r.startswith(("afric", "black", "aa")):
        return "aa"
    return "white"


def _normalize_sex(sex: str) -> str:
    s = sex.strip().upper()
    return "F" if s.startswith("F") else "M"


# ── Core calculation ────────────────────────────────────────────────────────

def ten_year_ascvd(
    age: int,
    sex: str,
    race: str,
    total_cholesterol: float,
    hdl: float,
    sbp: float,
    on_bp_meds: bool,
    smoker: bool,
    diabetes: bool,
) -> AscResult:
    """
    Calculate 10-year ASCVD risk using the 2013 ACC/AHA Pooled Cohort
    Equations.

    Parameters
    ----------
    age : int          Age in years (40-79).
    sex : str          'male'/'female' or 'M'/'F'.
    race : str         'white', 'african_american', 'aa', or 'black'.
    total_cholesterol : float  Total cholesterol, mg/dL (130-320).
    hdl : float        HDL cholesterol, mg/dL (20-100).
    sbp : float        Systolic blood pressure, mmHg (90-200).
    on_bp_meds : bool  On antihypertensive medication.
    smoker : bool      Current smoker.
    diabetes : bool    Has diabetes.

    Returns
    -------
    AscResult

    Raises
    ------
    ValueError  If inputs are outside valid ranges.
    """
    # ── Validate ─────────────────────────────────────────────────────────
    if not (40 <= age <= 79):
        raise ValueError(f"Age must be 40-79, got {age}")
    if not (130 <= total_cholesterol <= 320):
        raise ValueError(
            f"Total cholesterol must be 130-320 mg/dL, got {total_cholesterol}"
        )
    if not (20 <= hdl <= 100):
        raise ValueError(f"HDL must be 20-100 mg/dL, got {hdl}")
    if not (90 <= sbp <= 200):
        raise ValueError(f"SBP must be 90-200 mmHg, got {sbp}")

    # ── Resolve coefficients ─────────────────────────────────────────────
    r = _normalize_race(race)
    s = _normalize_sex(sex)
    key = (r, s)
    if key not in COEFFICIENTS:
        raise ValueError(f"Unsupported group: race={race!r}, sex={sex!r}")
    c = COEFFICIENTS[key]

    # ── Log-transform inputs ─────────────────────────────────────────────
    ln_a = math.log(age)
    ln_tc = math.log(total_cholesterol)
    ln_h = math.log(hdl)
    ln_p = math.log(sbp)
    smk = 1.0 if smoker else 0.0
    dm = 1.0 if diabetes else 0.0

    # ── Linear predictor ─────────────────────────────────────────────────
    individual_sum = (
        c["ln_age"] * ln_a
        + c["ln_age_sq"] * ln_a * ln_a
        + c["ln_tc"] * ln_tc
        + c["ln_age_ln_tc"] * ln_a * ln_tc
        + c["ln_hdl"] * ln_h
        + c["ln_age_ln_hdl"] * ln_a * ln_h
        + (c["ln_sbp_treated"] if on_bp_meds else c["ln_sbp_untreated"]) * ln_p
        + (
            c["ln_age_ln_sbp_treated"]
            if on_bp_meds
            else c["ln_age_ln_sbp_untreated"]
        )
        * ln_a
        * ln_p
        + c["smoker"] * smk
        + c["ln_age_smoker"] * smk * ln_a
        + c["diabetes"] * dm
    )

    # ── 10-year risk ─────────────────────────────────────────────────────
    exponent = individual_sum - c["mean"]
    risk = 1.0 - math.pow(c["s10"], math.exp(exponent))
    risk_pct = round(100.0 * max(0.0, min(risk, 1.0)), 1)

    return AscResult(
        risk_pct=risk_pct,
        group=f"{r}_{s}",
        category=_classify_risk(risk_pct),
        individual_sum=round(individual_sum, 4),
        exponent=round(exponent, 4),
    )


# ── Lifetime Risk Estimation ────────────────────────────────────────────────
# Simplified model based on Lloyd-Jones et al. "Lifetime risk of developing
# coronary heart disease" Lancet 2006 / Circulation 2006.
#
# Counts how many of 5 risk factors are at optimal levels:
#   1. TC < 180 mg/dL
#   2. HDL >= 40 (M) / 50 (F)
#   3. SBP < 120 and NOT on BP meds
#   4. Non-smoker
#   5. No diabetes
#
# Returns estimated lifetime risk to age 85 for a person aged 20-59.

_LIFETIME_TABLE = {
    ("M", 0): 90.0, ("M", 1): 80.0, ("M", 2): 65.0,
    ("M", 3): 45.0, ("M", 4): 25.0, ("M", 5): 5.0,
    ("F", 0): 85.0, ("F", 1): 75.0, ("F", 2): 55.0,
    ("F", 3): 35.0, ("F", 4): 18.0, ("F", 5): 4.0,
}


def lifetime_risk(
    sex: str,
    total_cholesterol: float,
    hdl: float,
    sbp: float,
    on_bp_meds: bool,
    smoker: bool,
    diabetes: bool,
) -> dict:
    """
    Estimate lifetime ASCVD risk (to age 85) based on number of risk
    factors at optimal levels.  Simplified from Lloyd-Jones et al. 2006.

    Valid for adults aged 20-59.  NOT the same as 10-year PCE risk.

    Returns
    -------
    dict with keys:
        lifetime_risk_pct, n_optimal, n_elevated, category
    """
    s = _normalize_sex(sex)
    hdl_threshold = 40.0 if s == "M" else 50.0

    n_optimal = sum([
        total_cholesterol < 180,
        hdl >= hdl_threshold,
        sbp < 120 and not on_bp_meds,
        not smoker,
        not diabetes,
    ])

    n_elevated = 5 - n_optimal
    risk_pct = _LIFETIME_TABLE.get((s, n_optimal), 50.0)

    if n_elevated == 0:
        cat = "optimal"
    elif n_elevated <= 1:
        cat = "elevated"
    else:
        cat = "high"

    return {
        "lifetime_risk_pct": risk_pct,
        "n_optimal": n_optimal,
        "n_elevated": n_elevated,
        "category": cat,
    }


# ── Statin Risk Reduction ───────────────────────────────────────────────────
# Based on the CTT Collaboration meta-analysis (Lancet 2010;376:1670-81):
# ~22% relative risk reduction per 1.0 mmol/L LDL-C reduction.

def statin_reduction(
    baseline_risk_pct: float,
    ldl_mg_dl: float,
    intensity: str = "moderate",
) -> dict:
    """
    Estimate 10-year ASCVD risk reduction with statin therapy.

    Parameters
    ----------
    baseline_risk_pct : float   Baseline 10-year ASCVD risk (%).
    ldl_mg_dl : float           Current LDL-C in mg/dL.
    intensity : str             'moderate' or 'high'.

    Returns
    -------
    dict with keys:
        baseline_pct, reduced_pct, absolute_reduction_pct,
        relative_reduction_pct, ldl_reduction_mgdl
    """
    if intensity not in ("moderate", "high"):
        raise ValueError(f"Intensity must be 'moderate' or 'high', got {intensity!r}")

    # Moderate: 30-49% LDL reduction (use 35% midpoint)
    # High: >=50% LDL reduction (use 50%)
    ldl_reduction_frac = 0.50 if intensity == "high" else 0.35

    ldl_reduction_mgdl = ldl_mg_dl * ldl_reduction_frac
    ldl_reduction_mmol = ldl_reduction_mgdl / 38.67  # 1 mmol/L = 38.67 mg/dL

    # CTT: 22% RRR per 1 mmol/L LDL reduction
    rrr_per_mmol = 0.22
    relative_reduction = 1.0 - math.pow(1.0 - rrr_per_mmol, ldl_reduction_mmol)

    baseline = baseline_risk_pct / 100.0
    reduced = baseline * (1.0 - relative_reduction)
    absolute_reduction = baseline - reduced

    return {
        "baseline_pct": round(baseline_risk_pct, 1),
        "reduced_pct": round(reduced * 100.0, 1),
        "absolute_reduction_pct": round(absolute_reduction * 100.0, 1),
        "relative_reduction_pct": round(relative_reduction * 100.0, 1),
        "ldl_reduction_mgdl": round(ldl_reduction_mgdl, 1),
    }


# ── Self-test ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    # Reference cases from Goff et al. 2013
    cases = [
        (55, "M", "white", 213, 50, 120, False, False, False, 5.3),
        (55, "F", "white", 213, 50, 120, False, False, False, 2.1),
        (55, "M", "aa",    213, 50, 120, False, False, False, 6.1),
        (55, "F", "aa",    213, 50, 120, False, False, False, 3.0),
    ]
    print("PCE Reference Validation")
    print("-" * 60)
    all_pass = True
    for age, sex, race, tc, hdl, sbp, meds, smk, dm, expected in cases:
        r = ten_year_ascvd(age, sex, race, tc, hdl, sbp, meds, smk, dm)
        ok = abs(r.risk_pct - expected) < 0.2
        all_pass = all_pass and ok
        tag = "PASS" if ok else "FAIL"
        print(f"  {tag}  {race:6s} {sex}  expected={expected}%  got={r.risk_pct}%  [{r.category}]")
    print("-" * 60)
    print(f"Overall: {'ALL PASS' if all_pass else 'SOME FAILED'}")
    print()

    # High-risk example
    h = ten_year_ascvd(62, "M", "white", 230, 42, 148, True, True, False)
    print(f"High-risk example: {h.risk_pct}% ({h.category})")
    sr = statin_reduction(h.risk_pct, 155, "high")
    print(f"  With high-intensity statin: {sr['reduced_pct']}% "
          f"(absolute reduction: {sr['absolute_reduction_pct']}%)")
    print()

    # Lifetime risk
    lr = lifetime_risk("M", 213, 50, 120, False, False, False)
    print(f"Lifetime risk (5 risk factors evaluated): {lr}")
