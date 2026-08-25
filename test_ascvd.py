#!/usr/bin/env python3
"""
Tests for the ACC/AHA 2013 Pooled Cohort Equations ASCVD Risk Calculator.

Reference values from Goff et al. Circulation 2014;129(25 Suppl 2):S49-S73
and the AHA online PCE calculator.

Run:  python -m pytest test_ascvd.py -v
  or: python test_ascvd.py
"""

import csv
import math
import sys
import unittest

from ascvd import (
    ten_year_ascvd,
    lifetime_risk,
    statin_reduction,
    _classify_risk,
    _normalize_race,
    _normalize_sex,
    COEFFICIENTS,
)


class TestPCEReferenceCases(unittest.TestCase):
    """
    The four canonical reference cases from the PCE paper:
    55-year-old, TC 213, HDL 50, untreated SBP 120,
    non-smoker, no diabetes.
    """

    def test_white_male(self):
        r = ten_year_ascvd(55, "M", "white", 213, 50, 120, False, False, False)
        self.assertAlmostEqual(r.risk_pct, 5.3, delta=0.2)
        self.assertEqual(r.category, "borderline")
        self.assertEqual(r.group, "white_M")

    def test_white_female(self):
        r = ten_year_ascvd(55, "F", "white", 213, 50, 120, False, False, False)
        self.assertAlmostEqual(r.risk_pct, 2.1, delta=0.2)
        self.assertEqual(r.category, "low")
        self.assertEqual(r.group, "white_F")

    def test_aa_male(self):
        r = ten_year_ascvd(55, "M", "aa", 213, 50, 120, False, False, False)
        self.assertAlmostEqual(r.risk_pct, 6.1, delta=0.2)
        self.assertEqual(r.category, "borderline")
        self.assertEqual(r.group, "aa_M")

    def test_aa_female(self):
        r = ten_year_ascvd(55, "F", "aa", 213, 50, 120, False, False, False)
        self.assertAlmostEqual(r.risk_pct, 3.0, delta=0.2)
        self.assertEqual(r.category, "low")
        self.assertEqual(r.group, "aa_F")


class TestRiskCategories(unittest.TestCase):
    """Verify the four ACC/AHA risk category thresholds."""

    def test_low(self):
        self.assertEqual(_classify_risk(2.0), "low")
        self.assertEqual(_classify_risk(4.9), "low")

    def test_borderline(self):
        self.assertEqual(_classify_risk(5.0), "borderline")
        self.assertEqual(_classify_risk(7.4), "borderline")

    def test_intermediate(self):
        self.assertEqual(_classify_risk(7.5), "intermediate")
        self.assertEqual(_classify_risk(19.9), "intermediate")

    def test_high(self):
        self.assertEqual(_classify_risk(20.0), "high")
        self.assertEqual(_classify_risk(35.0), "high")


class TestInputValidation(unittest.TestCase):
    """Out-of-range inputs must raise ValueError."""

    def test_age_too_low(self):
        with self.assertRaises(ValueError):
            ten_year_ascvd(39, "M", "white", 213, 50, 120, False, False, False)

    def test_age_too_high(self):
        with self.assertRaises(ValueError):
            ten_year_ascvd(80, "M", "white", 213, 50, 120, False, False, False)

    def test_tc_too_low(self):
        with self.assertRaises(ValueError):
            ten_year_ascvd(55, "M", "white", 129, 50, 120, False, False, False)

    def test_tc_too_high(self):
        with self.assertRaises(ValueError):
            ten_year_ascvd(55, "M", "white", 321, 50, 120, False, False, False)

    def test_hdl_too_low(self):
        with self.assertRaises(ValueError):
            ten_year_ascvd(55, "M", "white", 213, 19, 120, False, False, False)

    def test_hdl_too_high(self):
        with self.assertRaises(ValueError):
            ten_year_ascvd(55, "M", "white", 213, 101, 120, False, False, False)

    def test_sbp_too_low(self):
        with self.assertRaises(ValueError):
            ten_year_ascvd(55, "M", "white", 213, 50, 89, False, False, False)

    def test_sbp_too_high(self):
        with self.assertRaises(ValueError):
            ten_year_ascvd(55, "M", "white", 213, 50, 201, False, False, False)

    def test_age_boundary_40(self):
        r = ten_year_ascvd(40, "M", "white", 213, 50, 120, False, False, False)
        self.assertGreaterEqual(r.risk_pct, 0)

    def test_age_boundary_79(self):
        r = ten_year_ascvd(79, "M", "white", 213, 50, 120, False, False, False)
        self.assertLessEqual(r.risk_pct, 100)


class TestRiskFactorEffects(unittest.TestCase):
    """Verify that risk increases with worse risk factors."""

    def test_smoking_increases_risk(self):
        base = ten_year_ascvd(55, "M", "white", 213, 50, 120, False, False, False)
        smk = ten_year_ascvd(55, "M", "white", 213, 50, 120, False, True, False)
        self.assertGreater(smk.risk_pct, base.risk_pct)

    def test_diabetes_increases_risk(self):
        base = ten_year_ascvd(55, "M", "white", 213, 50, 120, False, False, False)
        dm = ten_year_ascvd(55, "M", "white", 213, 50, 120, False, False, True)
        self.assertGreater(dm.risk_pct, base.risk_pct)

    def test_higher_sbp_increases_risk(self):
        low = ten_year_ascvd(55, "M", "white", 213, 50, 120, False, False, False)
        high = ten_year_ascvd(55, "M", "white", 213, 50, 160, False, False, False)
        self.assertGreater(high.risk_pct, low.risk_pct)

    def test_higher_tc_increases_risk(self):
        low = ten_year_ascvd(55, "M", "white", 180, 50, 120, False, False, False)
        high = ten_year_ascvd(55, "M", "white", 280, 50, 120, False, False, False)
        self.assertGreater(high.risk_pct, low.risk_pct)

    def test_lower_hdl_increases_risk(self):
        high_hdl = ten_year_ascvd(55, "M", "white", 213, 70, 120, False, False, False)
        low_hdl = ten_year_ascvd(55, "M", "white", 213, 35, 120, False, False, False)
        self.assertGreater(low_hdl.risk_pct, high_hdl.risk_pct)

    def test_older_age_increases_risk(self):
        young = ten_year_ascvd(45, "M", "white", 213, 50, 120, False, False, False)
        old = ten_year_ascvd(65, "M", "white", 213, 50, 120, False, False, False)
        self.assertGreater(old.risk_pct, young.risk_pct)

    def test_bp_meds_vs_untreated(self):
        """Treated SBP at same level should give slightly different risk."""
        untreated = ten_year_ascvd(55, "M", "white", 213, 50, 140, False, False, False)
        treated = ten_year_ascvd(55, "M", "white", 213, 50, 140, True, False, False)
        # The coefficients differ, so risk should differ
        self.assertNotEqual(untreated.risk_pct, treated.risk_pct)


class TestHighRiskScenario(unittest.TestCase):
    """A 62-year-old male smoker with high lipids and treated hypertension."""

    def test_high_risk_white_male(self):
        r = ten_year_ascvd(62, "M", "white", 230, 42, 148, True, True, False)
        self.assertGreaterEqual(r.risk_pct, 20.0)
        self.assertEqual(r.category, "high")


class TestNormalization(unittest.TestCase):
    """Test input normalization helpers."""

    def test_race_white_variants(self):
        self.assertEqual(_normalize_race("white"), "white")
        self.assertEqual(_normalize_race("White"), "white")
        self.assertEqual(_normalize_race("WHITE"), "white")

    def test_race_aa_variants(self):
        self.assertEqual(_normalize_race("african_american"), "aa")
        self.assertEqual(_normalize_race("African American"), "aa")
        self.assertEqual(_normalize_race("aa"), "aa")
        self.assertEqual(_normalize_race("AA"), "aa")
        self.assertEqual(_normalize_race("black"), "aa")
        self.assertEqual(_normalize_race("Black"), "aa")

    def test_sex_male_variants(self):
        self.assertEqual(_normalize_sex("male"), "M")
        self.assertEqual(_normalize_sex("M"), "M")
        self.assertEqual(_normalize_sex("m"), "M")

    def test_sex_female_variants(self):
        self.assertEqual(_normalize_sex("female"), "F")
        self.assertEqual(_normalize_sex("F"), "F")
        self.assertEqual(_normalize_sex("f"), "F")


class TestCoefficientStructure(unittest.TestCase):
    """Verify all four race-sex groups have complete coefficient sets."""

    EXPECTED_KEYS = {
        "ln_age", "ln_age_sq", "ln_tc", "ln_age_ln_tc",
        "ln_hdl", "ln_age_ln_hdl",
        "ln_sbp_treated", "ln_sbp_untreated",
        "ln_age_ln_sbp_treated", "ln_age_ln_sbp_untreated",
        "smoker", "ln_age_smoker", "diabetes",
        "mean", "s10",
    }

    def test_all_groups_present(self):
        expected_groups = {("white", "M"), ("white", "F"), ("aa", "M"), ("aa", "F")}
        self.assertEqual(set(COEFFICIENTS.keys()), expected_groups)

    def test_all_keys_present(self):
        for group, coeffs in COEFFICIENTS.items():
            with self.subTest(group=group):
                self.assertEqual(set(coeffs.keys()), self.EXPECTED_KEYS)

    def test_s10_in_valid_range(self):
        for group, coeffs in COEFFICIENTS.items():
            with self.subTest(group=group):
                self.assertGreater(coeffs["s10"], 0.0)
                self.assertLess(coeffs["s10"], 1.0)


class TestLifetimeRisk(unittest.TestCase):
    """Test the simplified lifetime risk estimator."""

    def test_optimal_male(self):
        lr = lifetime_risk("M", 170, 60, 110, False, False, False)
        self.assertEqual(lr["n_optimal"], 5)
        self.assertEqual(lr["lifetime_risk_pct"], 5.0)
        self.assertEqual(lr["category"], "optimal")

    def test_all_elevated_male(self):
        lr = lifetime_risk("M", 250, 30, 150, True, True, True)
        self.assertEqual(lr["n_optimal"], 0)
        self.assertEqual(lr["lifetime_risk_pct"], 90.0)
        self.assertEqual(lr["category"], "high")

    def test_optimal_female(self):
        lr = lifetime_risk("F", 170, 60, 110, False, False, False)
        self.assertEqual(lr["n_optimal"], 5)
        self.assertEqual(lr["lifetime_risk_pct"], 4.0)
        self.assertEqual(lr["category"], "optimal")

    def test_partial_optimal(self):
        lr = lifetime_risk("M", 213, 50, 120, False, False, False)
        # TC=213 (not optimal), HDL=50 (optimal), SBP=120 (not <120),
        # non-smoker (optimal), no diabetes (optimal) => 3 optimal
        self.assertEqual(lr["n_optimal"], 3)
        self.assertEqual(lr["lifetime_risk_pct"], 45.0)


class TestStatinReduction(unittest.TestCase):
    """Test statin risk reduction estimation."""

    def test_moderate_statin(self):
        sr = statin_reduction(10.0, 150, "moderate")
        self.assertEqual(sr["baseline_pct"], 10.0)
        self.assertLess(sr["reduced_pct"], 10.0)
        self.assertGreater(sr["absolute_reduction_pct"], 0.0)
        self.assertGreater(sr["relative_reduction_pct"], 0.0)
        self.assertAlmostEqual(sr["ldl_reduction_mgdl"], 52.5, delta=0.1)

    def test_high_statin(self):
        sr = statin_reduction(10.0, 150, "high")
        self.assertLess(sr["reduced_pct"], 10.0)
        self.assertAlmostEqual(sr["ldl_reduction_mgdl"], 75.0, delta=0.1)

    def test_high_statin_greater_benefit(self):
        mod = statin_reduction(10.0, 150, "moderate")
        high = statin_reduction(10.0, 150, "high")
        self.assertGreater(
            high["absolute_reduction_pct"], mod["absolute_reduction_pct"]
        )

    def test_invalid_intensity(self):
        with self.assertRaises(ValueError):
            statin_reduction(10.0, 150, "extreme")

    def test_zero_baseline(self):
        sr = statin_reduction(0.0, 150, "moderate")
        self.assertEqual(sr["reduced_pct"], 0.0)
        self.assertEqual(sr["absolute_reduction_pct"], 0.0)


class TestFormulaIntegrity(unittest.TestCase):
    """Verify the mathematical formula produces bounded, sensible results."""

    def test_risk_bounded_0_to_100(self):
        """Risk should always be between 0 and 100."""
        for age in [40, 55, 79]:
            for sex in ["M", "F"]:
                for race in ["white", "aa"]:
                    r = ten_year_ascvd(age, sex, race, 200, 50, 120, False, False, False)
                    self.assertGreaterEqual(r.risk_pct, 0.0)
                    self.assertLessEqual(r.risk_pct, 100.0)

    def test_extreme_values_still_bounded(self):
        """Even with extreme risk factors, result stays in [0, 100]."""
        r = ten_year_ascvd(79, "M", "aa", 320, 20, 200, True, True, True)
        self.assertGreaterEqual(r.risk_pct, 0.0)
        self.assertLessEqual(r.risk_pct, 100.0)

    def test_result_is_frozen_dataclass(self):
        """AscResult should be immutable."""
        r = ten_year_ascvd(55, "M", "white", 213, 50, 120, False, False, False)
        with self.assertRaises(AttributeError):
            r.risk_pct = 99.9


class TestCLI(unittest.TestCase):
    """Test CLI argument parsing and output (no subprocess, just import)."""

    def test_cli_calculate_runs(self):
        """CLI calculate command should not crash."""
        from cli import main
        # Should return 0 and print output
        ret = main([
            "calculate",
            "--age", "55", "--sex", "male", "--race", "white",
            "--tc", "213", "--hdl", "50", "--sbp", "120",
        ])
        self.assertEqual(ret, 0)

    def test_cli_calculate_with_statin(self):
        from cli import main
        ret = main([
            "calculate",
            "--age", "62", "--sex", "male", "--race", "white",
            "--tc", "230", "--hdl", "42", "--sbp", "148",
            "--bp-meds", "--smoker",
            "--ldl", "155", "--statin-intensity", "high",
        ])
        self.assertEqual(ret, 0)

    def test_cli_calculate_with_lifetime(self):
        from cli import main
        ret = main([
            "calculate",
            "--age", "55", "--sex", "female", "--race", "aa",
            "--tc", "213", "--hdl", "50", "--sbp", "120",
            "--lifetime",
        ])
        self.assertEqual(ret, 0)

    def test_cli_calculate_with_json(self):
        from cli import main
        ret = main([
            "calculate",
            "--age", "55", "--sex", "male", "--race", "white",
            "--tc", "213", "--hdl", "50", "--sbp", "120",
            "--json",
        ])
        self.assertEqual(ret, 0)

    def test_cli_no_command_returns_1(self):
        from cli import main
        ret = main([])
        self.assertEqual(ret, 1)


class TestBatchCSV(unittest.TestCase):
    """Test batch CSV processing."""

    def test_batch_creates_output(self):
        import os
        import tempfile
        from cli import main

        csv_content = (
            "age,sex,race,total_cholesterol,hdl,sbp,bp_meds,smoker,diabetes\n"
            "55,M,white,213,50,120,false,false,false\n"
            "55,F,aa,213,50,120,false,false,false\n"
        )
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".csv", delete=False, newline=""
        ) as inp:
            inp.write(csv_content)
            inp_path = inp.name

        out_path = inp_path.replace(".csv", "_out.csv")

        try:
            ret = main(["batch", "-i", inp_path, "-o", out_path])
            self.assertEqual(ret, 0)
            self.assertTrue(os.path.exists(out_path))

            with open(out_path, newline="") as f:
                reader = csv.DictReader(f)
                rows = list(reader)
            self.assertEqual(len(rows), 2)
            self.assertIn("ascvd_10yr_pct", rows[0])
            self.assertEqual(rows[0]["ascvd_category"], "borderline")
            self.assertEqual(rows[1]["ascvd_category"], "low")
        finally:
            os.unlink(inp_path)
            if os.path.exists(out_path):
                os.unlink(out_path)


if __name__ == "__main__":
    unittest.main()
