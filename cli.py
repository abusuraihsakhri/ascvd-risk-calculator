#!/usr/bin/env python3
"""
CLI for the ACC/AHA 2013 Pooled Cohort Equations ASCVD Risk Calculator.

Usage:
    python cli.py calculate --age 55 --sex male --race white \
        --tc 213 --hdl 50 --sbp 120

    python cli.py calculate --age 62 --sex male --race white \
        --tc 230 --hdl 42 --sbp 148 --bp-meds --smoker \
        --ldl 155 --statin-intensity high

    python cli.py batch -i patients.csv -o results.csv

Stdlib only.
"""

import argparse
import csv
import json
import sys

from ascvd import ten_year_ascvd, lifetime_risk, statin_reduction


def _add_patient_args(parser):
    """Add the common patient-input arguments to a parser."""
    parser.add_argument("--age", type=int, required=True, help="Age (40-79)")
    parser.add_argument(
        "--sex", required=True, choices=["male", "female", "M", "F"],
        help="Biological sex",
    )
    parser.add_argument(
        "--race", required=True,
        choices=["white", "african_american", "aa", "black"],
        help="Race (PCE groups: white or african_american)",
    )
    parser.add_argument("--tc", type=float, required=True,
                        help="Total cholesterol (mg/dL, 130-320)")
    parser.add_argument("--hdl", type=float, required=True,
                        help="HDL cholesterol (mg/dL, 20-100)")
    parser.add_argument("--sbp", type=float, required=True,
                        help="Systolic BP (mmHg, 90-200)")
    parser.add_argument("--bp-meds", action="store_true",
                        help="On antihypertensive medication")
    parser.add_argument("--smoker", action="store_true",
                        help="Current smoker")
    parser.add_argument("--diabetes", action="store_true",
                        help="Has diabetes")


def cmd_calculate(args):
    """Run a single ASCVD risk calculation."""
    result = ten_year_ascvd(
        age=args.age,
        sex=args.sex,
        race=args.race,
        total_cholesterol=args.tc,
        hdl=args.hdl,
        sbp=args.sbp,
        on_bp_meds=args.bp_meds,
        smoker=args.smoker,
        diabetes=args.diabetes,
    )

    # ── Print report ─────────────────────────────────────────────────────
    print("=" * 60)
    print("  10-YEAR ASCVD RISK ASSESSMENT")
    print("  ACC/AHA 2013 Pooled Cohort Equations")
    print("=" * 60)
    print(f"  Age:              {args.age}")
    print(f"  Sex:              {args.sex}")
    print(f"  Race:             {args.race}")
    print(f"  Total Cholesterol:{args.tc} mg/dL")
    print(f"  HDL Cholesterol:  {args.hdl} mg/dL")
    print(f"  Systolic BP:      {args.sbp} mmHg")
    print(f"  BP Medication:    {'Yes' if args.bp_meds else 'No'}")
    print(f"  Smoker:           {'Yes' if args.smoker else 'No'}")
    print(f"  Diabetes:         {'Yes' if args.diabetes else 'No'}")
    print("-" * 60)
    print(f"  10-Year Risk:     {result.risk_pct}%")
    print(f"  Risk Category:    {result.category}")
    print(f"  PCE Group:        {result.group}")

    # ── Statin recommendation ────────────────────────────────────────────
    if args.statin_intensity:
        if not args.ldl:
            print("\n  [!] --ldl required for statin reduction estimate")
        else:
            sr = statin_reduction(result.risk_pct, args.ldl, args.statin_intensity)
            print("-" * 60)
            print(f"  Statin Intensity: {args.statin_intensity}")
            print(f"  Current LDL:      {args.ldl} mg/dL")
            print(f"  Est. LDL Reduced: {sr['ldl_reduction_mgdl']} mg/dL")
            print(f"  Risk on Statin:   {sr['reduced_pct']}%")
            print(f"  Absolute Benefit: {sr['absolute_reduction_pct']}%")

    # ── Lifetime risk ────────────────────────────────────────────────────
    if args.lifetime:
        lr = lifetime_risk(
            sex=args.sex,
            total_cholesterol=args.tc,
            hdl=args.hdl,
            sbp=args.sbp,
            on_bp_meds=args.bp_meds,
            smoker=args.smoker,
            diabetes=args.diabetes,
        )
        print("-" * 60)
        print(f"  Lifetime Risk:    {lr['lifetime_risk_pct']}% (to age 85)")
        print(f"  Optimal Factors:  {lr['n_optimal']}/5")
        print(f"  Lifetime Category:{lr['category']}")

    print("=" * 60)

    # ── JSON output ──────────────────────────────────────────────────────
    if args.json:
        out = {
            "risk_pct": result.risk_pct,
            "category": result.category,
            "group": result.group,
        }
        if args.statin_intensity and args.ldl:
            out["statin"] = statin_reduction(
                result.risk_pct, args.ldl, args.statin_intensity
            )
        if args.lifetime:
            out["lifetime"] = lifetime_risk(
                args.sex, args.tc, args.hdl, args.sbp,
                args.bp_meds, args.smoker, args.diabetes,
            )
        print(json.dumps(out, indent=2))

    return 0


def cmd_batch(args):
    """Process a CSV file of patients."""
    with open(args.input, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)

    out_fields = fieldnames + [
        "ascvd_10yr_pct", "ascvd_category", "ascvd_group",
    ]
    out_rows = []

    for row in rows:
        try:
            age = int(row["age"])
            sex = row["sex"]
            race = row["race"]
            tc = float(row["total_cholesterol"])
            hdl_val = float(row["hdl"])
            sbp_val = float(row["sbp"])
            on_meds = row.get("bp_meds", "").strip().lower() in ("1", "true", "yes")
            is_smoker = row.get("smoker", "").strip().lower() in ("1", "true", "yes")
            has_dm = row.get("diabetes", "").strip().lower() in ("1", "true", "yes")

            r = ten_year_ascvd(age, sex, race, tc, hdl_val, sbp_val,
                               on_meds, is_smoker, has_dm)
            merged = dict(row)
            merged["ascvd_10yr_pct"] = str(r.risk_pct)
            merged["ascvd_category"] = r.category
            merged["ascvd_group"] = r.group
        except (ValueError, KeyError) as e:
            merged = dict(row)
            merged["ascvd_10yr_pct"] = f"ERROR: {e}"
            merged["ascvd_category"] = ""
            merged["ascvd_group"] = ""
        out_rows.append(merged)

    with open(args.output, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=out_fields)
        writer.writeheader()
        writer.writerows(out_rows)

    print(f"Processed {len(out_rows)} records -> {args.output}")
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="ascvd",
        description="ACC/AHA 2013 Pooled Cohort Equations — ASCVD Risk Calculator",
    )
    sub = parser.add_subparsers(dest="command")

    # ── calculate ────────────────────────────────────────────────────────
    p_calc = sub.add_parser("calculate", help="Single patient risk calculation")
    _add_patient_args(p_calc)
    p_calc.add_argument("--ldl", type=float, default=None,
                        help="LDL cholesterol (mg/dL) for statin estimate")
    p_calc.add_argument(
        "--statin-intensity", choices=["moderate", "high"], default=None,
        help="Statin intensity for risk reduction estimate",
    )
    p_calc.add_argument("--lifetime", action="store_true",
                        help="Include lifetime risk estimate")
    p_calc.add_argument("--json", action="store_true",
                        help="Also print JSON output")

    # ── batch ────────────────────────────────────────────────────────────
    p_batch = sub.add_parser("batch", help="Batch-process CSV of patients")
    p_batch.add_argument("-i", "--input", required=True, help="Input CSV path")
    p_batch.add_argument("-o", "--output", default="results.csv",
                         help="Output CSV path")

    args = parser.parse_args(argv)

    if args.command == "calculate":
        return cmd_calculate(args)
    if args.command == "batch":
        return cmd_batch(args)

    parser.print_help()
    return 1


if __name__ == "__main__":
    sys.exit(main())
