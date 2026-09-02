# Ascvd Risk Calculator

> **Domain:** Clinical Decision Support & Biomedical Computing  
> **Reference Guidelines & Standards:** `Standard Clinical Formulations & ISO/IEC Quality Frameworks`

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-0.111-009688.svg?logo=fastapi&logoColor=white)
![Audit Trail](https://img.shields.io/badge/Audit-HMAC--SHA256_Tamper--Evident-brightgreen.svg)
![Zero-PHI Guard](https://img.shields.io/badge/Guard-Zero--PHI_Outbound-blue.svg)
![Docker](https://img.shields.io/badge/Docker-Ready-2496ED.svg?logo=docker&logoColor=white)

</div>

---

## 📖 What It Does

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

ACC/AHA 2013 Pooled Cohort Equations - verified implementation.

Race/sex-specific coefficients from Goff et al. Circulation 2013 / NHLBI
risk-assessment report, validated against the published worked examples
(55 y white male 5.3%, white female 2.1%, African-American female 3.0%,
African-American male 6.1% at TC 213, HDL 50, untreated SBP 120).
Plus: statin eligibility per 2018 AHA/ACC cholesterol guideline and
risk-enhancing factor review.
Stdlib only.

---

## ⚙️ Key Capabilities & Algorithmic Modules

### 🔬 Core Algorithmic & Evaluation Engines

- **`AscResult`**: Result of a 10-year ASCVD risk calculation.
- **`AscvdResult`** — dedicated module for ascvd result evaluation and state verification.

---

## 📐 Mathematical Formulation & Logic

```text
  Formula:
  10-year risk = 1 - S0(10) ^ exp(individual_sum - mean)
  Calculate 10-year ASCVD risk using the 2013 ACC/AHA Pooled Cohort
  risk = 1.0 - math.pow(c["s10"], math.exp(exponent))
```

---

## 💻 CLI Quickstart & Usage

### 1. Guided Interactive Mode
```bash
python cli.py
```

### 2. Direct Parameterized Evaluation
```bash
python cli.py --age <value> --sex <value> --race <value> --tc <value>
```

### Parameter Reference
- `--age`: Specifies input measurement or parameter value.
- `--sex`: Specifies input measurement or parameter value.
- `--race`: Specifies input measurement or parameter value.
- `--tc`: Specifies input measurement or parameter value.
- `--hdl`: Specifies input measurement or parameter value.
- `--sbp`: Specifies input measurement or parameter value.
- `--bp-meds`: Specifies input measurement or parameter value.
- `--smoker`: Specifies input measurement or parameter value.
- `--ldl`: Specifies input measurement or parameter value.
- `--statin-intensity`: Specifies input measurement or parameter value.

### Input Data Schema

| Field | Description | Requirement |
|:------|:------------|:------------|
| `id` | Parameter / observation metric | Required |
| `value` | Parameter / observation metric | Required |
| `qty` | Parameter / observation metric | Required |

---

## 🛡️ Security & Enterprise Architecture

* **Zero-PHI Outbound Interceptor:** Active AST and regex inspection blocking SSNs, MRNs, phone numbers, and patient identifiers.
* **Tamper-Evident HMAC-SHA256 Audit Trail:** Chained, cryptographically signed logs for every evaluation and state transition.
* **Air-Gapped LLM Reasoning Adapter:** Agnostic integration for local Ollama instances (`llama3`, `mistral`), Claude 3.5 Sonnet, GPT-4o, and deterministic test mocks.
* **Active Learning Bayesian Calibration:** Dynamic tracker updating worker reliability weights and monitoring Brier calibration drift.
* **FastAPI & Prometheus Telemetry:** Exposes OpenAPI 3.1 REST endpoints and operational Prometheus metrics (`/metrics`).

---

## 🧪 Testing & Verification

Run the automated test suite:

```bash
pytest -v
```

Execute high-throughput batch simulation benchmarks:

```bash
python simulator.py --tasks 1000 --concurrency 8
```

---

## 🐳 Container Deployment

```bash
docker build -t ascvd-risk-calculator .
docker run -p 8000:8000 ascvd-risk-calculator
```
