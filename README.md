# apn-gb-search

**Gröbner Basis Search for New Quadratic APN Functions in Dimension 8**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![arXiv](https://img.shields.io/badge/arXiv-XXXX.XXXXX-b31b1b.svg)](https://arxiv.org/abs/XXXX.XXXXX)

*Oleksandr Kuznetsov — eCampus University (Italy) & V.N. Karazin Kharkiv National University (Ukraine)*

---

## Key Result

**4 new CCZ-equivalence classes** of quadratic APN functions over F₂⁸, verified absent from the Beierle et al. 2025 database of 3,775,599 known instances. A total of **566 APN functions in 6 CCZ-classes** were found by Gröbner basis search within the Class-22 self-equivalence subspace V₂₂ (dim = 40), covering 428 of 65,536 possible NL=4 hyperplane slices (0.65%).

---

## Background

APN (Almost Perfect Nonlinear) functions over F₂ⁿ provide optimal resistance against differential attacks and are essential for S-box design in block ciphers. The complete classification of quadratic APN functions for n=8 remains open despite decades of effort.

This project searches within the **self-equivalence subspace V_A** defined by the order-5 linear automorphism A = block_diag(C(q), C(q)), where q = X⁴+X³+X²+X+1. This subspace has dimension 40 over F₂, giving an APN density of ≈10⁻⁴ — far higher than the global density ≈10⁻⁵⁵. Prior work (Beierle–Brinkmann–Leander 2021) found no APN functions in this subspace via recursive tree search; we demonstrate that Gröbner basis enumeration finds four new CCZ-classes.

---

## Results

| Class   | Size | Ortho-sig (full)         | In V_A | Status                         |
|---------|-----:|--------------------------|:------:|--------------------------------|
| CLASS-A |  362 | `9d95d9c4a5e2dfdd`       |   ✓    | **New** — absent from Beierle 2025 |
| CLASS-B |   72 | `c5f52f5659346f9f`       |   ✗    | **New** — absent from Beierle 2025 |
| CLASS-C |   36 | `74a30c1a17d45761`       |   ✗    | **New** — absent from Beierle 2025 |
| CLASS-D |   30 | `024ba500f4bac35b`       |   ✗    | **New** — absent from Beierle 2025 |
| CLASS-E |   36 | `11cc72af6d61925d`       |   ✓    | Known: Gold x³ (positive control) |
| CLASS-F |   30 | `d15183d1bdf9a1b4`       |   ✓    | Known: Gold x⁹ (positive control) |

**All 566 functions**: differential uniformity = 2, algebraic degree = 2, non-permutation.

### Slice-type structure (key discovery)

The 428 explored slices decompose into exactly three types:

| Type | Slices | APN/slice | Center class | Neighbors found |
|:----:|-------:|:---------:|:------------:|-----------------|
| I    |    362 |     1     | CLASS-A      | none            |
| II   |     36 |     4     | CLASS-E (Gold x³) | 2×CLASS-B + CLASS-C |
| III  |     30 |     2     | CLASS-F (Gold x⁹) | CLASS-D         |

**CLASS-B, C, D lie entirely outside V_A** — they are inaccessible to any search within V_A and are found exclusively by the Gröbner basis step.

### Comparative experiments

| Experiment | Centers | In V_A? | Slices | New neighbors |
|------------|---------|:-------:|-------:|:-------------:|
| Class-22 (Gold-centered) | APN from V_A | ✓ | 66 | **138** |
| Beierle 2025 dataset | APN, outside V_A | ✗ | 532 | **0** |
| Random functions | random (non-APN) | ✗ | 20 | **0** |

Among 3.8 million known APN functions used as GB-slice centers, **zero** new neighbors were found. Gold functions inside V_A are uniquely productive seeds.

---

## Verification of Novelty

The `verification/` directory contains ortho-derivative signature comparison
reports for all 566 found functions against two reference databases:
- **Beierle et al. 2021** (`apn_8bit.txt`, 12,921 functions): 0 matches
- **Beierle et al. 2025** (`new_apns.txt`, 3,775,599 functions): 0 matches

Gold functions (CLASS-E, CLASS-F) return NO_MATCH against `new_apns.txt`
by design — that database contains only non-Gold functions.

---

## Repository Structure

```
apn-gb-search/
│
├── src/                              # Python pipeline scripts
│   ├── find_apn_centers.py           # Step 1: brute-force APN center search in V_A
│   ├── gen_batch.py                  # Step 2a: generate Magma .m files from centers
│   ├── gen_batch_from_dataset.py     # Step 2a alt: use Beierle dataset as centers
│   ├── collect_results.py            # Step 2b: parse Magma output, deduplicate
│   ├── merge_batches.py              # Step 2c: merge results across batches
│   ├── classify_apn.py               # Step 3: CCZ classification via sboxU
│   ├── check_vs_beierle_db.py        # Step 4: streaming comparison vs Beierle 2025 DB
│   ├── analyze_dataset.py            # Step 5: Class-22 membership test
│   ├── verify_vs_known.py            # Step 5b: pairwise CCZ vs Gold controls
│   ├── exact_pairwise.py             # Step 5c: exact pairwise CCZ test (inter-class)
│   ├── make_representatives.py       # Extract one S-box per class
│   └── apn_registry.py               # Provenance registry management
│
├── magma/
│   └── run_batch.ps1                 # PowerShell parallel queue runner (Windows)
│
├── data/
│   ├── class22_basis.json            # Automorphism A, subspace parameters, search summary
│   ├── apn_all.json                  # All 566 found functions (hash, S-box, batch, class)
│   ├── found_apn_classes_report.json # CCZ-class report with ortho-derivative signatures
│   ├── found_apn_classes_report.txt  # Human-readable class summary
│   ├── classification_report.json    # Full classification output
│   ├── ccz_verify_representatives.m  # Magma: pairwise CCZ test for 6 representatives
│   ├── class22_membership_batch003.json  # V_A membership for batch_003 (174 functions)
│   ├── class22_membership_batch004.csv   # V_A membership for batch_004 (392 functions)
│   ├── gold_test.json                # Gold x³, x⁹, x³³ S-boxes (positive controls)
│   └── seeds/
│       ├── batch_003_centers.json    # 128 APN center seeds used in batch_003
│       └── batch_004_centers.json    # 300 APN center seeds used in batch_004
│
├── tests/
│   ├── ccz_test_dual_fast.m          # Magma: fast dual-graph CCZ test
│   ├── ea_sanity_check.m             # Magma: EA-equivalence sanity check
│   ├── gen_test_magma.py             # Generate test Magma files
│   ├── load_sboxes.py                # Load and verify S-box data
│   └── run_all_tests.sh              # Run full test suite
│
├── comparative/
│   ├── sample_batch_rnd_log.txt      # Sample log: random center → dim=-1 (empty)
│   └── sample_beierle_center_log.txt # Sample log: Beierle center → 1 solution only
│
├── .gitignore
├── CITATION.cff
├── LICENSE                           # Code: MIT; Data: CC BY 4.0
├── README.md
└── requirements.txt
```

---

## How It Works

The search uses a two-phase pipeline:

### Phase 1 — Find APN Centers in V_A (Python, fast)

Random 40-bit vectors are projected onto V_A via the RREF parameterization (`fc22 + sol22`): the 184 pivot variables of M_A are determined by the 40 free variables, giving any point in V_A. Each point is evaluated for the APN property (≈0.6 ms/check, ≈600 APN centers/core-hour).

```bash
# Generate 300 APN centers for a new batch
python src/find_apn_centers.py --count 300 --out-dir batch/batch_005
```

### Phase 2 — Gröbner Basis NL=4 Slice (Magma)

For each APN center F₀, fix the 200 ANF coefficients corresponding to pairs (p,q) with p≤4, leaving 24 free variables (pairs (5,6),(5,7),(6,7)). Build the APN ideal over `BooleanPolynomialRing(224)` (≈8–15 min/slice) and solve (< 0.1 s). All APN functions in the 24-dimensional affine hyperplane are found, including those outside V_A.

```bash
# Generate Magma files from centers
python src/gen_batch.py --centers batch/batch_005/apn_centers.json --out-dir batch/batch_005

# Run batch (Windows PowerShell, 8 parallel cores)
cd batch\batch_005
powershell -ExecutionPolicy Bypass -File ..\..\magma\run_batch.ps1 -MaxParallel 8
```

### Phase 3 — Collect and Classify

```bash
# Parse Magma output, deduplicate
python src/collect_results.py --batch-dir batch/batch_005

# Merge with previous results
python src/merge_batches.py --batches batch/batch_003 batch/batch_004 batch/batch_005 \
    --out data/apn_all.json

# CCZ classification (requires SageMath + sboxU)
conda activate sage
python src/classify_apn.py --found data/apn_all.json

# Check V_A membership
python src/analyze_dataset.py --found data/apn_all.json

# Compare against Beierle 2025 database (requires new_apns.txt)
python src/check_vs_beierle_db.py --found data/apn_all.json --known data/new_apns.txt
```

---

## Requirements

| Component | Version | Notes |
|-----------|---------|-------|
| Python | 3.9+ | Pipeline scripts |
| numpy | any | Fast APN center search |
| Magma | V2.20+ | Gröbner basis; commercial license |
| SageMath + sboxU | latest | CCZ classification |
| PowerShell | 5+ | Batch runner (Windows) |
| Beierle 2025 DB | — | Optional; see below |

### Install sboxU

```bash
conda create -n sage sage python=3.11
conda activate sage
pip install sboxU
```

### Download Beierle 2025 database (optional, for Step 4 only)

Download `new_apns.txt` from [Zenodo record 16752428](https://doi.org/10.5281/zenodo.16752428) and place at `data/new_apns.txt`. This file (~3.7M functions) is required only for database comparison and not included in this repository.

---

## Reproducing the Results

The complete computation from seeds:

```bash
# 1. Have the seeds already:
#    data/seeds/batch_003_centers.json  (128 centers)
#    data/seeds/batch_004_centers.json  (300 centers)

# 2. Generate Magma files
python src/gen_batch.py --centers data/seeds/batch_003_centers.json \
    --out-dir batch/batch_003
python src/gen_batch.py --centers data/seeds/batch_004_centers.json \
    --out-dir batch/batch_004

# 3. Run Magma (requires ~57 CPU-hours on Ryzen 7840HS)
cd batch\batch_003
powershell -ExecutionPolicy Bypass -File ..\..\magma\run_batch.ps1 -MaxParallel 8
cd ..\batch_004
powershell -ExecutionPolicy Bypass -File ..\..\magma\run_batch.ps1 -MaxParallel 8

# 4. Collect and classify
python src/collect_results.py --batch-dir batch/batch_003
python src/collect_results.py --batch-dir batch/batch_004
python src/merge_batches.py --batches batch/batch_003 batch/batch_004 --out data/apn_all.json
conda activate sage
python src/classify_apn.py --found data/apn_all.json

# 5. Verify CCZ pairwise (Magma, fast)
magma data/ccz_verify_representatives.m
# Expected: all 15 inter-class pairs return false
```

Pre-computed results are in `data/` and can be used directly without re-running Magma.

---

## Search Coverage

| Parameter | Value |
|-----------|-------|
| Total NL=4 slices in V_A | 65,536 |
| Slices explored | 428 (0.65%) |
| Slices with ≥1 APN | 428 (100%) |
| Functions in V_A (IN_CLASS22) | 428 |
| Functions outside V_A (OUT_CLASS22) | 138 |
| APN density in V_A | ≈1 per 11,000 |
| CPU-hours (Phase 2) | ≈57 h |
| Hardware | AMD Ryzen 7840HS, 8 cores |
| Magma version | V2.28-9 |

---

## References

1. C. Beierle, M. Brinkmann, G. Leander. Linearly self-equivalent APN permutations in small dimension. *IEEE Trans. Inf. Theory* 67(7):4863–4875, 2021. [doi:10.1109/TIT.2021.3071533](https://doi.org/10.1109/TIT.2021.3071533)

2. C. Beierle, G. Leander. New instances of quadratic APN functions. *IEEE Trans. Inf. Theory* 68(1):670–678, 2022. [doi:10.1109/TIT.2021.3120698](https://doi.org/10.1109/TIT.2021.3120698)

3. C. Beierle, P. Langevin, G. Leander, A. Polujan, S. Rasoolzadeh. Millions of inequivalent quadratic APN functions in eight variables. arXiv:2508.04644, 2025. Dataset: [doi:10.5281/zenodo.16752428](https://doi.org/10.5281/zenodo.16752428)

4. A. Canteaut, A. Couvreur, L. Perrin. Recovering or testing extended-affine equivalence. *IEEE Trans. Inf. Theory* 68(9):6187–6206, 2022. [doi:10.1109/TIT.2022.3166692](https://doi.org/10.1109/TIT.2022.3166692)

5. S. Yoshiara. Equivalences of power APN functions with power or quadratic APN functions. *J. Algebraic Combin.* 44(2):561–585, 2016. [doi:10.1007/s10801-016-0680-z](https://doi.org/10.1007/s10801-016-0680-z)

6. L. Perrin. sboxU: Tools for analysing S-boxes. <https://github.com/lpp-crypto/sboxU>

---

## Citation

If you use this code or data, please cite:

```bibtex
@software{kuznetsov2026apngb,
  author  = {Kuznetsov, Oleksandr},
  title   = {{apn-gb-search}: {G}r\"{o}bner Basis Search for New Quadratic {APN} Functions in Dimension 8},
  year    = {2026},
  url     = {https://github.com/KuznetsovKarazin/apn-gb-search},
  doi     = {10.5281/zenodo.XXXXXXX},
  version = {1.0.0}
}
```

---

## Contact

Oleksandr Kuznetsov  
eCampus University, Italy · V.N. Karazin Kharkiv National University, Ukraine  
oleksandr.kuznetsov@uniecampus.it  
ORCID: [0000-0003-2331-6326](https://orcid.org/0000-0003-2331-6326)
