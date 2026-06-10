#!/usr/bin/env bash
# run_all_tests.sh
# Run all Magma verification tests sequentially.
#
# Usage:
#   bash tests/run_all_tests.sh                    # uses 'magma' from PATH
#   bash tests/run_all_tests.sh /path/to/magma     # explicit magma binary
#
# Each test writes its output to tests/results/<testname>.txt
# Exit code: 0 if all passed, 1 if any test failed or timed out.
#
# Estimated runtime: 60–90 minutes total.

set -e

MAGMA=${1:-magma}
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_DIR="${SCRIPT_DIR}/results"
mkdir -p "${RESULTS_DIR}"

echo "============================================================"
echo "APN Gröbner Basis Search — Full Magma Test Suite"
echo "Magma binary: ${MAGMA}"
echo "Tests dir: ${SCRIPT_DIR}"
echo "Results dir: ${RESULTS_DIR}"
echo "Started: $(date)"
echo "============================================================"
echo ""

PASS=0
FAIL=0
TOTAL=0

run_test() {
    local name="$1"
    local script="${SCRIPT_DIR}/${name}.m"
    local result="${RESULTS_DIR}/${name}.txt"
    local timeout_sec="${2:-600}"   # default 10 min timeout

    TOTAL=$((TOTAL + 1))
    echo "--- [${TOTAL}] ${name} ---"

    if [ ! -f "${script}" ]; then
        echo "  SKIP: file not found: ${script}"
        return
    fi

    local start_time=$SECONDS
    if timeout "${timeout_sec}" "${MAGMA}" "${script}" > "${result}" 2>&1; then
        local elapsed=$((SECONDS - start_time))
        # Check for common failure keywords in output
        if grep -qiE "ERROR|Runtime error|Assertion failed" "${result}" 2>/dev/null; then
            echo "  FAIL (Magma error in output) — ${elapsed}s — see ${result}"
            FAIL=$((FAIL + 1))
        else
            echo "  PASS — ${elapsed}s — output: ${result}"
            PASS=$((PASS + 1))
        fi
    else
        local elapsed=$((SECONDS - start_time))
        echo "  FAIL (timeout or non-zero exit after ${elapsed}s) — see ${result}"
        FAIL=$((FAIL + 1))
    fi
    echo ""
}

# ---- Level 0: Python quick check (no Magma) ----
echo "=== Level 0: Python sanity check (no Magma) ==="
if python3 "${SCRIPT_DIR}/load_sboxes.py" 2>&1 | tee "${RESULTS_DIR}/python_sanity.txt"; then
    echo "  PASS"
    PASS=$((PASS + 1))
else
    echo "  FAIL"
    FAIL=$((FAIL + 1))
fi
TOTAL=$((TOTAL + 1))
echo ""

# ---- Level 1: Fast invariant tests (~5–10 min each) ----
echo "=== Level 1: Invariant tests ==="
run_test "ea_sanity_check"             120    # ~2 min
run_test "invariant_test"              600    # ~5 min

# ---- Level 2: CCZ tests (~20–40 min each) ----
echo "=== Level 2: CCZ verification ==="
run_test "ccz_correct_table_test"      600    # ~10 min
run_test "ccz_test_dual_fast"         2400    # ~30 min (pairwise on 6 sboxes)

# ---- Level 3: Additional structural tests ----
echo "=== Level 3: Structural tests ==="
run_test "graph_code_weight_invariants" 600
run_test "graph_structure_test"         600
run_test "image_fiber_test"             600

# ---- Summary ----
echo "============================================================"
echo "TEST SUMMARY"
echo "  Total : ${TOTAL}"
echo "  Pass  : ${PASS}"
echo "  Fail  : ${FAIL}"
echo "Finished: $(date)"
echo "============================================================"

if [ "${FAIL}" -gt 0 ]; then
    echo ""
    echo "Some tests FAILED. Check ${RESULTS_DIR}/ for details."
    exit 1
else
    echo ""
    echo "All tests PASSED."
    exit 0
fi
