#!/usr/bin/env python3
"""
gen_batch.py
================
Safe generator for APN Magma batch files in Class-22 and Class-30.

Main fixes compared with gen_batch_v3.py:
  * can generate NEW random seeds/slices and update all_seeds.json;
  * does not silently reuse old all_seeds.json unless --reuse is selected;
  * avoids duplicates by normalized free-bit signatures;
  * every generated Magma job writes to its own result file;
  * every generated Magma job ends with `quit;`, so PowerShell can detect completion;
  * pure ASCII Magma output.

Typical usage on Windows PowerShell:

  python gen_batch.py --new --c22 20 --c30 20
  powershell -ExecutionPolicy Bypass -File .\run_batch_v6.ps1
  python apn_registry.py add-seeds
  python apn_registry.py collect

If you really want to regenerate .m files from the current all_seeds.json without changing seeds:

  python gen_batch.py --reuse

Files generated:
  apn_c22_nl4_s01.m ... apn_c22_nl4_sNN.m
  apn_c30_rnd01.m  ... apn_c30_rndNN.m
  all_seeds.json
"""

import argparse
import json
import os
import random
import time
from pathlib import Path

n = 8
N = 256
pairs = [(p, q) for p in range(n) for q in range(p + 1, n)]
K = len(pairs)
nvars = n * K
pair_idx = {pq: i for i, pq in enumerate(pairs)}

# ---------------------------------------------------------------------------
# Path resolution: gen_batch.py lives in src/, repo root is one level up.
# Seeds file: <repo_root>/seeds/all_seeds.json
# Output .m files: <repo_root>/batch/  (created on demand)
# Results files written by Magma: <repo_root>/batch/results/
# ---------------------------------------------------------------------------
_SCRIPT_DIR = Path(__file__).resolve().parent
_REPO_ROOT   = _SCRIPT_DIR.parent
_DEFAULT_SEEDS_FILE  = _REPO_ROOT / "seeds" / "all_seeds.json"
_DEFAULT_OUT_DIR     = _REPO_ROOT / "batch"

ALL_SEEDS_FILE = _DEFAULT_SEEDS_FILE   # may be overridden by --seeds


def make_companion(poly, n_):
    M = [[0] * n_ for _ in range(n_)]
    for i in range(n_):
        M[i][n_ - 1] = poly[i]
    for i in range(1, n_):
        M[i][i - 1] = 1
    return M


def block_diag(blocks):
    size = sum(len(b) for b in blocks)
    M = [[0] * size for _ in range(size)]
    off = 0
    for b in blocks:
        k = len(b)
        for i in range(k):
            for j in range(k):
                M[off + i][off + j] = b[i][j]
        off += k
    return M


def selfequiv_eqs(A):
    eqs = []
    for j in range(n):
        for r, s in pairs:
            eq = {}
            for p, q in pairs:
                c = (A[p][r] * A[q][s]) ^ (A[p][s] * A[q][r])
                if c:
                    vi = j * K + pair_idx[(p, q)]
                    eq[vi] = eq.get(vi, 0) ^ 1
            for k in range(n):
                if A[j][k]:
                    vi = k * K + pair_idx[(r, s)]
                    eq[vi] = eq.get(vi, 0) ^ 1
            eq = {v: c for v, c in eq.items() if c % 2}
            if eq:
                eqs.append(eq)
    return eqs


def gf2_rref(eqs, nvars_):
    M = []
    for eq in eqs:
        row = [0] * nvars_
        for v, c in eq.items():
            row[v] = c % 2
        M.append(row[:])
    nr = len(M)
    pc = []
    cur = 0
    for col in range(nvars_):
        pivot = next((r for r in range(cur, nr) if M[r][col]), None)
        if pivot is None:
            continue
        M[cur], M[pivot] = M[pivot], M[cur]
        pc.append(col)
        for r in range(nr):
            if r != cur and M[r][col]:
                M[r] = [M[r][k] ^ M[cur][k] for k in range(nvars_)]
        cur += 1
        if cur >= nr:
            break
    fc = [c for c in range(nvars_) if c not in set(pc)]
    sol = {p: {f: M[i][f] for f in fc if M[i][f]} for i, p in enumerate(pc)}
    return pc, fc, sol


# Class definitions
C4 = make_companion([1, 1, 1, 1], 4)
A22 = block_diag([C4, C4])
pc22, fc22, sol22 = gf2_rref(selfequiv_eqs(A22), nvars)
NORM22 = 4
C22_FREE_VARS = [j * K + pair_idx[(p, q)] for (p, q) in pairs if p > NORM22 for j in range(n)]
C22_FIXED_VARS = [j * K + pair_idx[(p, q)] for (p, q) in pairs if p <= NORM22 for j in range(n)]

C2 = make_companion([1, 0], 2)
I2 = [[1, 0], [0, 1]]
A30 = block_diag([I2, C2, C2, C2])
pc30, fc30, sol30 = gf2_rref(selfequiv_eqs(A30), nvars)
N_FIX30 = 100
# In class 30 a slice is defined by fixing the first 100 free variables in fc30.
C30_SLICE_LEN = N_FIX30


def random_bits(k):
    return [random.getrandbits(1) for _ in range(k)]


def c22_signature(coeffs):
    # Deduplicate by FIXED vars (p<=4): these define the slice center.
    # Using free vars (p>4) was incorrect: different centers could share
    # the same free-var bits and produce identical Magma files.
    return "".join(str(coeffs[i]) for i in C22_FIXED_VARS)


def c30_signature(extra100):
    return "".join(str(x) for x in extra100)


def random_c22_coeffs():
    """Return a random point in the Class-22 self-equivalence subspace (dim=40).

    F lies in Class-22 subspace iff F∘A22 = A22∘F. This is a 40-dimensional
    linear subspace of the 224-dimensional quadratic space, parameterised by the
    40 free variables fc22 (RREF of self-equivalence equations).

    The old random_bits(224) generated points OUTSIDE the subspace with probability
    essentially 1, giving an inconsistent APN system -> dim=-1 always.
    """
    coeffs = [0] * nvars
    for f in fc22:
        coeffs[f] = random.getrandbits(1)
    for p, deps in sol22.items():
        val = 0
        for f, c in deps.items():
            val ^= coeffs[f] * c
        coeffs[p] = val % 2
    return coeffs


def random_c30_slice():
    return random_bits(C30_SLICE_LEN)


def load_seed_data():
    if not ALL_SEEDS_FILE.exists():
        return {"meta": {}, "class22_seeds": [], "class30_slices": []}
    with ALL_SEEDS_FILE.open("r", encoding="utf-8") as f:
        data = json.load(f)
    data.setdefault("meta", {})
    data.setdefault("class22_seeds", [])
    data.setdefault("class30_slices", [])
    return data


def save_seed_data(data):
    data.setdefault("meta", {})
    data["meta"]["updated_at"] = time.strftime("%Y-%m-%d %H:%M:%S")
    data["meta"]["generator"] = "gen_batch.py"
    ALL_SEEDS_FILE.parent.mkdir(parents=True, exist_ok=True)
    with ALL_SEEDS_FILE.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)


def extend_seeds(data, c22_count, c30_count):
    old22 = len(data["class22_seeds"])
    old30 = len(data["class30_slices"])
    sig22 = {c22_signature(x) for x in data["class22_seeds"]}
    sig30 = {c30_signature(x) for x in data["class30_slices"]}

    tries = 0
    while len(data["class22_seeds"]) < old22 + c22_count:
        tries += 1
        coeffs = random_c22_coeffs()
        sig = c22_signature(coeffs)
        if sig not in sig22:
            data["class22_seeds"].append(coeffs)
            sig22.add(sig)
        if tries > 100000:
            raise RuntimeError("Too many attempts while generating class-22 seeds")

    tries = 0
    while len(data["class30_slices"]) < old30 + c30_count:
        tries += 1
        sl = random_c30_slice()
        sig = c30_signature(sl)
        if sig not in sig30:
            data["class30_slices"].append(sl)
            sig30.add(sig)
        if tries > 100000:
            raise RuntimeError("Too many attempts while generating class-30 slices")

    return old22, old30, len(data["class22_seeds"]), len(data["class30_slices"])


MAGMA_COMMON = r'''
function DMatrix(a_int, n, pairs, K, sub)
    ab := [(a_int div 2^i) mod 2 : i in [0..n-1]];
    M  := [[R!0 : l in [0..n-1]] : j in [0..n-1]];
    for j := 0 to n-1 do
        jK := j*K;
        for l := 0 to n-1 do
            e := R!0;
            for pp := 0 to l-1 do
                if ab[pp+1] eq 1 then
                    e +:= sub[jK+Index(pairs,[pp,l])];
                end if;
            end for;
            for qq := l+1 to n-1 do
                if ab[qq+1] eq 1 then
                    e +:= sub[jK+Index(pairs,[l,qq])];
                end if;
            end for;
            M[j+1][l+1] := e;
        end for;
    end for;
    return M;
end function;

procedure SaveAPN(results_file, tag, sol_idx, n_total, sb, apn_ok, is_perm)
    sbox_str := "";
    for i := 1 to #sb do
        sbox_str cat:= IntegerToString(sb[i]);
        if i lt #sb then sbox_str cat:= ","; end if;
    end for;
    line := "APN|" cat tag cat "|sol=" cat IntegerToString(sol_idx)
          cat "|total=" cat IntegerToString(n_total)
          cat "|apn=" cat (apn_ok select "true" else "false")
          cat "|perm=" cat (is_perm select "true" else "false")
          cat "|sbox=" cat sbox_str;
    PrintFile(results_file, line);
end procedure;

procedure RunGB(eqs, norm_eqs, results_file, tag, n, N, pairs, K)
    all_eqs := eqs cat norm_eqs;
    t1 := Cputime();
    I := ideal<R | all_eqs>;
    G := GroebnerBasis(I);
    t_gb := Cputime(t1);
    dim_I := Dimension(I);
    print tag, "| GB=", t_gb, "s | dim=", dim_I, "|G|=", #G;
    if dim_I eq 0 then
        V := Variety(I);
        n_sol := #V;
        print tag, "| Solutions:", n_sol;
        for sol_idx := 1 to n_sol do
            sv := V[sol_idx];
            sb := [0 : x in [0..N-1]];
            for x := 0 to N-1 do
                val := 0;
                for j := 0 to n-1 do
                    fj := GF(2)!0;
                    for k := 1 to K do
                        pp := pairs[k][1];
                        qq := pairs[k][2];
                        fj +:= GF(2)!(((x div 2^pp) mod 2)*((x div 2^qq) mod 2))
                             * GF(2)!Integers()!sv[j*K+k];
                    end for;
                    val +:= Integers()!fj * 2^j;
                end for;
                sb[x+1] := val;
            end for;
            apn_ok := true;
            for a := 1 to N-1 do
                cnts := [0 : b in [0..N-1]];
                for x := 0 to N-1 do
                    cnts[BitwiseXor(sb[x+1], sb[BitwiseXor(x,a)+1])+1] +:= 1;
                end for;
                if Max(cnts) gt 2 then
                    apn_ok := false;
                    break;
                end if;
            end for;
            is_perm := #Set(sb) eq N;
            print "  sol", sol_idx, ": APN=", apn_ok, "perm=", is_perm;
            if apn_ok then
                SaveAPN(results_file, tag, sol_idx, n_sol, sb, apn_ok, is_perm);
            end if;
        end for;
    elif dim_I eq -1 then
        print tag, "| dim=-1: no APN in this slice";
    else
        print tag, "| positive dimension or unsupported dimension:", dim_I;
    end if;
end procedure;
'''

MAGMA_APN_BUILD = r'''
m := n-1;
subsets_sorted := [];
for ss in Subsets({1..n}, m) do
    Append(~subsets_sorted, Sort(Setseq(ss)));
end for;
n_ss := #subsets_sorted;

t0 := Cputime();
eqs := [];
for a := 1 to N-1 do
    M := DMatrix(a, n, pairs, K, sub);
    p := R!1;
    for ri := 1 to n_ss do
        rows := subsets_sorted[ri];
        for ci := 1 to n_ss do
            cols := subsets_sorted[ci];
            entries := [M[rows[i]][cols[j]] : i in [1..m], j in [1..m]];
            p *:= (R!1 + Determinant(Matrix(R, m, m, entries)));
        end for;
    end for;
    Append(~eqs, p);
    if a mod 50 eq 0 then
        print "  a=", a, "/255 t=", Cputime(t0), "s";
    end if;
end for;
t_build := Cputime(t0);
print tag, "| build=", t_build, "s | max_terms=", Max([#Terms(e): e in eqs]);
'''


def gen_header(class_name, seed_label, norm_level, free_desc, result_tag):
    return f'''/* auto-generated by gen_batch.py
 * Class: {class_name}  Seed: {seed_label}
 * NL={norm_level} ({free_desc})
 * Results -> apn_results_{result_tag}.txt
 */
n := 8;
N := 2^n;
pairs := [];
for p := 0 to n-1 do
    for q := p+1 to n-1 do
        Append(~pairs, [p,q]);
    end for;
end for;
K := #pairs;
nvars := n*K;
results_file := "apn_results_{result_tag}.txt";
tag := "{result_tag}";
R<[c]> := BooleanPolynomialRing(nvars, "grevlex");
sub := [c[i] : i in [1..nvars]];
norm_eqs := [];
print "=== " cat tag cat " ===";
'''


def gen_class22_nl4(seed_idx, coeffs):
    fixed = [(p, q) for p, q in pairs if p <= NORM22]
    label = f"c22_nl4_s{seed_idx:02d}"
    lines = [gen_header("Class22", f"seed{seed_idx:02d}", 4,
                        "free=24, pairs (5,6)(5,7)(6,7)", label)]
    lines.append("// Fix pairs p<=4 to seed values")
    for p, q in fixed:
        k_ = pair_idx[(p, q)]
        for j in range(n):
            vi = j * K + k_
            v = coeffs[vi]
            if v == 0:
                lines.append(f"sub[{vi+1}] := R!0; Append(~norm_eqs, c[{vi+1}]);")
            else:
                lines.append(f"sub[{vi+1}] := R!1; Append(~norm_eqs, c[{vi+1}]+R!1);")
    lines.append('print "norm_eqs=", #norm_eqs, " free=", nvars-#norm_eqs;')
    lines.append(MAGMA_COMMON)
    lines.append(MAGMA_APN_BUILD)
    lines.append("RunGB(eqs, norm_eqs, results_file, tag, n, N, pairs, K);")
    lines.append(f'print "DONE {label}";')
    lines.append('print "MAGMA_PROCESS_EXIT";')
    lines.append('quit;')
    return "\n".join(lines) + "\n"


def gen_class30_rnd(slice_idx, extra100):
    fixed_extra = {fc30[i]: extra100[i] for i in range(N_FIX30)}
    sym_extra_set = set(fc30[N_FIX30:])
    label = f"c30_rnd{slice_idx:02d}"
    lines = [gen_header("Class30", f"rnd{slice_idx:02d}", 4,
                        "free=16, random slice of fc30", label)]
    lines.append("// Fix 100 of 116 free vars (class-30 random slice)")
    for fc_i, v in sorted(fixed_extra.items()):
        if v == 0:
            lines.append(f"sub[{fc_i+1}] := R!0; Append(~norm_eqs, c[{fc_i+1}]);")
        else:
            lines.append(f"sub[{fc_i+1}] := R!1; Append(~norm_eqs, c[{fc_i+1}]+R!1);")
    lines.append("// Pivot substitutions (class-30 RREF)")
    for pc_i, expr in sol30.items():
        extra_deps = [fc for fc in expr if fc in fixed_extra]
        sym_deps = [fc for fc in expr if fc in sym_extra_set]
        const_val = 0
        for fc in extra_deps:
            const_val ^= fixed_extra[fc]
        if not sym_deps:
            if const_val == 0:
                lines.append(f"sub[{pc_i+1}] := R!0; Append(~norm_eqs, c[{pc_i+1}]);")
            else:
                lines.append(f"sub[{pc_i+1}] := R!1; Append(~norm_eqs, c[{pc_i+1}]+R!1);")
        else:
            sym_part = " + ".join([f"c[{fc+1}]" for fc in sorted(sym_deps)])
            rhs = ("R!1 + " if const_val else "") + sym_part
            lines.append(f"sub[{pc_i+1}] := {rhs};")
            lines.append(f"Append(~norm_eqs, c[{pc_i+1}] + sub[{pc_i+1}]);")
    lines.append('print "norm_eqs=", #norm_eqs, " free=", nvars-#norm_eqs;')
    lines.append(MAGMA_COMMON)
    lines.append(MAGMA_APN_BUILD)
    lines.append("RunGB(eqs, norm_eqs, results_file, tag, n, N, pairs, K);")
    lines.append(f'print "DONE {label}";')
    lines.append('print "MAGMA_PROCESS_EXIT";')
    lines.append('quit;')
    return "\n".join(lines) + "\n"


def write_ascii(path, text):
    clean = "".join(ch if ord(ch) < 128 else "_" for ch in text)
    Path(path).write_text(clean, encoding="ascii")


def remove_old_jobs(out_dir: Path):
    for pattern in ["apn_c22_nl4_s*.m", "apn_c30_rnd*.m"]:
        for p in out_dir.glob(pattern):
            try:
                p.unlink()
            except OSError:
                pass


def generate_magma_files(data, out_dir: Path, only_last=False, old22=0, old30=0):
    out_dir.mkdir(parents=True, exist_ok=True)
    results_dir = out_dir / "results"
    results_dir.mkdir(exist_ok=True)
    remove_old_jobs(out_dir)
    files_c22 = []
    files_c30 = []

    c22_items = data["class22_seeds"]
    c30_items = data["class30_slices"]
    if only_last:
        c22_start = old22
        c30_start = old30
    else:
        c22_start = 0
        c30_start = 0

    for idx, coeffs in enumerate(c22_items[c22_start:], c22_start + 1):
        fname = out_dir / f"apn_c22_nl4_s{idx:02d}.m"
        write_ascii(str(fname), gen_class22_nl4(idx, coeffs))
        files_c22.append(str(fname))

    for idx, extra100 in enumerate(c30_items[c30_start:], c30_start + 1):
        fname = out_dir / f"apn_c30_rnd{idx:02d}.m"
        write_ascii(str(fname), gen_class30_rnd(idx, extra100))
        files_c30.append(str(fname))

    return files_c22, files_c30


def main():
    parser = argparse.ArgumentParser(
        description="Generate Magma APN batch files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Paths (all relative to repo root, auto-detected from script location):
  Seeds file : <repo>/seeds/all_seeds.json   (override with --seeds)
  Output .m  : <repo>/batch/                 (override with --out)
  Results    : <repo>/batch/results/         (Magma writes here)

Examples:
  python src/gen_batch.py --reuse                     # regenerate .m from existing seeds
  python src/gen_batch.py --new --c22 20 --c30 0      # add 20 new C22 slices
  python src/gen_batch.py --reuse --out E:\\run\\batch  # custom output folder
""")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--new",   action="store_true",
                      help="append new random seeds to seeds/all_seeds.json and regenerate .m files")
    mode.add_argument("--reuse", action="store_true",
                      help="reuse existing seeds/all_seeds.json and regenerate .m files (default if neither flag given)")
    parser.add_argument("--c22", type=int, default=20,
                        help="number of new Class-22 seeds to add (--new mode, default 20)")
    parser.add_argument("--c30", type=int, default=20,
                        help="number of new Class-30 slices to add (--new mode, default 20)")
    parser.add_argument("--only-new-files", action="store_true",
                        help="generate .m files only for newly appended seeds")
    parser.add_argument("--random-seed", type=int, default=None,
                        help="Python PRNG seed for reproducibility")
    parser.add_argument("--seeds", type=Path, default=None,
                        help=f"path to all_seeds.json (default: {_DEFAULT_SEEDS_FILE})")
    parser.add_argument("--out", type=Path, default=None,
                        help=f"output directory for .m files (default: {_DEFAULT_OUT_DIR})")
    args = parser.parse_args()

    # Apply overrides
    global ALL_SEEDS_FILE
    if args.seeds:
        ALL_SEEDS_FILE = args.seeds.resolve()
    out_dir = args.out.resolve() if args.out else _DEFAULT_OUT_DIR

    if args.random_seed is not None:
        random.seed(args.random_seed)
    else:
        random.seed()

    data = load_seed_data()
    old22 = len(data["class22_seeds"])
    old30 = len(data["class30_slices"])

    if args.new:
        old22, old30, new22, new30 = extend_seeds(data, args.c22, args.c30)
        save_seed_data(data)
        print(f"Updated {ALL_SEEDS_FILE}: C22 {old22} -> {new22}, C30 {old30} -> {new30}")
    else:
        # --reuse or default
        print(f"Seeds: {ALL_SEEDS_FILE}  (C22={old22}, C30={old30})")

    if old22 == 0 and old30 == 0:
        print("ERROR: seeds file is empty or missing.")
        print(f"  Expected: {ALL_SEEDS_FILE}")
        print("  Run with --new to generate seeds, or check the path.")
        return

    files_c22, files_c30 = generate_magma_files(
        data, out_dir=out_dir,
        only_last=args.only_new_files, old22=old22, old30=old30)

    print(f"Output directory : {out_dir}")
    print(f"Results directory: {out_dir / 'results'}")
    print(f"Class 22 files   : {len(files_c22)}")
    print(f"Class 30 files   : {len(files_c30)}")
    print(f"Total .m files   : {len(files_c22) + len(files_c30)}")

    bad_files = [f for f in files_c22 + files_c30
                 if any(b > 127 for b in Path(f).read_bytes())]
    print(f"Non-ASCII files  : {len(bad_files)}" + (" (OK)" if not bad_files else f" -> {bad_files}"))
    print()
    print("Next step: run Magma from the output directory:")
    print(f"  cd {out_dir}")
    print(r"  powershell -ExecutionPolicy Bypass -File ..\magma_templates\run_batch.ps1")
    print("Then collect results:")
    print(f"  python src/collect_results.py --results {out_dir / 'results'}")


if __name__ == "__main__":
    main()
