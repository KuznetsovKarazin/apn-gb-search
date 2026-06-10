#!/usr/bin/env python3
"""
find_apn_centers.py
===================
Step 1: Fast Python search for APN functions inside V_22.
Step 2: Generate Magma NL=4 batch files centered on each found APN.

This is how the original 25 functions were found:
  - Python random search inside V_22 finds APN centers (~1/11000)
  - Each APN center guarantees GB finds at least itself (dim=0 always)
  - GB also finds neighbors -> potentially new CCZ-classes

Speed: ~1000 APN centers/hour on one CPU core.

Usage:
  python find_apn_centers.py --count 100 --out-dir E:/apn-gb-search/batch/batch_002
  python find_apn_centers.py --count 1000 --out-dir E:/apn-gb-search/batch/batch_003
  python find_apn_centers.py --count 100 --random-seed 42  # reproducible
"""

import argparse
import json
import os
import random
import time
from pathlib import Path

import numpy as np

#    Math                                                                       
n, N = 8, 256
pairs    = [(p, q) for p in range(n) for q in range(p+1, n)]
K        = len(pairs)
nvars    = n * K
pair_idx = {pq: i for i, pq in enumerate(pairs)}
NL       = 4   # fix pairs p <= NL

fixed_pairs = [(p, q) for p, q in pairs if p <= NL]  # 25 pairs * 8 = 200 vars

def make_companion(poly, sz):
    M = [[0]*sz for _ in range(sz)]
    for i in range(sz): M[i][sz-1] = poly[i]
    for i in range(1, sz): M[i][i-1] = 1
    return M

def block_diag(blocks):
    sz = sum(len(b) for b in blocks)
    M = [[0]*sz for _ in range(sz)]
    off = 0
    for b in blocks:
        k = len(b)
        for i in range(k):
            for j in range(k): M[off+i][off+j] = b[i][j]
        off += k
    return M

def selfequiv_eqs(A):
    eqs = []
    for j in range(n):
        for r, s in pairs:
            eq = {}
            for p, q in pairs:
                c = (A[p][r]*A[q][s]) ^ (A[p][s]*A[q][r])
                if c:
                    vi = j*K + pair_idx[(p, q)]
                    eq[vi] = eq.get(vi, 0) ^ 1
            for k in range(n):
                if A[j][k]:
                    vi = k*K + pair_idx[(r, s)]
                    eq[vi] = eq.get(vi, 0) ^ 1
            eq = {v: c for v, c in eq.items() if c % 2}
            if eq: eqs.append(eq)
    return eqs

def gf2_rref(eqs, nv):
    M = []
    for eq in eqs:
        row = [0]*nv
        for v, c in eq.items(): row[v] = c % 2
        M.append(row[:])
    nr = len(M); pc = []; cur = 0
    for col in range(nv):
        piv = next((r for r in range(cur, nr) if M[r][col]), None)
        if piv is None: continue
        M[cur], M[piv] = M[piv], M[cur]; pc.append(col)
        for r in range(nr):
            if r != cur and M[r][col]:
                M[r] = [M[r][k]^M[cur][k] for k in range(nv)]
        cur += 1
        if cur >= nr: break
    fc  = [c for c in range(nv) if c not in set(pc)]
    sol = {p: {f: M[i][f] for f in fc if M[i][f]} for i, p in enumerate(pc)}
    return pc, fc, sol

# Precompute Class-22
C4   = make_companion([1,1,1,1], 4)
A22  = block_diag([C4, C4])
pc22, fc22, sol22 = gf2_rref(selfequiv_eqs(A22), nvars)
nf   = len(fc22)  # 40

# Build fast sbox masks for numpy
def build_masks():
    base = np.zeros(N, dtype=np.int32)
    masks = []
    for fc_i in fc22:
        j_fc = fc_i // K; k_fc = fc_i % K
        p_fc, q_fc = pairs[k_fc]
        m = np.zeros(N, dtype=np.int32)
        for x in range(N):
            if ((x >> p_fc) & 1) and ((x >> q_fc) & 1):
                m[x] ^= (1 << j_fc)
        for pc_i, expr in sol22.items():
            if fc_i in expr:
                j_pc = pc_i // K; k_pc = pc_i % K
                p_pc, q_pc = pairs[k_pc]
                for x in range(N):
                    if ((x >> p_pc) & 1) and ((x >> q_pc) & 1):
                        m[x] ^= (1 << j_pc)
        masks.append(m)
    return base, np.array(masks)

def bits_to_sbox(bits_int, base, masks):
    sb = base.copy()
    for i in range(nf):
        if (bits_int >> i) & 1:
            sb ^= masks[i]
    return sb

def is_apn(sb):
    for a in range(1, N):
        diff = sb ^ sb[np.arange(N, dtype=np.int32) ^ a]
        if np.bincount(diff % N, minlength=N).max() > 2:
            return False
    return True

def bits_to_coeffs(bits_int):
    coeffs = [0] * nvars
    for i, fc in enumerate(fc22):
        coeffs[fc] = (bits_int >> i) & 1
    for pc_i, deps in sol22.items():
        val = 0
        for f, c in deps.items():
            val ^= coeffs[f] * c
        coeffs[pc_i] = val % 2
    return coeffs

def slice_sig(coeffs):
    return "".join(str(coeffs[j*K + pair_idx[(p,q)]])
                   for p, q in fixed_pairs for j in range(n))

#    Magma file generation                                                      
MAGMA_COMMON = r"""
function DMatrix(a_int, n, pairs, K, sub)
    ab := [(a_int div 2^i) mod 2 : i in [0..n-1]];
    M  := [[R!0 : l in [0..n-1]] : j in [0..n-1]];
    for j := 0 to n-1 do
        jK := j*K;
        for l := 0 to n-1 do
            e := R!0;
            for pp := 0 to l-1 do
                if ab[pp+1] eq 1 then e +:= sub[jK+Index(pairs,[pp,l])]; end if;
            end for;
            for qq := l+1 to n-1 do
                if ab[qq+1] eq 1 then e +:= sub[jK+Index(pairs,[l,qq])]; end if;
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
    line := "APN|" cat tag
          cat "|sol=" cat IntegerToString(sol_idx)
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
        V := Variety(I); n_sol := #V;
        print tag, "| Solutions:", n_sol;
        for sol_idx := 1 to n_sol do
            sv := V[sol_idx];
            sb := [0 : x in [0..N-1]];
            for x := 0 to N-1 do
                val := 0;
                for j := 0 to n-1 do
                    fj := GF(2)!0;
                    for k := 1 to K do
                        pp := pairs[k][1]; qq := pairs[k][2];
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
                    cnts[BitwiseXor(sb[x+1],sb[BitwiseXor(x,a)+1])+1] +:= 1;
                end for;
                if Max(cnts) gt 2 then apn_ok := false; break; end if;
            end for;
            is_perm := #Set(sb) eq N;
            print "  sol", sol_idx, ": APN=", apn_ok, "perm=", is_perm;
            SaveAPN(results_file, tag, sol_idx, n_sol, sb, apn_ok, is_perm);
        end for;
    elif dim_I eq -1 then
        print tag, "| dim=-1: no APN in this slice";
    end if;
end procedure;
"""

MAGMA_BUILD = r"""
m := n-1;
subsets_sorted := [];
for ss in Subsets({1..n}, m) do
    Append(~subsets_sorted, Sort(Setseq(ss)));
end for;
n_ss := #subsets_sorted;

t0 := Cputime(); eqs := [];
for a := 1 to N-1 do
    M := DMatrix(a, n, pairs, K, sub); p := R!1;
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
print tag, "| build=", Cputime(t0), "s";
"""

def gen_magma(idx, coeffs, tag, results_subdir="results"):
    lines = []
    lines.append(f"""/* APN center search: Class-22 NL=4
 * Center #{idx:04d}  tag={tag}
 * Center is APN -> GB guaranteed dim=0 (finds at least center itself)
 * Results -> {results_subdir}/apn_results_{tag}.txt
 */
n := 8; N := 2^n;
pairs := [];
for p := 0 to n-1 do
    for q := p+1 to n-1 do Append(~pairs, [p,q]); end for;
end for;
K := #pairs; nvars := n*K;
results_file := "{results_subdir}/apn_results_{tag}.txt";
tag := "{tag}";
R<[c]> := BooleanPolynomialRing(nvars, "grevlex");
sub := [c[i] : i in [1..nvars]];
norm_eqs := [];
print "=== " cat tag cat " ===";""")

    lines.append("// Fix pairs p<=4 to APN center values")
    for p, q in fixed_pairs:
        k_ = pair_idx[(p, q)]
        for j in range(n):
            vi = j*K + k_
            v  = coeffs[vi]
            if v == 0:
                lines.append(f"sub[{vi+1}] := R!0; Append(~norm_eqs, c[{vi+1}]);")
            else:
                lines.append(f"sub[{vi+1}] := R!1; Append(~norm_eqs, c[{vi+1}]+R!1);")

    lines.append('print "norm_eqs=", #norm_eqs, " free=", nvars-#norm_eqs;')
    lines.append(MAGMA_COMMON)
    lines.append(MAGMA_BUILD)
    lines.append("RunGB(eqs, norm_eqs, results_file, tag, n, N, pairs, K);")
    lines.append(f'print "DONE {tag}";')
    lines.append('print "MAGMA_PROCESS_EXIT";')
    lines.append('quit;')
    # NO quit; here - avoids ExitCode=1 issue with run_batch.ps1

    code = "\n".join(lines)
    return "".join(ch if ord(ch) < 128 else "_" for ch in code)

#    Main                                                                       
def main():
    parser = argparse.ArgumentParser(
        description="Find APN centers in V_22 via Python, generate Magma NL=4 files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    parser.add_argument("--count",       type=int,  default=100,
                        help="Number of APN centers to find (default: 100)")
    parser.add_argument("--out-dir",     type=Path, default=Path("batch_apn"),
                        help="Output directory for .m files (default: batch_apn)")
    parser.add_argument("--seeds-file",  type=Path, default=None,
                        help="Save found centers to this JSON file")
    parser.add_argument("--random-seed", type=int,  default=None,
                        help="Python PRNG seed for reproducibility")
    parser.add_argument("--tag-prefix",  default="c22a",
                        help="Tag prefix for Magma files (default: c22a)")
    parser.add_argument("--skip-magma",  action="store_true",
                        help="Only find centers, do not generate .m files")
    args = parser.parse_args()

    if args.random_seed is not None:
        random.seed(args.random_seed)
    else:
        random.seed()

    out_dir = Path(args.out_dir)
    if not args.skip_magma:
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "results").mkdir(exist_ok=True)

    print(f"Building sbox masks for V_22 ...")
    t0 = time.time()
    base, masks = build_masks()
    print(f"  Done in {time.time()-t0:.1f}s")
    print()
    print(f"Searching for {args.count} APN centers in V_22 ...")
    print(f"  Expected density: ~1/11000   ~{args.count*11000:,} trials needed")
    print()

    centers   = []
    seen_sigs = set()
    trials    = 0
    t_start   = time.time()
    last_print = t_start

    while len(centers) < args.count:
        trials += 1
        bits = random.randint(0, (1 << nf) - 1)
        sb   = bits_to_sbox(bits, base, masks)

        if not is_apn(sb):
            continue

        coeffs = bits_to_coeffs(bits)
        sig    = slice_sig(coeffs)
        if sig in seen_sigs:
            continue
        seen_sigs.add(sig)
        centers.append({"idx": len(centers)+1, "bits": bits, "coeffs": coeffs,
                        "sbox": sb.tolist(), "sig": sig})

        now = time.time()
        if now - last_print >= 10 or len(centers) % 50 == 0:
            elapsed = now - t_start
            rate    = len(centers) / elapsed * 3600
            eta     = (args.count - len(centers)) / (len(centers)/elapsed) if len(centers) else 0
            print(f"  Found {len(centers):4d}/{args.count}  trials={trials:,}  "
                  f"rate={rate:.0f}/hr  ETA={eta/60:.1f}min", flush=True)
            last_print = now

    elapsed = time.time() - t_start
    print()
    print(f"Found {len(centers)} APN centers in {elapsed:.1f}s "
          f"({elapsed/len(centers)*1000:.1f}ms/center, {len(centers)/elapsed*3600:.0f}/hr)")
    print(f"Trials: {trials:,}  Density: 1/{trials//len(centers):,}")
    print()

    # Save seeds JSON
    seeds_path = args.seeds_file or (out_dir / "apn_centers.json")
    seeds_path = Path(seeds_path)
    seeds_path.parent.mkdir(parents=True, exist_ok=True)
    with open(seeds_path, "w") as f:
        json.dump({
            "count":       len(centers),
            "random_seed": args.random_seed,
            "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
            "density_1_in": trials // len(centers),
            "centers": [
                {"idx": c["idx"], "bits": c["bits"],
                 "sbox_first8": c["sbox"][:8], "sig": c["sig"]}
                for c in centers
            ]
        }, f, indent=2)
    print(f"Centers saved: {seeds_path}")

    # Generate Magma files
    if not args.skip_magma:
        print(f"Generating {len(centers)} Magma files in {out_dir} ...")
        bad = 0
        for c in centers:
            tag   = f"{args.tag_prefix}_{c['idx']:04d}"
            code  = gen_magma(c["idx"], c["coeffs"], tag)
            fname = out_dir / f"apn_{tag}.m"
            fname.write_text(code, encoding="ascii")
            if any(b > 127 for b in fname.read_bytes()):
                bad += 1
        print(f"Generated {len(centers)} files, non-ASCII issues: {bad}")
        print()
        print("Next steps:")
        print(f"  cd {out_dir}")
        print(f"  powershell -ExecutionPolicy Bypass -File ..\\run_batch.ps1 -MaxParallel 8")
        print(f"  python ..\\collect_results.py")

if __name__ == "__main__":
    main()