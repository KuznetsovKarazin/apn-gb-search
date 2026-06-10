#!/usr/bin/env python3
"""
src/gen_batch_from_dataset.py
==============================
Generate NL=k Gröbner basis search jobs with centers taken from an
external APN dataset (e.g., Beierle 2025 new_apns.txt or apn_8bit.txt).

Idea
----
For any quadratic APN function F, its 224 ANF coefficients define a point
in the full quadratic space. Fixing the ANF coefficients of pairs (p,q)
with p <= NL (the "NL=k slice") leaves (n-k-1 choose 2) * n variables free.
At NL=4 (our standard): 3 pairs free * 8 = 24 free variables.

When F itself is APN, the GB system is guaranteed to have at least one
solution (F itself), and build time drops to ~15s (vs 5-15 min for random
slices without a guaranteed APN). This means we can afford many more slices.

Workflow
--------
1. Read functions from an external dataset file (one sbox per line, [v0,...,v255]).
2. For each function: compute quadratic ANF, extract the NL=k fixed part.
3. Deduplicate: skip functions whose fixed-part signature is already in seeds.
4. Generate Magma .m files (same format as gen_batch.py, writes to batch/).
5. Run with run_batch.ps1 as usual; collect with collect_results.py.

Usage
-----
  # Preview: how many unique slices from first 1000 functions?
  python src/gen_batch_from_dataset.py --file data/new_apns.txt --max-lines 1000 --dry-run

  # Generate 200 slice jobs from apn_8bit.txt
  python src/gen_batch_from_dataset.py --file data/apn_8bit.txt --count 200

  # Generate 500 jobs from the big dataset, starting at line 10000
  python src/gen_batch_from_dataset.py --file data/new_apns.txt --count 500 --skip 10000

  # Generate jobs only from functions with image_size in 168-176 (like our new classes)
  python src/gen_batch_from_dataset.py --file data/new_apns.txt --count 200 --filter-image 168 176

Output: batch/*.m files (ready for run_batch.ps1)
        seeds/all_seeds.json updated with new entries
"""

import argparse
import json
import random
import sys
import time
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths and constants
# ---------------------------------------------------------------------------
_REPO       = Path(__file__).resolve().parent.parent
DATA_DIR    = _REPO / "data"
SEEDS_FILE  = _REPO / "seeds" / "all_seeds.json"
BATCH_DIR   = _REPO / "batch"

N  = 256
n  = 8
NL = 4   # normalization level: fix pairs with p <= NL

pairs    = [(p, q) for p in range(n) for q in range(p + 1, n)]
K        = len(pairs)
nvars    = n * K
pair_idx = {pq: i for i, pq in enumerate(pairs)}

fixed_pairs = [(p, q) for p, q in pairs if p <= NL]
free_pairs  = [(p, q) for p, q in pairs if p >  NL]

# ---------------------------------------------------------------------------
# ANF extraction
# ---------------------------------------------------------------------------

def sbox_to_quadratic_anf(sbox):
    """
    Return the 224-bit ANF coefficient vector (list of 0/1 ints) for a
    quadratic function, ordered as: for j in 0..7, for (p,q) in pairs.
    Only quadratic (weight-2) coefficients are extracted; lower-degree
    terms are ignored (they don't affect the GB search).
    """
    coeffs = [0] * nvars
    for j in range(n):
        f = [(sbox[x] >> j) & 1 for x in range(N)]
        # Möbius transform
        for i in range(n):
            step = 1 << i
            for mask in range(N):
                if mask & step:
                    f[mask] ^= f[mask ^ step]
        # Extract weight-2 monomial coefficients
        for k_idx, (p, q) in enumerate(pairs):
            mono_mask = (1 << p) | (1 << q)
            coeffs[j * K + k_idx] = f[mono_mask]
    return coeffs


def nl4_signature(coeffs):
    """Signature of the fixed part (p<=NL pairs) — used for deduplication."""
    parts = []
    for p, q in fixed_pairs:
        k_idx = pair_idx[(p, q)]
        for j in range(n):
            parts.append(str(coeffs[j * K + k_idx]))
    return "".join(parts)

# ---------------------------------------------------------------------------
# Invariants (cheap, pure Python)
# ---------------------------------------------------------------------------

def image_size(sbox):
    return len(set(sbox))

def algebraic_degree(sbox):
    max_deg = 0
    for bit in range(n):
        f = [(sbox[x] >> bit) & 1 for x in range(N)]
        for i in range(n):
            step = 1 << i
            for mask in range(N):
                if mask & step:
                    f[mask] ^= f[mask ^ step]
        for mask in range(N):
            if f[mask]:
                w = bin(mask).count('1')
                if w > max_deg:
                    max_deg = w
    return max_deg

# ---------------------------------------------------------------------------
# File parsing
# ---------------------------------------------------------------------------

def parse_line(line):
    line = line.strip()
    if not line or not line.startswith('['):
        return None
    try:
        vals = json.loads(line)
        if isinstance(vals, list) and len(vals) == N and all(0 <= v < N for v in vals):
            return vals
    except Exception:
        return None

def iter_sboxes(path, skip=0, max_lines=None):
    count = 0
    with open(path, 'r', encoding='utf-8', errors='replace') as f:
        for line_no, line in enumerate(f, 1):
            if line_no <= skip:
                continue
            if max_lines and count >= max_lines:
                break
            sb = parse_line(line)
            if sb is not None:
                count += 1
                yield line_no, sb

# ---------------------------------------------------------------------------
# Magma generation (same format as gen_batch.py)
# ---------------------------------------------------------------------------

MAGMA_HEADER = '''\
/* auto-generated by gen_batch_from_dataset.py
 * Source: {source}  line {line_no}
 * NL={nl} (free={n_free} variables: pairs {free_pairs_str})
 * Center function: image_size={img_size}  sbox[0:4]={sbox_head}
 * Results -> results/apn_results_{tag}.txt
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
results_file := "results/apn_results_{tag}.txt";
tag := "{tag}";
R<[c]> := BooleanPolynomialRing(nvars, "grevlex");
sub := [c[i] : i in [1..nvars]];
norm_eqs := [];
print "=== " cat tag cat " ===";
'''

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
        print tag, "| positive dimension:", dim_I;
    end if;
end procedure;
'''

MAGMA_BUILD_AND_RUN = r'''
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
RunGB(eqs, norm_eqs, results_file, tag, n, N, pairs, K);
print "DONE " cat tag;
print "MAGMA_PROCESS_EXIT";
quit;
'''


def gen_slice_from_function(coeffs, tag, source_file, line_no, img_size, sbox_head):
    """Generate a Magma .m file for an NL=4 slice centered on a given function."""
    free_str = str(free_pairs).replace("(", "[").replace(")", "]").replace(" ", "")

    header = MAGMA_HEADER.format(
        source=Path(source_file).name,
        line_no=line_no,
        nl=NL,
        n_free=len(free_pairs) * n,
        free_pairs_str=str([(p, q) for p, q in free_pairs]),
        img_size=img_size,
        sbox_head=sbox_head,
        tag=tag,
    )

    # Fix-variable equations: for each fixed pair, fix all 8 j-coefficients
    fix_lines = ["// Fix ANF coefficients for pairs p<={NL} (center = source function)".format(NL=NL)]
    for p, q in fixed_pairs:
        k_idx = pair_idx[(p, q)]
        for j in range(n):
            vi = j * K + k_idx + 1   # 1-indexed for Magma
            v  = coeffs[j * K + k_idx]
            if v == 0:
                fix_lines.append(f"sub[{vi}] := R!0; Append(~norm_eqs, c[{vi}]);")
            else:
                fix_lines.append(f"sub[{vi}] := R!1; Append(~norm_eqs, c[{vi}]+R!1);")

    fix_lines.append('print "norm_eqs=", #norm_eqs, " free=", nvars-#norm_eqs;')

    script = header + "\n".join(fix_lines) + "\n" + MAGMA_COMMON + MAGMA_BUILD_AND_RUN

    # Ensure ASCII-only (Magma requirement)
    return "".join(ch if ord(ch) < 128 else "_" for ch in script)


# ---------------------------------------------------------------------------
# Seed registry update
# ---------------------------------------------------------------------------

def load_seeds():
    if not SEEDS_FILE.exists():
        return {"meta": {}, "class22_seeds": [], "class30_slices": [], "dataset_seeds": []}
    with open(SEEDS_FILE) as f:
        d = json.load(f)
    d.setdefault("dataset_seeds", [])
    return d

def save_seeds(data):
    data.setdefault("meta", {})
    data["meta"]["updated_at"] = time.strftime("%Y-%m-%d %H:%M:%S")
    SEEDS_FILE.parent.mkdir(parents=True, exist_ok=True)
    with open(SEEDS_FILE, "w") as f:
        json.dump(data, f, indent=2)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Generate NL=4 GB search jobs centered on functions from an APN dataset.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Strategy
--------
Each APN function in the dataset becomes the center of an NL=4 search slice.
When the center is itself APN, the GB system has at least one known solution
and runs much faster (~15s vs 5-15 min for random slices). Neighbors found in
the same slice that are NOT in the original dataset are potentially new.

Examples
--------
  # Preview: how many unique slices from first 1000 functions?
  python src/gen_batch_from_dataset.py --file data/new_apns.txt --count 1000 --dry-run

  # Generate 200 jobs from the small dataset (all 8157 functions):
  python src/gen_batch_from_dataset.py --file data/apn_8bit.txt --count 200

  # Generate 500 jobs from big dataset, with image_size filter (our new classes: 170-174):
  python src/gen_batch_from_dataset.py --file data/new_apns.txt --count 500 \\
      --filter-image 168 176

  # After generating, run from batch\\ directory:
  #   cd batch && powershell -ExecutionPolicy Bypass -File ..\\magma_templates\\run_batch.ps1
  #   cd .. && python src\\collect_results.py
""")
    parser.add_argument("--file", required=True, type=Path,
                        help="Dataset file (one sbox per line as [v0,...,v255])")
    parser.add_argument("--count", type=int, default=100,
                        help="Number of slice jobs to generate (default: 100)")
    parser.add_argument("--skip", type=int, default=0,
                        help="Skip first N lines of the dataset (default: 0)")
    parser.add_argument("--filter-image", type=int, nargs=2, metavar=("MIN", "MAX"),
                        default=None,
                        help="Only use functions with image_size in [MIN, MAX]")
    parser.add_argument("--dry-run", action="store_true",
                        help="Count unique slices without writing files")
    parser.add_argument("--tag-prefix", default="ds",
                        help="Prefix for job tags (default: 'ds')")
    parser.add_argument("--out", type=Path, default=BATCH_DIR,
                        help=f"Output directory for .m files (default: {BATCH_DIR})")
    args = parser.parse_args()

    if not args.file.exists():
        print(f"ERROR: file not found: {args.file}")
        print(f"Put the dataset in {DATA_DIR}/")
        sys.exit(1)

    out_dir = args.out.resolve()
    if not args.dry_run:
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "results").mkdir(exist_ok=True)

    # Load existing seed signatures to avoid duplicates
    seed_data = load_seeds()
    existing_sigs = set()
    for entry in seed_data.get("class22_seeds", []):
        # compute signature for existing seeds
        sig = "".join(str(entry[j * K + pair_idx[(p, q)]])
                      for p, q in fixed_pairs for j in range(n))
        existing_sigs.add(sig)
    for entry in seed_data.get("dataset_seeds", []):
        existing_sigs.add(entry.get("sig", ""))

    print(f"Dataset  : {args.file}")
    print(f"Skip     : {args.skip}")
    print(f"Target   : {args.count} jobs")
    if args.filter_image:
        print(f"Filter   : image_size in [{args.filter_image[0]}, {args.filter_image[1]}]")
    print(f"Existing : {len(existing_sigs)} known signatures")
    print()

    generated  = 0
    scanned    = 0
    skipped_dup = 0
    skipped_filter = 0
    new_seed_entries = []
    t0 = time.time()

    for line_no, sbox in iter_sboxes(args.file, skip=args.skip):
        scanned += 1

        # Optional image_size filter
        if args.filter_image:
            img = image_size(sbox)
            if not (args.filter_image[0] <= img <= args.filter_image[1]):
                skipped_filter += 1
                if scanned % 50_000 == 0:
                    print(f"  scanned={scanned:,}  generated={generated}  "
                          f"skipped_filter={skipped_filter:,}  elapsed={time.time()-t0:.0f}s")
                continue
        else:
            img = image_size(sbox)

        # Extract ANF and compute NL=4 signature
        coeffs = sbox_to_quadratic_anf(sbox)
        sig    = nl4_signature(coeffs)

        if sig in existing_sigs:
            skipped_dup += 1
            continue

        existing_sigs.add(sig)
        generated += 1
        tag = f"{args.tag_prefix}_l{line_no:07d}"

        if args.dry_run:
            if generated <= 5 or generated % 20 == 0:
                print(f"  [{generated:4d}] line={line_no:7d}  img={img}  tag={tag}")
        else:
            m_path = out_dir / f"apn_{tag}.m"
            script = gen_slice_from_function(
                coeffs, tag,
                source_file=str(args.file),
                line_no=line_no,
                img_size=img,
                sbox_head=sbox[:4],
            )
            m_path.write_text(script, encoding="ascii")

            new_seed_entries.append({
                "tag": tag,
                "source_file": args.file.name,
                "line_no": line_no,
                "image_size": img,
                "sbox_first8": sbox[:8],
                "sig": sig,
            })

        if generated >= args.count:
            break

    elapsed = time.time() - t0

    print()
    print(f"{'='*55}")
    print(f"SUMMARY")
    print(f"{'='*55}")
    print(f"Lines scanned    : {scanned:,}")
    print(f"Skipped (filter) : {skipped_filter:,}")
    print(f"Skipped (dup sig): {skipped_dup:,}")
    print(f"Jobs generated   : {generated}")
    print(f"Elapsed          : {elapsed:.1f}s  ({scanned/(elapsed+0.001):.0f} lines/s)")

    if args.dry_run:
        print()
        print("DRY RUN — no files written.")
        print(f"To generate: remove --dry-run flag")
    else:
        # Update seeds registry
        seed_data["dataset_seeds"].extend(new_seed_entries)
        save_seeds(seed_data)
        print()
        print(f"Files written to : {out_dir}/")
        print(f"Seeds updated    : {SEEDS_FILE}")
        print()
        print("Next steps:")
        print(f"  cd {out_dir.name}")
        print(r"  powershell -ExecutionPolicy Bypass -File ..\magma_templates\run_batch.ps1")
        print("  cd ..")
        print("  python src\\collect_results.py")

    # Estimate total if running on full dataset
    if scanned > 100 and not args.dry_run:
        rate = scanned / elapsed
        total_lines = 3_775_599  # Beierle 2025
        est_total_sec = total_lines / rate
        est_unique = int(total_lines * generated / scanned)
        print()
        print(f"Extrapolation to full Beierle 2025 dataset ({total_lines:,} lines):")
        print(f"  Estimated unique slices : {est_unique:,}")
        print(f"  Scan time               : {est_total_sec/60:.0f} min")


if __name__ == "__main__":
    main()
