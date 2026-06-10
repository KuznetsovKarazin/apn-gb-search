#!/usr/bin/env python3
"""
tests/gen_test_magma.py
=======================
Generate Magma test scripts from apn_found.json and run them,
saving all results to data/verification/.

Usage:
    # Generate .m files only (then run manually in Magma)
    python tests/gen_test_magma.py --gen

    # Generate and run immediately (requires magma in PATH or --magma flag)
    python tests/gen_test_magma.py --run
    python tests/gen_test_magma.py --run --magma /path/to/magma

    # Generate for a specific class only
    python tests/gen_test_magma.py --gen --class CLASS-004
    python tests/gen_test_magma.py --run --class CLASS-004

    # Show what would be generated without writing
    python tests/gen_test_magma.py --dry-run

Generated files (in tests/generated/):
    invariant_test_all25.m          -- all 25 functions, full invariants
    invariant_test_class004.m       -- CLASS-004 only (17 functions)
    invariant_test_class001.m .. _class003.m -- per-class
    ccz_pairwise_class004.m         -- pairwise CCZ within CLASS-004

Results saved to:
    data/verification/invariant_all25.json
    data/verification/invariant_all25_raw.txt
    data/verification/invariant_class004.json
    data/verification/invariant_class004_raw.txt
    data/verification/ccz_pairwise_class004_raw.txt
"""

import argparse
import json
import os
import subprocess
import sys
import textwrap
from datetime import datetime, timezone
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).parent
REPO_ROOT = SCRIPT_DIR.parent
DATA_DIR = REPO_ROOT / 'data'
VERIF_DIR = DATA_DIR / 'verification'
GEN_DIR = SCRIPT_DIR / 'generated'
FOUND_JSON = DATA_DIR / 'apn_found.json'
CLASS_REPORT = VERIF_DIR / 'found_apn_classes_report.json'

N = 256
MODPOLY = 285

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------

def load_functions(path=FOUND_JSON):
    with open(path, 'r', encoding='utf-8') as f:
        d = json.load(f)
    funcs = d.get('functions', {})
    result = {}
    for h, v in funcs.items():
        fid = v.get('id', h[:8])
        sb = v.get('sbox', [])
        if len(sb) == N:
            result[fid] = {
                'id': fid,
                'tag': v.get('tag'),
                'sbox': [int(x) for x in sb],
                'hash': v.get('hash', h[:16]),
            }
    return dict(sorted(result.items()))


def load_class_map(path=CLASS_REPORT):
    """Returns dict: apn_id -> class_id"""
    if not path.exists():
        return {}
    with open(path, 'r', encoding='utf-8') as f:
        d = json.load(f)
    mapping = {}
    for cls in d.get('classes', []):
        cls_id = cls['class_id']
        for m in cls.get('members', []):
            mapping[m['id']] = cls_id
    return mapping


def filter_by_class(funcs, class_map, class_id):
    return {fid: v for fid, v in funcs.items()
            if class_map.get(fid) == class_id}


# ---------------------------------------------------------------------------
# Magma code generation
# ---------------------------------------------------------------------------

MAGMA_LIBRARY = r"""
n := 8;
N := 2^n;
MODPOLY := 285; // x^8 + x^4 + x^3 + x^2 + 1 = 0x11D
COMPUTE_AUTOCORR := true;

function GF2MulInt(a, b, modpoly)
    aa := Integers()!a; bb := Integers()!b; res := 0;
    while bb gt 0 do
        if (bb mod 2) eq 1 then res := BitwiseXor(res, aa); end if;
        bb := bb div 2; aa := 2*aa;
        if aa ge 256 then aa := BitwiseXor(aa, modpoly); end if;
    end while;
    return res mod 256;
end function;

function GF2PowInt(a, e, modpoly)
    result := 1; base := Integers()!a; ee := Integers()!e;
    while ee gt 0 do
        if (ee mod 2) eq 1 then result := GF2MulInt(result, base, modpoly); end if;
        base := GF2MulInt(base, base, modpoly); ee := ee div 2;
    end while;
    return result;
end function;

function GoldSboxPower3()
    return [ GF2PowInt(x, 3, MODPOLY) : x in [0..255] ];
end function;

function ParityInt(z)
    p := 0; zz := Integers()!z;
    while zz gt 0 do p := BitwiseXor(p, zz mod 2); zz := zz div 2; end while;
    return p;
end function;

function HammingWeightInt(z)
    w := 0; zz := Integers()!z;
    while zz gt 0 do w +:= zz mod 2; zz := zz div 2; end while;
    return w;
end function;

function IsPermutationSbox(S) return #Seqset(S) eq #S; end function;
function ImageSize(S) return #Seqset(S); end function;

function HammingDistanceSboxes(A, B)
    return #[ i : i in [1..#A] | A[i] ne B[i] ];
end function;

function SpectrumFromList(L)
    spec := AssociativeArray(Integers());
    for v in L do
        if IsDefined(spec, v) then spec[v] +:= 1; else spec[v] := 1; end if;
    end for;
    keys := Sort(Setseq(Keys(spec)));
    return [ <k, spec[k]> : k in keys ];
end function;

function DifferentialUniformity(S)
    maxdu := 0;
    for a in [1..N-1] do
        counts := [0 : i in [1..N]];
        for x in [0..N-1] do
            b := BitwiseXor(S[x+1], S[BitwiseXor(x,a)+1]);
            counts[b+1] +:= 1;
        end for;
        ma := Maximum(counts);
        if ma gt maxdu then maxdu := ma; end if;
    end for;
    return maxdu;
end function;

function DifferentialSpectrum(S)
    vals := [];
    for a in [1..N-1] do
        counts := [0 : i in [1..N]];
        for x in [0..N-1] do
            b := BitwiseXor(S[x+1], S[BitwiseXor(x,a)+1]);
            counts[b+1] +:= 1;
        end for;
        vals cat:= counts;
    end for;
    return SpectrumFromList(vals);
end function;

function WalshAbsSpectrum(S)
    vals := [];
    for u in [0..N-1] do
        for v in [1..N-1] do
            sum := 0;
            for x in [0..N-1] do
                parity := ParityInt(BitwiseXor(BitwiseAnd(u,x), BitwiseAnd(v,S[x+1])));
                if parity eq 0 then sum +:= 1; else sum -:= 1; end if;
            end for;
            Append(~vals, Abs(sum));
        end for;
    end for;
    return SpectrumFromList(vals);
end function;

function AutocorrelationAbsSpectrum(S)
    vals := [];
    for a in [1..N-1] do
        for v in [1..N-1] do
            sum := 0;
            for x in [0..N-1] do
                parity := ParityInt(BitwiseAnd(v, BitwiseXor(S[x+1], S[BitwiseXor(x,a)+1])));
                if parity eq 0 then sum +:= 1; else sum -:= 1; end if;
            end for;
            Append(~vals, Abs(sum));
        end for;
    end for;
    return SpectrumFromList(vals);
end function;

function ComponentAlgebraicDegrees(S)
    degs := [];
    for bit in [0..n-1] do
        f := [ (S[x+1] div 2^bit) mod 2 : x in [0..N-1] ];
        for i in [0..n-1] do
            step := 2^i;
            for mask in [0..N-1] do
                if BitwiseAnd(mask, step) ne 0 then
                    f[mask+1] := (f[mask+1] + f[mask-step+1]) mod 2;
                end if;
            end for;
        end for;
        d := -1;
        for mask in [0..N-1] do
            if f[mask+1] ne 0 then
                w := HammingWeightInt(mask);
                if w gt d then d := w; end if;
            end if;
        end for;
        Append(~degs, d);
    end for;
    return degs;
end function;

function FullInvariantRecord(S)
    rec := AssociativeArray();
    rec["perm"]       := IsPermutationSbox(S);
    rec["image"]      := ImageSize(S);
    rec["DU"]         := DifferentialUniformity(S);
    rec["DDT"]        := DifferentialSpectrum(S);
    rec["WalshAbs"]   := WalshAbsSpectrum(S);
    rec["coord_degs"] := ComponentAlgebraicDegrees(S);
    rec["alg_deg"]    := Maximum(rec["coord_degs"]);
    if COMPUTE_AUTOCORR then
        rec["ACAbs"]  := AutocorrelationAbsSpectrum(S);
    end if;
    return rec;
end function;

procedure PrintRecord(label, rec, gold_rec, outfile)
    PrintFile(outfile, "-------------------------------------------------------");
    PrintFile(outfile, label);
    PrintFile(outfile, "  perm      : " cat (rec["perm"] select "true" else "false"));
    PrintFile(outfile, "  image_size: " cat IntegerToString(rec["image"]));
    PrintFile(outfile, "  DU        : " cat IntegerToString(rec["DU"]));
    PrintFile(outfile, "  alg_degree: " cat IntegerToString(rec["alg_deg"]));
    s := ""; for d in rec["coord_degs"] do s cat:= IntegerToString(d) cat " "; end for;
    PrintFile(outfile, "  coord_degs: " cat s);
    PrintFile(outfile, "  DDT       : " cat Sprint(rec["DDT"]));
    PrintFile(outfile, "  WalshAbs  : " cat Sprint(rec["WalshAbs"]));
    if COMPUTE_AUTOCORR then
        PrintFile(outfile, "  ACAbs     : " cat Sprint(rec["ACAbs"]));
    end if;
    PrintFile(outfile, "  same_Walsh_as_gold: " cat
              (rec["WalshAbs"] eq gold_rec["WalshAbs"] select "true" else "false"));
    PrintFile(outfile, "  same_image_as_gold: " cat
              (rec["image"] eq gold_rec["image"] select "true" else "false"));
end procedure;

GOLD_SBOX := GoldSboxPower3();
"""

MAGMA_MAIN_LOOP = r"""
print "=======================================================";
print "INVARIANT TEST: " cat test_label;
print "Reference: Gold x^3, modulus 285";
print "Functions tested: " cat IntegerToString(#FOUND_SBOXES);
print "Output file: " cat outfile;
print "=======================================================";

PrintFile(outfile, "=======================================================");
PrintFile(outfile, "INVARIANT TEST: " cat test_label);
PrintFile(outfile, "Modulus: 285 = 0x11D");
PrintFile(outfile, "Functions tested: " cat IntegerToString(#FOUND_SBOXES));
PrintFile(outfile, "=======================================================");

t0 := Cputime();

// Gold reference
gold_rec := FullInvariantRecord(GOLD_SBOX);
print "Gold done:", Cputime(t0), "sec";
PrintFile(outfile, "REFERENCE: Gold x^3");
PrintRecord("Gold x^3", gold_rec, gold_rec, outfile);

// All found functions
records := [];
for item in FOUND_SBOXES do
    label := item[1];
    S     := item[2];
    t := Cputime();
    rec := FullInvariantRecord(S);
    Append(~records, <label, rec>);
    PrintRecord(label, rec, gold_rec, outfile);
    msg := label cat " done in " cat Sprint(Cputime(t)) cat " sec";
    print msg;
    PrintFile(outfile, "  time: " cat Sprint(Cputime(t)) cat " sec");
end for;

// Pairwise image_size comparison (cheap discriminator)
PrintFile(outfile, "=======================================================");
PrintFile(outfile, "PAIRWISE INVARIANT EQUALITY");
PrintFile(outfile, "=======================================================");
for i in [1..#FOUND_SBOXES] do
    for j in [i+1..#FOUND_SBOXES] do
        ri := records[i][2]; rj := records[j][2];
        same := (ri["perm"] eq rj["perm"]) and (ri["image"] eq rj["image"]) and
                (ri["DU"] eq rj["DU"]) and (ri["DDT"] eq rj["DDT"]) and
                (ri["WalshAbs"] eq rj["WalshAbs"]) and
                (ri["alg_deg"] eq rj["alg_deg"]);
        if COMPUTE_AUTOCORR then same := same and (ri["ACAbs"] eq rj["ACAbs"]); end if;
        PrintFile(outfile, FOUND_SBOXES[i][1] cat " vs " cat FOUND_SBOXES[j][1] cat
                  " : all_invariants_equal=" cat (same select "true" else "false") cat
                  "  image=" cat IntegerToString(ri["image"]) cat
                  " vs " cat IntegerToString(rj["image"]));
    end for;
end for;

PrintFile(outfile, "=======================================================");
PrintFile(outfile, "DONE. Total time: " cat Sprint(Cputime(t0)) cat " sec");
PrintFile(outfile, "=======================================================");
print "DONE. Total time:", Cputime(t0), "sec";
print "Results written to:", outfile;
quit;
"""


def format_sbox_magma(fid, sbox, ccz_class=None):
    """Format one sbox entry for Magma FOUND_SBOXES list."""
    label = fid
    if ccz_class:
        label = f"{fid}_{ccz_class}"
    nums = ', '.join(str(x) for x in sbox)
    # Wrap at 16 per line for readability
    rows = []
    for i in range(0, 256, 16):
        rows.append('        ' + ', '.join(str(x) for x in sbox[i:i+16]))
    body = ',\n'.join(rows)
    return f'    <"{label}", [\n{body}\n    ]>'


def generate_invariant_test(funcs, label, outfile_rel):
    """Generate a complete Magma invariant test script."""
    entries = []
    for fid, v in funcs.items():
        # label = "APN-0006_CLASS-004" etc
        entries.append(format_sbox_magma(fid, v['sbox']))

    sboxes_block = 'FOUND_SBOXES := [\n' + ',\n'.join(entries) + '\n];\n'

    script = (
        f'/* Auto-generated by gen_test_magma.py\n'
        f' * Test: {label}\n'
        f' * Functions: {len(funcs)}\n'
        f' * Generated: {datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")}\n'
        f' * Results -> {outfile_rel}\n'
        f' */\n\n'
        f'outfile := "{outfile_rel}";\n'
        f'test_label := "{label}";\n\n'
    )
    script += MAGMA_LIBRARY + '\n'
    script += sboxes_block + '\n'
    script += MAGMA_MAIN_LOOP

    return script


def generate_ccz_pairwise(funcs, label, outfile_rel):
    """Generate pairwise CCZ test using dual graph-code isomorphism."""
    entries = []
    for fid, v in funcs.items():
        entries.append(format_sbox_magma(fid, v['sbox']))

    sboxes_block = 'FOUND_SBOXES := [\n' + ',\n'.join(entries) + '\n];\n'

    script = (
        f'/* Auto-generated by gen_test_magma.py\n'
        f' * Test: pairwise CCZ - {label}\n'
        f' * Functions: {len(funcs)}\n'
        f' * Generated: {datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")}\n'
        f' * Results -> {outfile_rel}\n'
        f' */\n\n'
        f'SetVerbose("Code", 0);\n'
        f'outfile := "{outfile_rel}";\n'
        f'test_label := "{label}";\n'
        f'n := 8;\nN := 2^n;\n\n'
    )

    script += r"""
function IntToBitsLE(a, n)
    bits := [];
    x := Integers()!a;
    for i in [1..n] do Append(~bits, x mod 2); x := x div 2; end for;
    return bits;
end function;

function GraphParityCheckMatrix(S)
    F2 := GF(2);
    H := ZeroMatrix(F2, 1 + 2*n, N);
    for x in [0..N-1] do
        col := x + 1;
        H[1, col] := F2!1;
        xb := IntToBitsLE(x, n);
        yb := IntToBitsLE(S[col], n);
        for i in [1..n] do
            H[1+i, col]   := F2!xb[i];
            H[1+n+i, col] := F2!yb[i];
        end for;
    end for;
    return H;
end function;

function GraphDualCode(S)
    return LinearCode(GraphParityCheckMatrix(S));
end function;

function WDCompact(C)
    WD := WeightDistribution(C);
    return [ <Integers()!t[1], Integers()!t[2]> : t in WD | Integers()!t[2] ne 0 ];
end function;

function FastCCZEquivalent(S, T)
    C1 := GraphDualCode(S); C2 := GraphDualCode(T);
    if WDCompact(C1) ne WDCompact(C2) then return false; end if;
    ok := IsIsomorphic(C1, C2);
    return ok;
end function;

"""
    script += sboxes_block + '\n'
    script += f"""
PrintFile(outfile, "=======================================================");
PrintFile(outfile, "PAIRWISE CCZ TEST: {label}");
PrintFile(outfile, "Method: dual graph-code isomorphism");
PrintFile(outfile, "Functions: " cat IntegerToString(#FOUND_SBOXES));
PrintFile(outfile, "=======================================================");

t0 := Cputime();
n_equiv := 0; n_tests := 0;
for i in [1..#FOUND_SBOXES] do
    for j in [i+1..#FOUND_SBOXES] do
        t := Cputime();
        ok := FastCCZEquivalent(FOUND_SBOXES[i][2], FOUND_SBOXES[j][2]);
        n_tests +:= 1;
        if ok then n_equiv +:= 1; end if;
        line := FOUND_SBOXES[i][1] cat " vs " cat FOUND_SBOXES[j][1] cat
                " : ccz_equivalent=" cat (ok select "true" else "false") cat
                "  time=" cat Sprint(Cputime(t)) cat "s";
        PrintFile(outfile, line);
        print line;
    end for;
end for;

PrintFile(outfile, "=======================================================");
PrintFile(outfile, "SUMMARY: " cat IntegerToString(n_tests) cat " tests, " cat
          IntegerToString(n_equiv) cat " equivalent pairs");
PrintFile(outfile, "Total time: " cat Sprint(Cputime(t0)) cat " sec");
PrintFile(outfile, "=======================================================");
print "DONE.", n_tests, "tests,", n_equiv, "equivalent pairs";
print "Total time:", Cputime(t0), "sec";
quit;
"""
    return script


# ---------------------------------------------------------------------------
# Parse results and save JSON
# ---------------------------------------------------------------------------

def parse_raw_results(raw_txt_path):
    """Parse invariant test raw output into structured dict."""
    if not os.path.exists(raw_txt_path):
        return None

    results = {'functions': {}, 'pairwise': {}}
    current = None

    with open(raw_txt_path, 'r', encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('---'):
                current = None
            # New function block
            elif line.strip() and not line.startswith('=') and not line.startswith('PAIRWISE') \
                    and not line.startswith('DONE') and not line.startswith('INVARIANT') \
                    and not line.startswith('Reference') and not line.startswith('Functions') \
                    and not line.startswith('Modulus') and not line.startswith('REFERENCE') \
                    and ':' not in line and 'vs' not in line:
                name = line.strip()
                if name:
                    current = name
                    results['functions'][name] = {}
            elif current and ':' in line:
                key, _, val = line.partition(':')
                key = key.strip().lstrip()
                val = val.strip()
                results['functions'][current][key] = val
            elif 'vs' in line and 'all_invariants_equal' in line:
                parts = line.split(':')
                if len(parts) >= 2:
                    pair = parts[0].strip()
                    equal = 'true' in parts[1]
                    results['pairwise'][pair] = equal
    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Generate and run Magma invariant/CCZ tests from apn_found.json')
    parser.add_argument('--gen', action='store_true',
                        help='Generate .m test files only')
    parser.add_argument('--run', action='store_true',
                        help='Generate and run tests (requires Magma)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Show what would be generated')
    parser.add_argument('--class', dest='class_filter',
                        help='Only generate for this CCZ class (e.g. CLASS-004)')
    parser.add_argument('--magma', default='magma',
                        help='Path to Magma binary (default: magma)')
    parser.add_argument('--data', default=str(FOUND_JSON),
                        help=f'Path to apn_found.json (default: {FOUND_JSON})')
    parser.add_argument('--ccz-pairwise', action='store_true',
                        help='Also generate pairwise CCZ test for filtered class')
    args = parser.parse_args()

    if not args.gen and not args.run and not args.dry_run:
        parser.print_help()
        sys.exit(1)

    # Load data
    funcs = load_functions(Path(args.data))
    class_map = load_class_map()
    print(f'Loaded {len(funcs)} functions from {args.data}')

    # Create output directories
    GEN_DIR.mkdir(exist_ok=True)
    VERIF_DIR.mkdir(parents=True, exist_ok=True)

    # Determine what to generate
    tasks = []

    if args.class_filter:
        subset = filter_by_class(funcs, class_map, args.class_filter)
        if not subset:
            print(f'ERROR: no functions found for class {args.class_filter}')
            sys.exit(1)
        cls_lower = args.class_filter.lower().replace('-', '')
        m_name = f'invariant_test_{cls_lower}.m'
        out_rel = f'data/verification/invariant_{cls_lower}.txt'
        label = f'Invariant test {args.class_filter} ({len(subset)} functions)'
        tasks.append(('invariant', subset, m_name, out_rel, label))

        if args.ccz_pairwise:
            m_name_ccz = f'ccz_pairwise_{cls_lower}.m'
            out_rel_ccz = f'data/verification/ccz_pairwise_{cls_lower}.txt'
            label_ccz = f'Pairwise CCZ {args.class_filter} ({len(subset)} functions)'
            tasks.append(('ccz', subset, m_name_ccz, out_rel_ccz, label_ccz))
    else:
        # Generate for all 25
        m_name = 'invariant_test_all25.m'
        out_rel = 'data/verification/invariant_all25.txt'
        label = f'Invariant test ALL ({len(funcs)} functions)'
        tasks.append(('invariant', funcs, m_name, out_rel, label))

        # Also generate per-class
        for cls_id in sorted(set(class_map.values())):
            subset = filter_by_class(funcs, class_map, cls_id)
            if not subset:
                continue
            cls_lower = cls_id.lower().replace('-', '')
            m_name_c = f'invariant_test_{cls_lower}.m'
            out_rel_c = f'data/verification/invariant_{cls_lower}.txt'
            label_c = f'Invariant test {cls_id} ({len(subset)} functions)'
            tasks.append(('invariant', subset, m_name_c, out_rel_c, label_c))

    # Execute
    for task_type, subset, m_name, out_rel, label in tasks:
        m_path = GEN_DIR / m_name
        out_abs = REPO_ROOT / out_rel

        if args.dry_run:
            print(f'  [DRY] Would generate: {m_path}')
            print(f'  [DRY] Results -> {out_abs}')
            print(f'  [DRY] Functions: {list(subset.keys())}')
            continue

        # Generate
        if task_type == 'invariant':
            script = generate_invariant_test(subset, label, out_rel)
        else:
            script = generate_ccz_pairwise(subset, label, out_rel)

        # Write ASCII-only (Magma requirement)
        clean = ''.join(ch if ord(ch) < 128 else '_' for ch in script)
        m_path.write_text(clean, encoding='ascii')
        print(f'Generated: {m_path}  ({len(subset)} functions)')
        print(f'Results -> {out_abs}')

        if args.run:
            print(f'Running: {args.magma} {m_path} ...')
            try:
                result = subprocess.run(
                    [args.magma, str(m_path)],
                    capture_output=True, text=True, timeout=7200
                )
                if result.returncode != 0 and result.stderr:
                    print(f'  Magma stderr: {result.stderr[:500]}')
                print(f'  Done. Exit code: {result.returncode}')

                # Save raw stdout alongside the results file
                raw_path = str(out_abs).replace('.txt', '_stdout.txt')
                Path(raw_path).write_text(result.stdout, encoding='utf-8')
                print(f'  Stdout saved to: {raw_path}')

            except subprocess.TimeoutExpired:
                print(f'  TIMEOUT after 2 hours')
            except FileNotFoundError:
                print(f'  ERROR: Magma not found at "{args.magma}"')
                print(f'  Use --magma /path/to/magma to specify the binary')

    if not args.dry_run:
        print()
        print('To run in Magma manually:')
        for _, _, m_name, out_rel, _ in tasks:
            print(f'  magma {GEN_DIR / m_name}')
        print()
        print('Results will be saved to:')
        for _, _, _, out_rel, _ in tasks:
            print(f'  {REPO_ROOT / out_rel}')


if __name__ == '__main__':
    main()
