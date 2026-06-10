#!/usr/bin/env python3
# check_apn_known_orthoderivative_stream.py
#
# Streaming comparison of found quadratic APN functions against a huge known APN database.
#
# Method:
#   For quadratic APN F, CCZ-equivalence coincides with EA-equivalence.
#   EA-equivalent quadratic APN functions have affine-equivalent ortho-derivatives.
#   The practical invariant used in Canteaut-Couvreur-Perrin/sboxU is:
#       differential spectrum of ortho-derivative
#       absolute Walsh spectrum of ortho-derivative
#
# This script does NOT run full CCZ on every pair.
# It streams known functions line by line and compares only the ortho-derivative invariant.
#
# Usage:
#   cd /mnt/e/Magma
#   conda activate sage
#
# Small test:
#   python check_apn_known_orthoderivative_stream.py --found apn_found.json --known new_apns.txt --max-known 1000 --progress 100
#
# Full run:
#   python check_apn_known_orthoderivative_stream.py --found apn_found.json --known new_apns.txt --progress 10000 --skip-pairwise
#
# Optional: if invariant hits are found, verify those few hits by exact code-equivalence:
#   python check_apn_known_orthoderivative_stream.py --found apn_found.json --known new_apns.txt --progress 10000 --exact-on-hit

import argparse
import hashlib
import json
import re
import sys
import time
from collections import defaultdict

import sboxU

N = 256


def sbox_hash(sb):
    return hashlib.sha256(bytes(sb)).hexdigest()[:16]


def normalize_sbox(xs):
    try:
        sb = [int(x) for x in xs]
    except Exception:
        return None
    if len(sb) == N and all(0 <= x < N for x in sb):
        return sb
    return None


def parse_sbox_line(line):
    # Expected: [0,0,...,234]
    nums = re.findall(r"\d+", line)
    return normalize_sbox(nums)


def load_found(path):
    data = json.load(open(path, "r", encoding="utf-8", errors="replace"))
    funcs = data.get("functions", {})
    out = []
    if isinstance(funcs, dict):
        for aid, rec in sorted(funcs.items()):
            sb = normalize_sbox(rec.get("sbox", []))
            if sb is not None:
                out.append((aid, sb))
    elif isinstance(funcs, list):
        for i, rec in enumerate(funcs, 1):
            sb = normalize_sbox(rec.get("sbox", []))
            if sb is not None:
                out.append((str(rec.get("id", f"APN-{i:04d}")), sb))
    return out


def stable(obj):
    # sboxU/Sage may return tuples, dict-like spectra, Sage integers.
    # Convert to stable comparable string.
    try:
        return json.dumps(obj, sort_keys=True, default=str)
    except Exception:
        return repr(obj)


def spectrum_as_stable(name, obj):
    return (name, stable(obj))


def ortho_signature(sb):
    """
    Main invariant:
      (differential spectrum of ortho-derivative,
       absolute Walsh spectrum of ortho-derivative)

    Returns stable tuple suitable for dict/set comparison.
    """
    pi = sboxU.ortho_derivative(sb)
    ds = sboxU.differential_spectrum(pi)
    ws = sboxU.absolute_walsh_spectrum(pi)
    return (
        spectrum_as_stable("diff_ortho", ds),
        spectrum_as_stable("walsh_ortho", ws),
    )


def short_sig(sig):
    # For printing compactly.
    h = hashlib.sha256(repr(sig).encode()).hexdigest()[:16]
    return h


def exact_ccz(a, b):
    return bool(sboxU.are_ccz_equivalent_from_code(a, b))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--found", default="apn_found.json")
    ap.add_argument("--known", required=True)
    ap.add_argument("--max-known", type=int, default=0)
    ap.add_argument("--progress", type=int, default=10000)
    ap.add_argument("--exact-on-hit", action="store_true")
    ap.add_argument("--skip-pairwise", action="store_true", help="Skip pairwise exact CCZ tests among found functions at the end")
    ap.add_argument("--report-prefix", default="orthoderivative_known")
    args = ap.parse_args()

    found = load_found(args.found)
    if not found:
        print("ERROR: no found functions loaded from", args.found, flush=True)
        sys.exit(1)

    print("Loaded found candidates:", len(found), flush=True)

    found_sig = {}
    found_by_sig = defaultdict(list)
    for name, sb in found:
        t = time.time()
        try:
            sig = ortho_signature(sb)
        except Exception as e:
            print(f"ERROR computing signature for {name}: {e}", flush=True)
            raise
        found_sig[name] = sig
        found_by_sig[sig].append((name, sb))
        print(
            f"  {name}: sbox_hash={sbox_hash(sb)} ortho_sig_hash={short_sig(sig)} time={time.time()-t:.3f}s",
            flush=True,
        )

    print("", flush=True)
    print("Streaming known database:", args.known, flush=True)
    print("Exact CCZ on invariant hit:", args.exact_on_hit, flush=True)
    print("Progress every:", args.progress, "known S-boxes", flush=True)
    print("", flush=True)

    matches = {name: [] for name, _ in found}
    exact_tests = {name: 0 for name, _ in found}

    parsed = 0
    bad = 0
    sig_errors = 0
    t0 = time.time()
    last_t = t0

    with open(args.known, "r", encoding="utf-8", errors="replace") as f:
        for line_no, line in enumerate(f, 1):
            if args.max_known and parsed >= args.max_known:
                break

            line = line.strip()
            if not line:
                continue

            sb = parse_sbox_line(line)
            if sb is None:
                bad += 1
                continue

            parsed += 1

            try:
                sig = ortho_signature(sb)
            except Exception as e:
                sig_errors += 1
                if sig_errors <= 10:
                    print(f"WARNING signature failed at line {line_no}: {e}", flush=True)
                continue

            if sig in found_by_sig:
                kh = sbox_hash(sb)
                for cname, csb in found_by_sig[sig]:
                    exact = None
                    if args.exact_on_hit:
                        exact_tests[cname] += 1
                        exact = exact_ccz(csb, sb)
                        if not exact:
                            # invariant collision but not equivalent by exact code criterion
                            pass

                    hit = {
                        "known_line": line_no,
                        "known_index": parsed,
                        "known_hash": kh,
                        "ortho_sig_hash": short_sig(sig),
                        "exact_ccz": exact,
                    }
                    matches[cname].append(hit)
                    print(
                        f"HIT: {cname} ~ known_line={line_no} known_index={parsed} known_hash={kh} exact={exact}",
                        flush=True,
                    )

            if parsed % args.progress == 0:
                now = time.time()
                total_dt = now - t0
                step_dt = now - last_t
                last_t = now
                rate = args.progress / step_dt if step_dt > 0 else 0.0
                print(
                    f"parsed={parsed} bad={bad} sig_errors={sig_errors} elapsed={total_dt:.1f}s rate={rate:.1f}/s",
                    flush=True,
                )
                for cname, _ in found:
                    print(
                        f"  {cname}: hits={len(matches[cname])} exact_tests={exact_tests[cname]}",
                        flush=True,
                    )
                print("", flush=True)

    # Pairwise among found by invariant first; exact only if signatures match.
    # For huge known-database runs this can be skipped with --skip-pairwise,
    # because exact code equivalence may take much longer than the streaming invariant scan.
    if args.skip_pairwise:
        print("Pairwise among found: skipped (--skip-pairwise)", flush=True)
        pairwise = []
        groups = [[name] for name, _ in found]
    else:
        print("Pairwise among found:", flush=True)
        pairwise = []
        parent = {name: name for name, _ in found}

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def union(a, b):
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[rb] = ra

        for i in range(len(found)):
            for j in range(i + 1, len(found)):
                n1, s1 = found[i]
                n2, s2 = found[j]
                same_sig = found_sig[n1] == found_sig[n2]
                if same_sig:
                    try:
                        eq = exact_ccz(s1, s2)
                        reason = "same_ortho_signature_exact_ccz"
                    except Exception as e:
                        eq = None
                        reason = f"same_ortho_signature_exact_failed:{e}"
                else:
                    eq = False
                    reason = "different_ortho_signature"
                print(f"  {n1} vs {n2}: {eq} ({reason})", flush=True)
                pairwise.append({"a": n1, "b": n2, "same_ortho_signature": same_sig, "ccz": eq, "reason": reason})
                if eq is True:
                    union(n1, n2)

        groups_dict = defaultdict(list)
        for name, _ in found:
            groups_dict[find(name)].append(name)
        groups = list(groups_dict.values())

    elapsed = time.time() - t0

    result = {
        "method": "ortho_derivative differential_spectrum + absolute_walsh_spectrum",
        "known_file": args.known,
        "parsed_known": parsed,
        "bad_lines": bad,
        "signature_errors": sig_errors,
        "elapsed_sec": elapsed,
        "exact_on_hit": args.exact_on_hit,
        "matches": matches,
        "exact_tests": exact_tests,
        "pairwise": pairwise,
        "groups": groups,
    }

    out_json = args.report_prefix + "_report.json"
    out_txt = args.report_prefix + "_report.txt"

    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)

    with open(out_txt, "w", encoding="utf-8") as f:
        f.write("ORTHO-DERIVATIVE INVARIANT CHECK AGAINST KNOWN APN DATABASE\n")
        f.write("=" * 80 + "\n")
        f.write("method: differential_spectrum(ortho_derivative(F)) + absolute_walsh_spectrum(ortho_derivative(F))\n")
        f.write(f"known_file: {args.known}\n")
        f.write(f"parsed_known: {parsed}\n")
        f.write(f"bad_lines: {bad}\n")
        f.write(f"signature_errors: {sig_errors}\n")
        f.write(f"elapsed_sec: {elapsed:.2f}\n")
        f.write(f"exact_on_hit: {args.exact_on_hit}\n\n")

        for cname, _ in found:
            status = "KNOWN_SIGNATURE_MATCH" if matches[cname] else "NO_SIGNATURE_MATCH_IN_DATABASE"
            f.write(f"{cname}: {status}\n")
            f.write(f"  hits: {len(matches[cname])}\n")
            f.write(f"  exact_tests: {exact_tests[cname]}\n")
            f.write(f"  matches: {matches[cname]}\n\n")

        f.write("Pairwise among found:\n")
        for p in pairwise:
            f.write(f"  {p['a']} vs {p['b']}: same_sig={p['same_ortho_signature']} ccz={p['ccz']} reason={p['reason']}\n")

        f.write("\nGroups:\n")
        for g in groups:
            f.write("  " + ", ".join(g) + "\n")

    print("", flush=True)
    print("DONE", flush=True)
    print("parsed_known:", parsed, flush=True)
    print("bad_lines:", bad, flush=True)
    print("signature_errors:", sig_errors, flush=True)
    print("elapsed_sec:", f"{elapsed:.2f}", flush=True)
    print("Report:", out_txt, flush=True)
    print("Report:", out_json, flush=True)
    print("", flush=True)
    for cname, _ in found:
        print(cname, "=>", "MATCH" if matches[cname] else "NO_MATCH", "hits:", len(matches[cname]), flush=True)


if __name__ == "__main__":
    main()
