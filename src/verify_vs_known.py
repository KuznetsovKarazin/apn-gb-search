#!/usr/bin/env sage -python
# verify_6classes_vs_known_controls_safe.py
#
# Safe verification:
#  - six class representatives pairwise;
#  - six representatives vs Gold/Kasami-style controls and optional first DB rows.
#
# It computes the ortho-derivative signature first and runs exact
# sboxU.are_ccz_equivalent_from_code only if signatures match.
#
# Usage:
#   cd /mnt/e/Magma
#   conda activate sage
#   sage -python verify_6classes_vs_known_controls_safe.py
#
# Optional:
#   sage -python verify_6classes_vs_known_controls_safe.py --also-db-first 20
#
# Dangerous full exact mode, can take hours:
#   sage -python verify_6classes_vs_known_controls_safe.py --force-exact

import argparse
import hashlib
import json
import os
import re
import time

import sboxU
from sage.all import GF

OUT_TXT = "verify_6classes_vs_known_controls_report.txt"
OUT_JSON = "verify_6classes_vs_known_controls_report.json"

CLASS_REPS = [
    ("CLASS-001", "APN-0001"),
    ("CLASS-002", "APN-0002"),
    ("CLASS-003", "APN-0004"),
    ("CLASS-004", "APN-0006"),
    ("CLASS-005", "APN-0007"),
    ("CLASS-006", "APN-0010"),
]


def h16(sb):
    return hashlib.sha256(bytes(sb)).hexdigest()[:16]


def stable(x):
    try:
        return json.dumps(x, sort_keys=True, default=str)
    except Exception:
        return repr(x)


def ortho_signature(sb):
    pi = sboxU.ortho_derivative(sb)
    return (
        stable(sboxU.differential_spectrum(pi)),
        stable(sboxU.absolute_walsh_spectrum(pi)),
    )


def sig_hash(sig):
    return hashlib.sha256(repr(sig).encode()).hexdigest()[:16]


def load_found(path):
    data = json.load(open(path, "r", encoding="utf-8", errors="replace"))
    funcs = data.get("functions", {})
    records = list(funcs.values()) if isinstance(funcs, dict) else funcs if isinstance(funcs, list) else []
    out = {}
    for i, rec in enumerate(records, 1):
        if not isinstance(rec, dict):
            continue
        sb = rec.get("sbox")
        if not isinstance(sb, list) or len(sb) != 256:
            continue
        aid = rec.get("id")
        if isinstance(aid, int):
            aid = f"APN-{aid:04d}"
        elif not aid:
            aid = f"APN-{i:04d}"
        sb = [int(x) for x in sb]
        out[aid] = {"id": aid, "sbox": sb, "hash": h16(sb)}
    return out


def gf256_power_sbox(exp):
    F = GF(2**8, name="a")
    a = F.gen()

    def int_to_field(x):
        y = F(0)
        p = F(1)
        for i in range(8):
            if (x >> i) & 1:
                y += p
            p *= a
        return y

    def field_to_int(y):
        v = list(y.vector()) if hasattr(y, "vector") else list(y._vector_())
        z = 0
        for i, b in enumerate(v):
            if int(b) & 1:
                z |= (1 << i)
        return z

    return [field_to_int(int_to_field(x) ** exp) for x in range(256)]


def parse_sbox_line(line):
    nums = [int(x) for x in re.findall(r"\d+", line)]
    if len(nums) == 256 and all(0 <= x <= 255 for x in nums):
        return nums
    return None


def load_db_first(path, limit):
    out = []
    if limit <= 0 or not os.path.exists(path):
        return out
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line_no, line in enumerate(f, 1):
            sb = parse_sbox_line(line)
            if sb is None:
                continue
            out.append((f"DB-line-{line_no}", sb))
            if len(out) >= limit:
                break
    return out


def info(name, sb):
    sig = ortho_signature(sb)
    return {
        "name": name,
        "hash": h16(sb),
        "du": int(sboxU.differential_uniformity(sb)),
        "degree": int(sboxU.algebraic_degree(sb)),
        "perm": len(set(sb)) == 256,
        "first8": sb[:8],
        "sig_hash": sig_hash(sig),
        "_sig": sig,
    }


def compare(na, sa, ia, nb, sb, ib, force_exact=False):
    same = ia["_sig"] == ib["_sig"]
    if not same and not force_exact:
        return {
            "a": na, "b": nb, "same_ortho_signature": False,
            "ccz": False, "exact_run": False,
            "reason": "different_ortho_signature", "time_sec": 0.0, "error": None,
        }

    t = time.time()
    try:
        eq = bool(sboxU.are_ccz_equivalent_from_code(sa, sb))
        err = None
    except Exception as e:
        eq = None
        err = repr(e)
    return {
        "a": na, "b": nb, "same_ortho_signature": same,
        "ccz": eq, "exact_run": True,
        "reason": "same_ortho_signature_exact_ccz" if same else "forced_exact",
        "time_sec": time.time() - t, "error": err,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--found", default="apn_found.json")
    ap.add_argument("--known-db", default="new_apns.txt")
    ap.add_argument("--also-db-first", type=int, default=0)
    ap.add_argument("--force-exact", action="store_true")
    args = ap.parse_args()

    t0 = time.time()
    lines = []

    found = load_found(args.found)
    missing = [aid for _, aid in CLASS_REPS if aid not in found]
    if missing:
        raise RuntimeError(f"Missing representatives: {missing}")

    reps = [(f"{cid}:{aid}", found[aid]["sbox"]) for cid, aid in CLASS_REPS]

    controls = []
    for k in [1, 3, 5, 7]:
        exp = 2**k + 1
        controls.append((f"Gold_k{k}_x^{exp}", gf256_power_sbox(exp)))
    for k in [1, 3, 5, 7]:
        exp = (2**(2*k) - 2**k + 1) % 255 or 255
        controls.append((f"Kasami_k{k}_x^{exp}", gf256_power_sbox(exp)))
    for name, exp in [
        ("Cube_x^3", 3),
        ("Power_x^9", 9),
        ("Power_x^57", 57),
        ("Power_x^113", 113),
        ("Inverse_x^254_control", 254),
    ]:
        controls.append((name, gf256_power_sbox(exp)))

    uniq = {}
    for name, sb in controls:
        uniq.setdefault(h16(sb), (name, sb))
    controls = list(uniq.values())
    controls += load_db_first(args.known_db, args.also_db_first)

    lines.append("SAFE VERIFICATION: 6 APN CLASSES VS KNOWN CONTROLS")
    lines.append("=" * 90)
    lines.append(f"force_exact: {args.force_exact}")
    lines.append(f"also_db_first: {args.also_db_first}")
    lines.append("")

    print("Computing representative signatures...")
    rep_info = {}
    lines.append("REPRESENTATIVES")
    for name, sb in reps:
        ii = info(name, sb)
        rep_info[name] = ii
        msg = f"{name}: hash={ii['hash']} DU={ii['du']} degree={ii['degree']} perm={ii['perm']} sig={ii['sig_hash']} first8={ii['first8']}"
        print(msg)
        lines.append(msg)

    print("\nComputing known-control signatures...")
    ctrl_info = {}
    lines.append("\nKNOWN CONTROLS")
    for name, sb in controls:
        ii = info(name, sb)
        ctrl_info[name] = ii
        msg = f"{name}: hash={ii['hash']} DU={ii['du']} degree={ii['degree']} perm={ii['perm']} sig={ii['sig_hash']} first8={ii['first8']}"
        print(msg)
        lines.append(msg)

    pairwise = []
    lines.append("\nPAIRWISE BETWEEN SIX REPRESENTATIVES")
    print("\nPairwise checks...")
    for i in range(len(reps)):
        for j in range(i + 1, len(reps)):
            na, sa = reps[i]
            nb, sb = reps[j]
            r = compare(na, sa, rep_info[na], nb, sb, rep_info[nb], args.force_exact)
            pairwise.append(r)
            msg = f"{na} vs {nb}: ccz={r['ccz']} same_sig={r['same_ortho_signature']} exact_run={r['exact_run']} reason={r['reason']} time={r['time_sec']:.2f}s"
            print(msg, flush=True)
            lines.append(msg)

    controls_res = []
    lines.append("\nSIX REPRESENTATIVES VS KNOWN CONTROLS")
    print("\nKnown-control checks...")
    for na, sa in reps:
        for nb, sb in controls:
            r = compare(na, sa, rep_info[na], nb, sb, ctrl_info[nb], args.force_exact)
            controls_res.append(r)
            msg = f"{na} vs {nb}: ccz={r['ccz']} same_sig={r['same_ortho_signature']} exact_run={r['exact_run']} reason={r['reason']} time={r['time_sec']:.2f}s"
            print(msg, flush=True)
            lines.append(msg)

    any_pair = any(r["ccz"] is True for r in pairwise)
    all_non = all(r["ccz"] is False for r in pairwise)
    any_control = any(r["ccz"] is True for r in controls_res)
    exact_runs = sum(1 for r in pairwise + controls_res if r["exact_run"])
    elapsed = time.time() - t0

    lines.append("\nSUMMARY")
    lines.append("=" * 90)
    lines.append(f"Any equivalent pair among six: {any_pair}")
    lines.append(f"All six pairwise non-equivalent: {all_non}")
    lines.append(f"Any known-control match: {any_control}")
    lines.append(f"Exact CCZ runs performed: {exact_runs}")
    lines.append(f"Elapsed seconds: {elapsed:.2f}")

    report = {
        "method": "ortho_signature prefilter; exact CCZ only on same signature unless --force-exact",
        "force_exact": args.force_exact,
        "also_db_first": args.also_db_first,
        "representatives": [{k: v for k, v in rep_info[n].items() if k != "_sig"} for n, _ in reps],
        "known_controls": [{k: v for k, v in ctrl_info[n].items() if k != "_sig"} for n, _ in controls],
        "pairwise_between_six": pairwise,
        "against_known_controls": controls_res,
        "any_equivalent_pair_among_six": any_pair,
        "all_six_pairwise_non_equivalent": all_non,
        "any_known_control_match": any_control,
        "exact_runs": exact_runs,
        "elapsed_sec": elapsed,
    }

    with open(OUT_TXT, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    with open(OUT_JSON, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)

    print("\nDONE")
    print(f"Any equivalent pair among six: {any_pair}")
    print(f"All six pairwise non-equivalent: {all_non}")
    print(f"Any known-control match: {any_control}")
    print(f"Exact CCZ runs performed: {exact_runs}")
    print(f"Report TXT : {OUT_TXT}")
    print(f"Report JSON: {OUT_JSON}")


if __name__ == "__main__":
    main()
