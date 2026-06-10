#!/usr/bin/env python3
"""
classify_apn.py
===============
Classify found APN functions into CCZ-classes using sboxU ortho-derivative.

For quadratic APN: CCZ = EA-equivalence (Yoshiara 2010).
Ortho-derivative signature is an exact CCZ-invariant.
Exact CCZ test (are_ccz_equivalent_from_code) only runs within same-sig buckets.

Usage:
  conda activate sage
  python src/classify_apn.py                                         # auto-detect
  python src/classify_apn.py --found batch/batch_004/results/apn_found.json
  python src/classify_apn.py --found results/merged/apn_all.json

Output (same folder as input file):
  found_apn_classes_report.txt
  found_apn_classes_report.json
"""

import argparse
import hashlib
import json
import time
from collections import defaultdict
from pathlib import Path

import sboxU

N = 256


#    Data loading                                                               

def load_functions(path):
    """
    Load APN functions from JSON. Supports all formats:
      - list:  [{"sbox": [...], "tag": ..., "total": ...}, ...]
      - dict:  {"APN-0001": {"sbox": [...], ...}, ...}
    Returns list of dicts each with keys: id, tag, sol, total, perm, hash, sbox
    """
    data = json.load(open(path, encoding="utf-8", errors="replace"))
    fns  = data.get("functions", {})

    if isinstance(fns, dict):
        records = list(fns.values())
    elif isinstance(fns, list):
        records = fns
    else:
        return []

    out = []
    for i, rec in enumerate(records, 1):
        if not isinstance(rec, dict):
            continue
        sb = rec.get("sbox", [])
        if not isinstance(sb, list) or len(sb) != N:
            continue
        sb = [int(x) for x in sb]

        # Resolve ID
        aid = rec.get("apn_id") or rec.get("id")
        if isinstance(aid, int):
            aid = f"APN-{aid:04d}"
        elif not aid:
            aid = f"APN-{i:04d}"

        out.append({
            "id":    aid,
            "tag":   rec.get("tag") or rec.get("seed") or "",
            "sol":   rec.get("sol", -1),
            "total": rec.get("total", rec.get("total_in_slice", -1)),
            "perm":  rec.get("perm", False),
            "batch": rec.get("batch", rec.get("_batch", "")),
            "hash":  rec.get("hash") or rec.get("sbox_hash") or
                     hashlib.sha256(bytes(sb)).hexdigest()[:16],
            "sbox":  sb,
        })

    out.sort(key=lambda r: r["id"])
    return out


#    sboxU helpers                                                              

def stable(obj):
    try:
        return json.dumps(obj, sort_keys=True, default=str)
    except Exception:
        return repr(obj)

def ortho_signature(sb):
    pi = sboxU.ortho_derivative(sb)
    return (
        stable(sboxU.differential_spectrum(pi)),
        stable(sboxU.absolute_walsh_spectrum(pi)),
    )

def sig_hash(sig):
    return hashlib.sha256(repr(sig).encode()).hexdigest()[:16]

def exact_ccz(a, b):
    return bool(sboxU.are_ccz_equivalent_from_code(a, b))

def union_find_groups(items, edges):
    parent = {x: x for x in items}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb: parent[rb] = ra
    for a, b, eq in edges:
        if eq: union(a, b)
    groups = defaultdict(list)
    for x in items:
        groups[find(x)].append(x)
    return list(groups.values())


#    Main                                                                       

def find_latest_apn_found(repo_root):
    for b in reversed(sorted(Path(repo_root, "batch").glob("batch_*"))):
        f = b / "results" / "apn_found.json"
        if f.exists(): return f
    return None

def main():
    repo = Path(__file__).resolve().parent.parent

    parser = argparse.ArgumentParser(description="Classify APN functions via sboxU.")
    parser.add_argument("--found", type=Path, default=None,
                        help="Path to apn_found.json or apn_all.json")
    parser.add_argument("--skip-intra-bucket", action="store_true",
                        help="Skip exact CCZ inside buckets (fast; safe for quadratic APN)")
    args = parser.parse_args()

    if args.found:
        found_path = args.found.resolve()
    else:
        found_path = find_latest_apn_found(repo)
        if not found_path:
            print("No apn_found.json found. Run collect_results.py first.")
            return
        print(f"Auto-detected: {found_path}")

    out_dir  = found_path.parent
    out_txt  = str(out_dir / "found_apn_classes_report.txt")
    out_json = str(out_dir / "found_apn_classes_report.json")

    #    Load                                                       
    t0    = time.time()
    funcs = load_functions(found_path)

    if not funcs:
        print(f"No valid functions found in {found_path}")
        print("Check that 'functions' key exists and sbox arrays have 256 elements.")
        return

    print(f"Loaded {len(funcs)} APN functions from {found_path}")
    print()

    #    Signatures                                                 
    print("Computing ortho-derivative signatures...")
    sig_to_ids = defaultdict(list)
    info       = {}

    for rec in funcs:
        aid = rec["id"]
        sb  = rec["sbox"]
        du  = int(sboxU.differential_uniformity(sb))
        deg = int(sboxU.algebraic_degree(sb))
        sig = ortho_signature(sb)
        sh  = sig_hash(sig)

        info[aid] = {**rec, "du": du, "degree": deg,
                     "ortho_sig_hash": sh, "sbox_first8": sb[:8]}
        sig_to_ids[sh].append(aid)
        print(f"  {aid}: DU={du} deg={deg} sig={sh}", flush=True)

    print()
    print(f"Signature buckets: {len(sig_to_ids)}")
    for sh, ids in sorted(sig_to_ids.items(), key=lambda x: -len(x[1])):
        print(f"  {sh}: {len(ids)} function(s)")
    print()

    #    Exact CCZ inside buckets                                    
    all_pairwise = []
    exact_tests  = 0

    if args.skip_intra_bucket:
        print("Skipping exact CCZ (--skip-intra-bucket).")
        print("Quadratic APN: same ortho-sig => same CCZ-class (Yoshiara 2010).")
        for sh, ids in sig_to_ids.items():
            for i in range(1, len(ids)):
                all_pairwise.append((ids[0], ids[i], True))
    else:
        print("Running exact CCZ inside same-signature buckets...")
        for sh, ids in sorted(sig_to_ids.items()):
            if len(ids) <= 1:
                continue
            for i in range(len(ids)):
                for j in range(i + 1, len(ids)):
                    a, b = ids[i], ids[j]
                    print(f"  {a} vs {b} ...", end=" ", flush=True)
                    eq = exact_ccz(info[a]["sbox"], info[b]["sbox"])
                    exact_tests += 1
                    all_pairwise.append((a, b, eq))
                    print(eq, flush=True)

    #    Classes                                                    
    all_ids     = [rec["id"] for rec in funcs]
    groups      = union_find_groups(all_ids, all_pairwise)
    groups.sort(key=lambda g: g[0])

    #    Output                                                     
    elapsed     = time.time() - t0
    lines       = []
    lines.append("FOUND APN CLASSIFICATION USING sboxU")
    lines.append("=" * 80)
    lines.append(f"Input    : {found_path}")
    lines.append(f"Loaded   : {len(funcs)} functions")
    lines.append(f"Buckets  : {len(sig_to_ids)} signature buckets")
    lines.append(f"Classes  : {len(groups)} final CCZ-classes")
    lines.append(f"Exact CCZ: {exact_tests} tests performed")
    lines.append(f"Elapsed  : {elapsed:.1f}s")
    lines.append("")

    lines.append("SIGNATURE BUCKETS")
    lines.append("-" * 80)
    for sh, ids in sorted(sig_to_ids.items(), key=lambda x: -len(x[1])):
        lines.append(f"  {sh}: {len(ids)} functions")
    lines.append("")

    lines.append("FINAL CCZ-CLASSES")
    lines.append("=" * 80)
    class_records = []
    for idx, group in enumerate(groups, 1):
        sh  = info[group[0]]["ortho_sig_hash"]
        deg = info[group[0]]["degree"]
        du  = info[group[0]]["du"]
        lines.append(f"CLASS-{idx:03d}  size={len(group)}  sig={sh}  DU={du}  deg={deg}")
        for aid in group:
            r = info[aid]
            lines.append(f"  {aid}: tag={r['tag']}  batch={r['batch']}  "
                         f"sol={r['sol']}/{r['total']}  first8={r['sbox_first8']}")
        lines.append("")

        class_records.append({
            "class_id": f"CLASS-{idx:03d}",
            "size": len(group),
            "ortho_sig_hash": sh,
            "du": du,
            "degree": deg,
            "members": [
                {"id": aid, "tag": info[aid]["tag"], "batch": info[aid]["batch"],
                 "sol": info[aid]["sol"], "total": info[aid]["total"],
                 "hash": info[aid]["hash"], "first8": info[aid]["sbox_first8"]}
                for aid in group
            ],
        })

    result = {
        "input": str(found_path),
        "functions_loaded": len(funcs),
        "signature_buckets": len(sig_to_ids),
        "final_classes": len(groups),
        "exact_ccz_tests": exact_tests,
        "elapsed_sec": round(elapsed, 2),
        "classes": class_records,
        "pairwise": [{"a": a, "b": b, "ccz": eq} for a, b, eq in all_pairwise],
    }

    Path(out_txt).write_text("\n".join(lines) + "\n", encoding="utf-8")
    Path(out_json).write_text(json.dumps(result, indent=2), encoding="utf-8")

    print()
    print("=" * 55)
    print(f"Functions : {len(funcs)}")
    print(f"Sig buckets: {len(sig_to_ids)}")
    print(f"CCZ classes: {len(groups)}")
    print(f"Report TXT : {out_txt}")
    print(f"Report JSON: {out_json}")

if __name__ == "__main__":
    main()
