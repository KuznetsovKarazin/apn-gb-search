#!/usr/bin/env python3
"""
merge_batches.py
================
Merges apn_found.json files from multiple batches into one
global deduplicated collection.

Usage:
  # Merge all batches automatically:
  python src/merge_batches.py

  # Specific batches:
  python src/merge_batches.py --batches batch_002 batch_003 batch_004

  # Output to custom dir:
  python src/merge_batches.py --out results/merged

Output:
  results/merged/apn_all.json      - all unique functions
  results/merged/apn_summary.txt   - human-readable summary
"""

import argparse, json, hashlib, glob, os
from pathlib import Path
from collections import defaultdict

N = 256

def sbox_hash(sb):
    return hashlib.sha256(bytes(sb)).hexdigest()[:16]

def load_batch(found_json):
    data = json.load(open(found_json, encoding='utf-8', errors='replace'))
    fns = data.get("functions", [])
    if isinstance(fns, dict): fns = list(fns.values())
    result = []
    for r in fns:
        sb = [int(x) for x in r.get("sbox", [])]
        if len(sb) != N: continue
        result.append({
            "tag":   r.get("tag", ""),
            "sol":   r.get("sol", -1),
            "total": r.get("total", -1),
            "perm":  r.get("perm", False),
            "hash":  r.get("hash") or sbox_hash(sb),
            "sbox":  sb,
            "batch": data.get("batch", "?"),
        })
    return result

def find_all_batches(repo_root):
    batches = sorted(Path(repo_root, "batch").glob("batch_*"))
    result = []
    for b in batches:
        f = b / "results" / "apn_found.json"
        if f.exists():
            result.append((b.name, str(f)))
    return result

def main():
    repo = Path(__file__).resolve().parent.parent

    parser = argparse.ArgumentParser(description="Merge APN results from multiple batches.")
    parser.add_argument("--batches", nargs="+", default=None,
                        help="Batch names to merge (default: all)")
    parser.add_argument("--out", type=Path, default=None,
                        help="Output directory (default: results/merged)")
    args = parser.parse_args()

    # Find batch files
    if args.batches:
        batch_files = []
        for b in args.batches:
            f = repo / "batch" / b / "results" / "apn_found.json"
            if f.exists():
                batch_files.append((b, str(f)))
            else:
                print(f"  WARNING: {f} not found, skipping")
    else:
        batch_files = find_all_batches(repo)

    if not batch_files:
        print("No apn_found.json files found."); return

    print(f"Merging {len(batch_files)} batch(es):")
    for name, path in batch_files:
        print(f"  {name}: {path}")
    print()

    # Merge with deduplication
    all_records = []
    seen_hashes = set()
    batch_stats = {}

    for batch_name, path in batch_files:
        records = load_batch(path)
        n_before = len(seen_hashes)
        new_in_batch = 0
        dup_in_batch = 0
        rich_in_batch = 0

        for r in records:
            h = r["hash"]
            if h not in seen_hashes:
                seen_hashes.add(h)
                all_records.append(r)
                new_in_batch += 1
                if r["total"] > 1:
                    rich_in_batch += 1
            else:
                dup_in_batch += 1

        batch_stats[batch_name] = {
            "total": len(records),
            "new_unique": new_in_batch,
            "duplicates": dup_in_batch,
            "with_neighbors": rich_in_batch,
        }
        print(f"  {batch_name}: {len(records)} records -> "
              f"{new_in_batch} new unique, {dup_in_batch} duplicates")

    total = len(all_records)
    with_nbrs = [r for r in all_records if r["total"] > 1]
    centers   = [r for r in all_records if r["total"] == 1]
    print()
    print(f"{'='*50}")
    print(f"MERGED TOTAL")
    print(f"{'='*50}")
    print(f"Unique functions   : {total}")
    print(f"  total=1 (centers): {len(centers)}")
    print(f"  total>1 (w/ nbrs): {len(with_nbrs)}")

    # Cross-batch duplicates (same function found in multiple batches)
    cross_dups = sum(s["duplicates"] for s in batch_stats.values())
    if cross_dups:
        print(f"Cross-batch dups   : {cross_dups} (same APN found in multiple batches)")
    print()

    # Save output
    out_dir = (args.out or repo / "results" / "merged").resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    # Compact JSON
    out_json = str(out_dir / "apn_all.json")
    fns_out = []
    for i, r in enumerate(all_records, 1):
        fns_out.append({
            "id": i, "batch": r["batch"], "tag": r["tag"],
            "sol": r["sol"], "total": r["total"],
            "perm": r["perm"], "hash": r["hash"], "sbox": r["sbox"]
        })

    lines = ["{",
             f'  "merged_batches": {json.dumps([b for b,_ in batch_files])},',
             f'  "total_unique": {total},',
             f'  "with_neighbors": {len(with_nbrs)},',
             '  "functions": [']
    for i, fn in enumerate(fns_out):
        sb  = json.dumps(fn["sbox"], separators=(",",":"))
        c   = "," if i < len(fns_out)-1 else ""
        lines += [
            "    {",
            f'      "id": {fn["id"]}, "batch": "{fn["batch"]}",',
            f'      "tag": "{fn["tag"]}", "sol": {fn["sol"]}, "total": {fn["total"]},',
            f'      "perm": {str(fn["perm"]).lower()}, "hash": "{fn["hash"]}",',
            f'      "sbox": {sb}',
            "    }" + c
        ]
    lines += ["  ]", "}"]
    open(out_json, "w").write("\n".join(lines) + "\n")
    print(f"Merged JSON -> {out_json}")

    # Summary TXT
    out_txt = str(out_dir / "apn_summary.txt")
    with open(out_txt, "w") as f:
        f.write("MERGED APN COLLECTION\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Total unique functions : {total}\n")
        f.write(f"  Centers (total=1)    : {len(centers)}\n")
        f.write(f"  With neighbors (>1)  : {len(with_nbrs)}\n\n")
        f.write("Per-batch stats:\n")
        for bn, s in batch_stats.items():
            f.write(f"  {bn}: {s['total']} records, "
                    f"{s['new_unique']} unique, "
                    f"{s['duplicates']} cross-dups, "
                    f"{s['with_neighbors']} rich\n")
        f.write("\n")
        f.write("Next steps:\n")
        f.write("  # Classify merged collection:\n")
        f.write(f"  python src/classify_apn_standalone.py --found {out_json}\n\n")
        f.write("  # Check vs Beierle 2025 (with sboxU):\n")
        f.write("  conda activate sage\n")
        f.write(f"  python src/check_vs_beierle_db.py --found {out_json} "
                f"--known data/new_apns.txt --skip-pairwise\n")
    print(f"Summary TXT -> {out_txt}")
    print()
    print("Next steps:")
    print(f"  python src/classify_apn_standalone.py --found {out_json}")
    print(f"  conda activate sage")
    print(f"  python src/check_vs_beierle_db.py --found {out_json} "
          f"--known data/new_apns.txt --skip-pairwise")

if __name__ == "__main__":
    main()
