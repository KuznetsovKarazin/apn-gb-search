#!/usr/bin/env python3
"""
collect_results.py  (v2 — auto-detects latest batch, filters centers)

Usage:
  # Auto-detect latest batch:
  python src/collect_results.py

  # Specific batch:
  python src/collect_results.py --batch batch_002

  # Specific results dir:
  python src/collect_results.py --results E:/apn-gb-search/batch/batch_002/results
"""
import argparse, json, csv, os, glob, hashlib
from pathlib import Path
from collections import defaultdict

N = 256

def join_wrapped(t):
    return t.replace('\\\r\n','').replace('\\\n','')

def parse_line(line):
    line = line.strip()
    if not line.startswith("APN|"): return None
    parts = line.split("|")
    if len(parts) < 6: return None
    try:
        tag = parts[1]; sol = int(parts[2].split("=")[1])
        i = 3; total = -1
        if parts[i].startswith("total="):
            total = int(parts[i].split("=")[1]); i += 1
        apn  = parts[i].split("=")[1] == "true"
        perm = parts[i+1].split("=")[1] == "true"
        raw  = parts[i+2].split("=",1)[1].strip()
        sbox = list(map(int, raw.split("," if "," in raw else None)))
        if len(sbox) != N: return None
        return {"tag":tag,"sol":sol,"total":total,"apn":apn,"perm":perm,"sbox":sbox}
    except: return None

def sbox_hash(sb):
    return hashlib.sha256(bytes(sb)).hexdigest()[:16]

def load_centers(centers_json):
    """Load center sboxes from apn_centers.json for filtering."""
    if not os.path.exists(centers_json): return set()
    try:
        data = json.load(open(centers_json))
        centers = data.get("centers", [])
        # We only have first 8 bytes in centers file, need to match differently
        # Just return empty — we'll filter by total=1 AND sol=1 heuristic
        return set()
    except: return set()

def find_latest_batch(repo_root):
    """Find the most recently modified batch_NNN directory."""
    batches = sorted(Path(repo_root, "batch").glob("batch_*"))
    batches = [b for b in batches if b.is_dir()]
    if not batches: return None
    # Return latest by name (highest number)
    return batches[-1]

def main():
    repo = Path(__file__).resolve().parent.parent

    parser = argparse.ArgumentParser(description="Collect Magma APN results.")
    parser.add_argument("--batch",   type=str, default=None,
                        help="Batch name, e.g. batch_002 (default: latest)")
    parser.add_argument("--results", type=Path, default=None,
                        help="Direct path to results directory")
    parser.add_argument("--out",     type=Path, default=None,
                        help="Output directory (default: same as results)")
    args = parser.parse_args()

    # Determine results directory
    if args.results:
        results_dir = args.results.resolve()
        batch_dir   = results_dir.parent
    elif args.batch:
        batch_dir   = (repo / "batch" / args.batch).resolve()
        results_dir = batch_dir / "results"
    else:
        batch_dir = find_latest_batch(repo)
        if not batch_dir:
            print("No batch directories found."); return
        results_dir = batch_dir / "results"
        print(f"Auto-detected batch: {batch_dir.name}")

    out_dir = (args.out or results_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    centers_json = str(batch_dir / "apn_centers.json")

    # Find result files
    files = sorted(str(f) for f in results_dir.glob("apn_results_*.txt"))
    if not files:
        print(f"No apn_results_*.txt in: {results_dir}"); return

    print(f"Batch dir   : {batch_dir}")
    print(f"Results dir : {results_dir}")
    print(f"Files found : {len(files)}")
    print()

    # Parse all records
    all_records, seen = [], set()
    for fname in files:
        raw = join_wrapped(open(fname, encoding="utf-8", errors="replace").read())
        for line in raw.splitlines():
            r = parse_line(line)
            if not r: continue
            key = tuple(r["sbox"])
            r["duplicate"] = key in seen
            r["hash"] = sbox_hash(list(key))
            seen.add(key)
            all_records.append(r)

    unique = [r for r in all_records if not r["duplicate"]]

    # Identify centers: sol=1 AND total=1 (only thing found = center itself)
    # Neighbors: total > 1 (multiple solutions in slice)
    centers_rec  = [r for r in unique if r["total"] == 1]
    neighbor_rec = [r for r in unique if r["total"] > 1]

    # Among multi-solution slices: which is the center?
    # The center is always present in the output. We can't reliably identify
    # it without comparing against apn_centers.json sboxes.
    # Strategy: mark ALL as "candidates" — user runs classify_apn.py to sort out.

    print(f"{'='*50}")
    print(f"SUMMARY")
    print(f"{'='*50}")
    print(f"Total records      : {len(all_records)}")
    print(f"Unique sboxes      : {len(unique)}")
    print(f"  total=1 (centers): {len(centers_rec)}")
    print(f"  total>1 (w/ nbrs): {len(neighbor_rec)}")
    print(f"Duplicates removed : {len(all_records)-len(unique)}")
    print()

    # Tags with neighbors
    tags_with_nbrs = defaultdict(list)
    for r in neighbor_rec:
        tags_with_nbrs[r["tag"]].append(r)
    if tags_with_nbrs:
        print("Slices with neighbors (total>1):")
        for tag, recs in sorted(tags_with_nbrs.items()):
            print(f"  {tag}: total={recs[0]['total']}, {len(recs)} records saved")
    else:
        print("No slices with neighbors found.")
    print()

    # Save ALL unique to JSON (compact sbox)
    out_json = str(out_dir / "apn_found.json")
    fns = [{"id": i+1, "tag": r["tag"], "sol": r["sol"], "total": r["total"],
            "perm": r["perm"], "hash": r["hash"], "sbox": r["sbox"]}
           for i, r in enumerate(unique)]
    lines_out = ["{",
                 f'  "batch": "{batch_dir.name}",',
                 f'  "total_unique": {len(unique)},',
                 f'  "with_neighbors": {len(neighbor_rec)},',
                 '  "functions": [']
    for i, fn in enumerate(fns):
        sb = json.dumps(fn["sbox"], separators=(",",":"))
        c  = "," if i < len(fns)-1 else ""
        lines_out += [
            "    {",
            f'      "id": {fn["id"]},',
            f'      "tag": "{fn["tag"]}",',
            f'      "sol": {fn["sol"]}, "total": {fn["total"]},',
            f'      "perm": {str(fn["perm"]).lower()},',
            f'      "hash": "{fn["hash"]}",',
            f'      "sbox": {sb}',
            "    }" + c
        ]
    lines_out += ["  ]", "}"]
    open(out_json, "w").write("\n".join(lines_out)+"\n")
    print(f"All unique  -> {out_json}  ({len(unique)} functions)")

    # Save CSV with neighbors only (the interesting ones)
    out_csv = str(out_dir / "apn_new.csv")
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id","tag","sol","total","perm","hash","sbox_first8"])
        for r in neighbor_rec:
            w.writerow([r.get("id",""),r["tag"],r["sol"],r["total"],r["perm"],
                        r["hash"]," ".join(str(x) for x in r["sbox"][:8])])
    print(f"Neighbors   -> {out_csv}  ({len(neighbor_rec)} functions)")
    print()
    print("Next steps:")
    print("  1. Run CCZ classification:  python src/classify_apn.py")
    print("  2. Check vs Beierle DB:     python src/check_vs_beierle_db.py \\")
    print("                                --found apn_found.json --known data/new_apns.txt \\")
    print("                                --skip-pairwise")

if __name__ == "__main__":
    main()
