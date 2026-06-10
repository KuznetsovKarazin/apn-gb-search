#!/usr/bin/env python3
"""
check_class22_membership.py
===========================

Check whether S-boxes satisfy the Class-22 self-equivalence:

    F(A x) = A F(x)    for all x in F_2^8,

where A = block_diag(C4, C4), and C4 is the same companion matrix as in
find_apn_centers.py:

    C4 = make_companion([1,1,1,1], 4)

Run examples:

  python src/check_class22_membership.py --input batch/batch_003/results/apn_found.json

  python src/check_class22_membership.py ^
      --input batch/batch_003/results/apn_found.json ^
      --out-csv batch/batch_003/results/class22_membership.csv ^
      --out-json batch/batch_003/results/class22_membership.json

  python src/check_class22_membership.py --input batch/batch_003/results/apn_new.csv

  python src/check_class22_membership.py --input "batch/batch_003/results/apn_results_*.txt"
"""

import argparse
import csv
import glob
import hashlib
import json
import re
import sys
from pathlib import Path


N = 256
n = 8


def make_companion(poly, sz):
    M = [[0] * sz for _ in range(sz)]
    for i in range(sz):
        M[i][sz - 1] = poly[i] & 1
    for i in range(1, sz):
        M[i][i - 1] = 1
    return M


def block_diag(blocks):
    sz = sum(len(b) for b in blocks)
    M = [[0] * sz for _ in range(sz)]
    off = 0
    for b in blocks:
        k = len(b)
        for i in range(k):
            for j in range(k):
                M[off + i][off + j] = b[i][j] & 1
        off += k
    return M


C4 = make_companion([1, 1, 1, 1], 4)
A22 = block_diag([C4, C4])


def mat_vec_mul_int(A, x):
    y = 0
    for j, row in enumerate(A):
        bit = 0
        for i, aji in enumerate(row):
            if aji:
                bit ^= (x >> i) & 1
        if bit:
            y |= (1 << j)
    return y


A22_TABLE = [mat_vec_mul_int(A22, x) for x in range(N)]


def parse_int_list(text):
    nums = [int(z) for z in re.findall(r"-?\d+", text)]
    if len(nums) != N:
        raise ValueError(f"Expected 256 integers, got {len(nums)}")
    if any(v < 0 or v > 255 for v in nums):
        raise ValueError("S-box values must be in [0,255]")
    return nums


def normalize_sbox(obj):
    if isinstance(obj, list):
        if len(obj) == N and all(isinstance(x, int) for x in obj):
            return obj
        return None

    if isinstance(obj, dict):
        for key in ("sbox", "table", "lookup", "values", "truth_table"):
            if key in obj:
                val = obj[key]
                if isinstance(val, list):
                    return parse_int_list(" ".join(map(str, val)))
                if isinstance(val, str):
                    return parse_int_list(val)
        return None

    if isinstance(obj, str):
        return parse_int_list(obj)

    return None


def iter_json_records(path):
    data = json.load(open(path, "r", encoding="utf-8"))

    sb = normalize_sbox(data)
    if sb is not None:
        yield {"source": str(path), "label": path.stem, "sbox": sb}
        return

    if isinstance(data, dict):
        for key in ("records", "functions", "found", "found_apn", "apn",
                    "sboxes", "items", "data", "centers"):
            if key in data and isinstance(data[key], list):
                for i, rec in enumerate(data[key], 1):
                    sb = normalize_sbox(rec)
                    if sb is not None:
                        if isinstance(rec, dict):
                            label = (rec.get("id") or rec.get("label") or
                                     rec.get("tag") or rec.get("name") or
                                     rec.get("idx") or str(i))
                        else:
                            label = str(i)
                        yield {"source": str(path), "label": str(label), "sbox": sb}
                return

        for key, val in data.items():
            sb = normalize_sbox(val)
            if sb is not None:
                yield {"source": str(path), "label": str(key), "sbox": sb}
                return

    if isinstance(data, list):
        for i, rec in enumerate(data, 1):
            sb = normalize_sbox(rec)
            if sb is not None:
                if isinstance(rec, dict):
                    label = (rec.get("id") or rec.get("label") or rec.get("tag") or
                             rec.get("name") or rec.get("idx") or str(i))
                else:
                    label = str(i)
                yield {"source": str(path), "label": str(label), "sbox": sb}


def iter_csv_records(path):
    with open(path, "r", encoding="utf-8", newline="") as f:
        sample = f.read(4096)
        f.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters=",;\t|")
        reader = csv.DictReader(f, dialect=dialect)
        fields = reader.fieldnames or []

        sbox_col = None
        for c in fields:
            if c.lower() in ("sbox", "table", "lookup", "values", "truth_table"):
                sbox_col = c
                break

        if sbox_col is None:
            raise ValueError(f"No S-box column found in CSV {path}. Columns: {fields}")

        for i, row in enumerate(reader, 1):
            sb = parse_int_list(row[sbox_col])
            label = row.get("id") or row.get("label") or row.get("tag") or row.get("name") or str(i)
            yield {"source": str(path), "label": str(label), "sbox": sb}


def iter_txt_records(path):
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        buf = ""
        active = False
        line_no = 0

        def emit_from_record(rec, label_hint):
            m = re.search(r"sbox=([0-9,\s\\-]+)", rec)
            if not m:
                return None
            sb = parse_int_list(m.group(1))
            tag = label_hint
            mt = re.search(r"APN\|([^|]+)\|", rec)
            if mt:
                tag = mt.group(1)
            ms = re.search(r"\|sol=([^|]+)", rec)
            if ms:
                tag = f"{tag}:sol={ms.group(1)}"
            return {"source": str(path), "label": tag, "sbox": sb}

        for raw in f:
            line_no += 1
            line = raw.strip()
            if line.startswith("APN|"):
                if active and buf:
                    rec = emit_from_record(buf, f"line{line_no}")
                    if rec:
                        yield rec
                buf = line
                active = True
            elif active:
                buf += " " + line.replace("\\", "")

        if active and buf:
            rec = emit_from_record(buf, "last")
            if rec:
                yield rec


def iter_records(patterns):
    for pattern in patterns:
        matches = sorted(glob.glob(pattern)) or [pattern]
        for p0 in matches:
            path = Path(p0)
            if not path.exists():
                print(f"WARNING: missing file: {path}", file=sys.stderr)
                continue
            suffix = path.suffix.lower()
            if suffix == ".json":
                yield from iter_json_records(path)
            elif suffix == ".csv":
                yield from iter_csv_records(path)
            else:
                yield from iter_txt_records(path)


def sbox_hash(sb):
    return hashlib.sha256(bytes(sb)).hexdigest()[:16]


def check_class22(sb, max_examples=5):
    violations = []
    for x in range(N):
        ax = A22_TABLE[x]
        lhs = sb[ax]             # F(Ax)
        rhs = A22_TABLE[sb[x]]   # A(Fx)
        if lhs != rhs:
            violations.append((x, ax, lhs, sb[x], rhs))
            if len(violations) >= max_examples:
                break
    return len(violations) == 0, violations


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", nargs="+", required=True,
                    help="Input JSON/CSV/TXT files or glob patterns.")
    ap.add_argument("--out-csv", default=None)
    ap.add_argument("--out-json", default=None)
    ap.add_argument("--max-examples", type=int, default=5)
    args = ap.parse_args()

    rows = []
    total = inside = outside = bad = 0
    seen_hashes = set()

    for rec in iter_records(args.input):
        total += 1
        label = rec["label"]
        source = rec["source"]

        try:
            sb = rec["sbox"]
            h = sbox_hash(sb)
            dup = h in seen_hashes
            seen_hashes.add(h)

            ok, violations = check_class22(sb, args.max_examples)
            status = "IN_CLASS22" if ok else "OUT_CLASS22"
            inside += int(ok)
            outside += int(not ok)

            row = {
                "index": total,
                "label": label,
                "source": source,
                "hash16": h,
                "duplicate_hash": dup,
                "status": status,
                "num_violation_examples": len(violations),
                "violations": violations,
            }
            rows.append(row)

            if total <= 20 or not ok:
                print(f"{total:5d} {label:30s} {h} {status}")

        except Exception as e:
            bad += 1
            rows.append({
                "index": total,
                "label": label,
                "source": source,
                "hash16": "",
                "duplicate_hash": False,
                "status": "PARSE_OR_CHECK_ERROR",
                "error": str(e),
                "violations": [],
            })
            print(f"{total:5d} {label:30s} ERROR: {e}")

    print()
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Records checked : {total}")
    print(f"Unique hashes   : {len(seen_hashes)}")
    print(f"IN_CLASS22      : {inside}")
    print(f"OUT_CLASS22     : {outside}")
    print(f"Errors          : {bad}")

    if outside:
        print()
        print("First OUT_CLASS22 examples:")
        k = 0
        for r in rows:
            if r["status"] == "OUT_CLASS22":
                k += 1
                print(f"  {r['index']} {r['label']} hash={r['hash16']}")
                for v in r["violations"][:args.max_examples]:
                    x, ax, fax, fx, afx = v
                    print(f"    x={x:3d} A(x)={ax:3d} F(Ax)={fax:3d} F(x)={fx:3d} A(Fx)={afx:3d}")
                if k >= 10:
                    break

    if args.out_csv:
        out = Path(args.out_csv)
        out.parent.mkdir(parents=True, exist_ok=True)
        with open(out, "w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=["index", "label", "source", "hash16",
                            "duplicate_hash", "status", "num_violation_examples"]
            )
            writer.writeheader()
            for r in rows:
                writer.writerow({
                    "index": r.get("index"),
                    "label": r.get("label"),
                    "source": r.get("source"),
                    "hash16": r.get("hash16"),
                    "duplicate_hash": r.get("duplicate_hash"),
                    "status": r.get("status"),
                    "num_violation_examples": r.get("num_violation_examples", 0),
                })
        print(f"CSV report: {out}")

    if args.out_json:
        out = Path(args.out_json)
        out.parent.mkdir(parents=True, exist_ok=True)
        with open(out, "w", encoding="utf-8") as f:
            json.dump({
                "matrix": "A22 = block_diag(C4,C4), C4=make_companion([1,1,1,1],4)",
                "criterion": "F(Ax) == A(F(x)) for all x in 0..255",
                "summary": {
                    "records_checked": total,
                    "unique_hashes": len(seen_hashes),
                    "in_class22": inside,
                    "out_class22": outside,
                    "errors": bad,
                },
                "records": rows,
            }, f, indent=2)
        print(f"JSON report: {out}")


if __name__ == "__main__":
    main()
