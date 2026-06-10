#!/usr/bin/env python3
# apn_registry.py
# Robust APN registry.
#
# Commands:
#   python apn_registry.py collect
#   python apn_registry.py show
#   python apn_registry.py add-seeds
#
# Handles both registry layouts:
#   old: {"functions": [ ... ]}
#   new: {"functions": {"hash": {...}}}
#
# Handles result files with comma/space separated sboxes and multiline backslashes.

import glob
import hashlib
import json
import os
import re
import sys
from datetime import datetime

N = 256
REGISTRY = "apn_found.json"
SEEDS = "all_seeds.json"


def sbox_hash(sb):
    return hashlib.sha256(bytes(sb)).hexdigest()


def short_hash(sb):
    return hashlib.sha256(bytes(sb)).hexdigest()[:16]


def parse_bool(x):
    return str(x).strip().lower() in ("true", "1", "yes")


def read_records_from_file(path):
    records = []
    cur = None
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.rstrip("\n\r")
            if line.startswith("APN|"):
                if cur is not None:
                    records.append(cur)
                cur = line
            elif cur is not None:
                # Do not add spaces: numbers may be split across line-continuation backslash.
                cur += line.strip()
    if cur is not None:
        records.append(cur)
    return records


def parse_record(rec):
    if not rec.startswith("APN|") or "|sbox=" not in rec:
        return None

    head, sbox_text = rec.split("|sbox=", 1)
    sbox_text = sbox_text.replace("\\", "")

    parts = head.split("|")
    if len(parts) < 2:
        return None

    tag = parts[1].strip()
    meta = {}
    for p in parts[2:]:
        if "=" in p:
            k, v = p.split("=", 1)
            meta[k.strip()] = v.strip()

    nums = [int(x) for x in re.findall(r"\d+", sbox_text)]
    if len(nums) != N:
        return None
    if any(x < 0 or x > 255 for x in nums):
        return None

    total_raw = meta.get("total", "")
    total = int(total_raw) if str(total_raw).isdigit() else None

    return {
        "tag": tag,
        "sol": int(meta.get("sol", "0")),
        "total": total,
        "apn": parse_bool(meta.get("apn", "true")),
        "perm": parse_bool(meta.get("perm", "false")),
        "sbox": nums,
        "hash": sbox_hash(nums),
        "short_hash": short_hash(nums),
    }


def empty_registry():
    return {
        "created_at": datetime.utcnow().isoformat() + "Z",
        "seeds": {},
        "functions": {}
    }


def normalize_registry(data):
    if not isinstance(data, dict):
        data = empty_registry()

    data.setdefault("seeds", {})
    funcs = data.get("functions", {})

    # Convert old list layout to hash-keyed dict.
    if isinstance(funcs, list):
        new_funcs = {}
        for i, rec in enumerate(funcs, 1):
            if not isinstance(rec, dict):
                continue
            sb = rec.get("sbox")
            if not isinstance(sb, list) or len(sb) != N:
                continue
            sb = [int(x) for x in sb]
            h = sbox_hash(sb)
            fid = rec.get("id")
            if isinstance(fid, int):
                fid = f"APN-{fid:04d}"
            elif not fid:
                fid = f"APN-{len(new_funcs)+1:04d}"

            new_funcs[h] = {
                "id": fid,
                "tag": rec.get("tag"),
                "sol": rec.get("sol"),
                "total": rec.get("total"),
                "perm": rec.get("perm", False),
                "hash": rec.get("hash", short_hash(sb)),
                "source_file": rec.get("source_file"),
                "ccz": rec.get("ccz", "pending"),
                "sbox": sb,
            }
        data["functions"] = new_funcs

    elif not isinstance(funcs, dict):
        data["functions"] = {}

    # Normalize records inside dict.
    fixed = {}
    for h, rec in data["functions"].items():
        if not isinstance(rec, dict):
            continue
        sb = rec.get("sbox")
        if not isinstance(sb, list) or len(sb) != N:
            continue
        sb = [int(x) for x in sb]
        hh = sbox_hash(sb)
        fid = rec.get("id") or f"APN-{len(fixed)+1:04d}"
        fixed[hh] = {
            "id": fid,
            "tag": rec.get("tag"),
            "sol": rec.get("sol"),
            "total": rec.get("total"),
            "perm": rec.get("perm", False),
            "hash": rec.get("hash", short_hash(sb)),
            "source_file": rec.get("source_file"),
            "ccz": rec.get("ccz", "pending"),
            "sbox": sb,
        }
    data["functions"] = fixed
    return data


def load_registry():
    if not os.path.exists(REGISTRY):
        return empty_registry()
    try:
        data = json.load(open(REGISTRY, "r", encoding="utf-8", errors="replace"))
    except Exception:
        data = empty_registry()
    return normalize_registry(data)


def save_registry(data):
    data = normalize_registry(data)
    data["updated_at"] = datetime.utcnow().isoformat() + "Z"
    with open(REGISTRY, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)


def iter_seed_sboxes_from_all_seeds():
    if not os.path.exists(SEEDS):
        return []

    try:
        data = json.load(open(SEEDS, "r", encoding="utf-8", errors="replace"))
    except Exception:
        return []

    out = []

    def visit(obj, label="seed"):
        if isinstance(obj, dict):
            if "sbox" in obj and isinstance(obj["sbox"], list) and len(obj["sbox"]) == N:
                try:
                    out.append((str(obj.get("id", obj.get("tag", label))), [int(x) for x in obj["sbox"]], obj))
                except Exception:
                    pass
            for k, v in obj.items():
                visit(v, str(k))
        elif isinstance(obj, list):
            if len(obj) == N and all(isinstance(x, int) and 0 <= x <= 255 for x in obj):
                out.append((label, obj, {}))
            else:
                for i, v in enumerate(obj):
                    visit(v, f"{label}_{i}")

    visit(data)
    return out


def add_seeds():
    reg = load_registry()
    before = len(reg["seeds"])

    added = 0
    for label, sb, meta in iter_seed_sboxes_from_all_seeds():
        h = sbox_hash(sb)
        if h in reg["seeds"]:
            continue
        sid = f"SEED-{len(reg['seeds'])+1:04d}"
        reg["seeds"][h] = {
            "id": sid,
            "label": label,
            "hash": h[:16],
            "sbox": sb,
        }
        added += 1

    save_registry(reg)
    print(f"Added {added} seeds (was {before}, now {len(reg['seeds'])})")


def collect():
    reg = load_registry()
    seed_hashes = set(reg.get("seeds", {}).keys())

    files = sorted(glob.glob("apn_results_*.txt"))
    parsed = []
    failed = 0
    seen_records = 0

    for path in files:
        for rec_text in read_records_from_file(path):
            seen_records += 1
            rec = parse_record(rec_text)
            if rec is None:
                failed += 1
            else:
                rec["source_file"] = os.path.basename(path)
                parsed.append(rec)

    existing_hashes = set(reg.get("functions", {}).keys())
    added = 0
    skipped_seed = 0
    skipped_existing = 0

    for rec in parsed:
        h = rec["hash"]
        if h in seed_hashes:
            skipped_seed += 1
            continue
        if h in existing_hashes:
            skipped_existing += 1
            continue

        fid = f"APN-{len(reg['functions'])+1:04d}"
        reg["functions"][h] = {
            "id": fid,
            "tag": rec["tag"],
            "sol": rec["sol"],
            "total": rec["total"],
            "perm": rec["perm"],
            "hash": rec["short_hash"],
            "source_file": rec["source_file"],
            "ccz": "pending",
            "sbox": rec["sbox"],
        }
        existing_hashes.add(h)
        added += 1

    save_registry(reg)
    print(f"Seen APN result records: {seen_records}")
    print(f"Parsed APN result records: {len(parsed)}, failed: {failed}")
    print(f"Skipped seeds: {skipped_seed}")
    print(f"Skipped already registered: {skipped_existing}")
    print(f"Collected {added} new APN functions")
    print(f"Total in registry: {len(reg['functions'])}")


def show():
    reg = load_registry()
    print(f"SEEDS:     {len(reg.get('seeds', {}))} registered")
    print(f"FOUND APN: {len(reg.get('functions', {}))} total")
    print("Found APN:")
    items = sorted(reg.get("functions", {}).values(), key=lambda r: r.get("id", ""))
    for r in items:
        print(f"  {r.get('id')}  seed={r.get('tag')}  sol={r.get('sol')}/{r.get('total')}  perm={r.get('perm')}  ccz={r.get('ccz')}")
        print(f"         hash: {r.get('hash')}")
        print(f"         sbox[0:8]: {r.get('sbox', [])[:8]}")


def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else "show"
    if cmd == "collect":
        collect()
    elif cmd == "show":
        show()
    elif cmd == "add-seeds":
        add_seeds()
    else:
        print("Usage:")
        print("  python apn_registry.py collect")
        print("  python apn_registry.py show")
        print("  python apn_registry.py add-seeds")


if __name__ == "__main__":
    main()
