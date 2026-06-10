#!/usr/bin/env python3
# make_6_class_representatives_txt.py
#
# Creates apn_6_classes.txt in the same one-line-per-S-box format as new_apns.txt.
#
# Representatives:
#   CLASS-001 -> APN-0001
#   CLASS-002 -> APN-0002
#   CLASS-003 -> APN-0004
#   CLASS-004 -> APN-0006
#   CLASS-005 -> APN-0007
#   CLASS-006 -> APN-0010
#
# Usage:
#   cd E:\Magma
#   python make_6_class_representatives_txt.py
#
# Output:
#   apn_6_classes.txt
#   apn_6_classes_meta.json

import json
import hashlib

INPUT = "apn_found.json"
OUT_TXT = "apn_6_classes.txt"
OUT_JSON = "apn_6_classes_meta.json"

REPS = [
    ("CLASS-001", "APN-0001"),
    ("CLASS-002", "APN-0002"),
    ("CLASS-003", "APN-0004"),
    ("CLASS-004", "APN-0006"),
    ("CLASS-005", "APN-0007"),
    ("CLASS-006", "APN-0010"),
]


def h16(sb):
    return hashlib.sha256(bytes(sb)).hexdigest()[:16]


def load_functions(path):
    data = json.load(open(path, "r", encoding="utf-8", errors="replace"))
    funcs = data.get("functions", {})

    if isinstance(funcs, dict):
        records = list(funcs.values())
    elif isinstance(funcs, list):
        records = funcs
    else:
        records = []

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
        out[aid] = {
            "id": aid,
            "tag": rec.get("tag") or rec.get("seed"),
            "sol": rec.get("sol"),
            "total": rec.get("total"),
            "perm": rec.get("perm", False),
            "hash": h16(sb),
            "sbox": sb,
        }

    return out


def main():
    funcs = load_functions(INPUT)

    missing = [aid for _, aid in REPS if aid not in funcs]
    if missing:
        raise RuntimeError(f"Missing APN representatives in {INPUT}: {missing}")

    selected = []
    with open(OUT_TXT, "w", encoding="utf-8", newline="\n") as f:
        for class_id, aid in REPS:
            rec = funcs[aid]
            sb = rec["sbox"]
            f.write("[" + ",".join(str(x) for x in sb) + "]\n")
            selected.append({
                "class_id": class_id,
                "apn_id": aid,
                "tag": rec.get("tag"),
                "sol": rec.get("sol"),
                "total": rec.get("total"),
                "hash": h16(sb),
                "first8": sb[:8],
            })

    out = {
        "description": "Six representatives of the six APN classes found",
        "format": "one S-box per line, same style as new_apns.txt",
        "txt_file": OUT_TXT,
        "representatives": selected,
    }

    with open(OUT_JSON, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)

    print(f"Created {OUT_TXT}")
    print(f"Created {OUT_JSON}")
    print("")
    for r in selected:
        print(f"{r['class_id']} -> {r['apn_id']}: hash={r['hash']} tag={r['tag']} sol={r['sol']}/{r['total']} first8={r['first8']}")


if __name__ == "__main__":
    main()
