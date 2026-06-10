#!/usr/bin/env sage -python
# exact_pairwise_6_classes.py
#
# Full exact pairwise CCZ/code-equivalence test for six representatives stored in:
#   apn_6_classes.txt
#
# Input format:
#   one S-box per line:
#   [0,0,0,...,234]
#
# Usage:
#   cd /mnt/e/Magma
#   conda activate sage
#   sage -python exact_pairwise_6_classes.py
#
# Output:
#   exact_pairwise_6_classes_report.txt
#   exact_pairwise_6_classes_report.json

import hashlib
import json
import re
import time

import sboxU

INPUT_TXT = "apn_6_classes.txt"
INPUT_META = "apn_6_classes_meta.json"
OUT_TXT = "exact_pairwise_6_classes_report.txt"
OUT_JSON = "exact_pairwise_6_classes_report.json"

DEFAULT_NAMES = ["APN-0001", "APN-0002", "APN-0004", "APN-0006", "APN-0007", "APN-0010"]


def h16(sb):
    return hashlib.sha256(bytes(sb)).hexdigest()[:16]


def parse_line(line):
    nums = [int(x) for x in re.findall(r"\d+", line)]
    if len(nums) != 256:
        return None
    if any(x < 0 or x > 255 for x in nums):
        return None
    return nums


def load_sboxes(path):
    sboxes = []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line_no, line in enumerate(f, 1):
            sb = parse_line(line)
            if sb is None:
                raise RuntimeError(f"Cannot parse S-box at line {line_no}; expected 256 integers")
            sboxes.append(sb)
    return sboxes


def load_names(n):
    try:
        meta = json.load(open(INPUT_META, "r", encoding="utf-8", errors="replace"))
        reps = meta.get("representatives", [])
        names = []
        for r in reps:
            names.append(f"{r.get('class_id')}:{r.get('apn_id')}")
        if len(names) == n:
            return names
    except Exception:
        pass
    return DEFAULT_NAMES[:n]


def basic_info(name, sb):
    return {
        "name": name,
        "hash": h16(sb),
        "du": int(sboxU.differential_uniformity(sb)),
        "degree": int(sboxU.algebraic_degree(sb)),
        "perm": len(set(sb)) == 256,
        "first8": sb[:8],
    }


def main():
    t0 = time.time()

    sboxes = load_sboxes(INPUT_TXT)
    names = load_names(len(sboxes))

    if len(sboxes) != 6:
        raise RuntimeError(f"Expected exactly 6 S-boxes in {INPUT_TXT}, got {len(sboxes)}")

    lines = []
    lines.append("EXACT PAIRWISE CCZ CHECK FOR SIX APN CLASS REPRESENTATIVES")
    lines.append("=" * 90)
    lines.append(f"Input: {INPUT_TXT}")
    lines.append("")

    print("Loaded 6 representatives:")
    infos = []
    lines.append("REPRESENTATIVES")
    lines.append("-" * 90)

    for name, sb in zip(names, sboxes):
        info = basic_info(name, sb)
        infos.append(info)
        msg = (
            f"{name}: hash={info['hash']} DU={info['du']} "
            f"degree={info['degree']} perm={info['perm']} first8={info['first8']}"
        )
        print("  " + msg)
        lines.append(msg)

    lines.append("")
    lines.append("PAIRWISE EXACT CCZ")
    lines.append("-" * 90)
    print("")
    print("Running exact pairwise CCZ tests...")

    results = []
    for i in range(len(sboxes)):
        for j in range(i + 1, len(sboxes)):
            ni, nj = names[i], names[j]
            print(f"Testing {ni} vs {nj} ...", flush=True)
            ts = time.time()
            try:
                eq = bool(sboxU.are_ccz_equivalent_from_code(sboxes[i], sboxes[j]))
                err = None
            except Exception as e:
                eq = None
                err = repr(e)
            dt = time.time() - ts

            msg = f"{ni} vs {nj}: ccz={eq} time={dt:.2f}s"
            if err:
                msg += f" error={err}"
            print("  " + msg, flush=True)
            lines.append(msg)

            results.append({
                "a": ni,
                "b": nj,
                "ccz": eq,
                "time_sec": dt,
                "error": err,
            })

    any_equiv = any(r["ccz"] is True for r in results)
    all_false = all(r["ccz"] is False for r in results)

    lines.append("")
    lines.append("SUMMARY")
    lines.append("=" * 90)
    lines.append(f"Representatives tested: {len(sboxes)}")
    lines.append(f"Pairwise tests: {len(results)}")
    lines.append(f"Any equivalent pair: {any_equiv}")
    lines.append(f"All pairwise non-equivalent: {all_false}")
    lines.append(f"Elapsed seconds: {time.time() - t0:.2f}")

    report = {
        "input": INPUT_TXT,
        "representatives": infos,
        "pairwise_results": results,
        "any_equivalent_pair": any_equiv,
        "all_pairwise_non_equivalent": all_false,
        "elapsed_sec": time.time() - t0,
    }

    with open(OUT_TXT, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

    with open(OUT_JSON, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)

    print("")
    print("DONE")
    print(f"Any equivalent pair: {any_equiv}")
    print(f"All pairwise non-equivalent: {all_false}")
    print(f"Report TXT : {OUT_TXT}")
    print(f"Report JSON: {OUT_JSON}")


if __name__ == "__main__":
    main()
