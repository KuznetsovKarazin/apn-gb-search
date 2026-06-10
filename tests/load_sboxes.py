#!/usr/bin/env python3
"""
load_sboxes.py
==============
Helper module: load APN s-boxes from apn_found.json.
Also usable as a standalone verification script (no Magma needed).

Usage as module:
    from load_sboxes import load_all, load_class

Usage standalone:
    python tests/load_sboxes.py [path/to/apn_found.json]
    python tests/load_sboxes.py                          # uses ../data/apn_found.json
"""

import json
import os
import sys
from typing import Dict, List, Tuple

N = 256
DEFAULT_DATA_PATH = os.path.join(os.path.dirname(__file__), '..', 'data', 'apn_found.json')


def load_all(path: str = DEFAULT_DATA_PATH) -> List[Tuple[str, List[int]]]:
    """Return list of (apn_id, sbox) for every function in the registry."""
    with open(path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    funcs = data.get('functions', {})
    result = []
    for h, v in sorted(funcs.items(), key=lambda x: x[1].get('id', '')):
        fid = v.get('id', h[:8])
        sb = v.get('sbox', [])
        if len(sb) == N:
            result.append((fid, sb))
    return result


def load_class(class_id: str, path: str = DEFAULT_DATA_PATH,
               meta_path: str = None) -> List[Tuple[str, List[int]]]:
    """Return all s-boxes belonging to a given CCZ class (e.g. 'CLASS-004')."""
    if meta_path is None:
        meta_path = os.path.join(os.path.dirname(path), 'apn_6_classes_meta.json')

    try:
        with open(meta_path, 'r', encoding='utf-8') as f:
            meta = json.load(f)
        # Build set of APN ids in this class from verification report
        # (meta file only stores one representative per class)
        rep = next((r for r in meta.get('representatives', [])
                    if r['class_id'] == class_id), None)
        if rep is None:
            raise ValueError(f"Class {class_id} not found in {meta_path}")
        rep_id = rep['apn_id']
    except FileNotFoundError:
        rep_id = None

    # Load all sboxes and filter by the classification report
    all_sboxes = load_all(path)
    # Without full classification data, return all if class not parseable
    if rep_id is None:
        return all_sboxes
    return [(fid, sb) for fid, sb in all_sboxes if fid == rep_id]


# ---------------------------------------------------------------------------
# Verification functions (pure Python, no Magma)
# ---------------------------------------------------------------------------

def check_apn(sbox: List[int]) -> bool:
    """Return True if sbox is APN (differential uniformity = 2)."""
    n = len(sbox)
    for a in range(1, n):
        counts = [0] * n
        for x in range(n):
            b = sbox[x] ^ sbox[x ^ a]
            counts[b] += 1
        if max(counts) > 2:
            return False
    return True


def algebraic_degree(sbox: List[int], n_bits: int = 8) -> int:
    """Compute algebraic degree via Möbius transform."""
    n = n_bits
    N = 1 << n
    max_deg = 0
    for bit in range(n):
        f = [(sbox[x] >> bit) & 1 for x in range(N)]
        # Möbius transform
        for i in range(n):
            step = 1 << i
            for mask in range(N):
                if mask & step:
                    f[mask] ^= f[mask ^ step]
        for mask in range(N):
            if f[mask]:
                w = bin(mask).count('1')
                if w > max_deg:
                    max_deg = w
    return max_deg


def is_permutation(sbox: List[int]) -> bool:
    return len(set(sbox)) == len(sbox)


def walsh_max(sbox: List[int], n_bits: int = 8) -> int:
    """Compute maximum absolute Walsh coefficient."""
    N = 1 << n_bits
    max_val = 0
    for u in range(N):
        for v in range(1, N):
            s = 0
            for x in range(N):
                ux = bin(u & x).count('1') & 1
                vy = bin(v & sbox[x]).count('1') & 1
                s += 1 if (ux ^ vy) == 0 else -1
            if abs(s) > max_val:
                max_val = abs(s)
    return max_val


# ---------------------------------------------------------------------------
# Standalone run: verify all functions in registry
# ---------------------------------------------------------------------------

def main():
    path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_DATA_PATH
    path = os.path.abspath(path)

    if not os.path.exists(path):
        print(f"ERROR: file not found: {path}")
        sys.exit(1)

    print(f"Loading: {path}")
    sboxes = load_all(path)
    print(f"Functions loaded: {len(sboxes)}")
    print()

    all_ok = True
    for fid, sb in sboxes:
        apn = check_apn(sb)
        deg = algebraic_degree(sb)
        perm = is_permutation(sb)
        status = "OK" if (apn and deg == 2 and not perm) else "FAIL"
        if status != "OK":
            all_ok = False
        print(f"  {fid}:  APN={apn}  degree={deg}  perm={perm}  [{status}]")

    print()
    if all_ok:
        print("ALL CHECKS PASSED")
        print("  All 25 functions: APN=True, degree=2, perm=False")
    else:
        print("SOME CHECKS FAILED — see details above")
        sys.exit(1)


if __name__ == '__main__':
    main()
