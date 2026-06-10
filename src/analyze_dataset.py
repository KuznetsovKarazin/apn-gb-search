#!/usr/bin/env python3
"""
src/analyze_dataset.py
======================
Analyze external APN datasets against our search subspace.

Tasks:
  1. Count how many functions in a dataset lie in the Class-22 self-equivalence
     subspace (i.e., satisfy F∘A = A∘F for the Class-22 matrix A).
  2. Compute ortho-derivative signatures and compare against our found classes.
  3. Estimate coverage: what fraction of the dataset could our method find.

Usage:
  # Scan all of new_apns.txt for Class-22 membership (streaming, ~30 sec):
  python src/analyze_dataset.py --file data/new_apns.txt --task membership

  # Scan and also compute ortho-signatures for Class-22 members:
  python src/analyze_dataset.py --file data/new_apns.txt --task signatures

  # Quick sample (first N lines):
  python src/analyze_dataset.py --file data/new_apns.txt --task membership --max-lines 10000

  # Full analysis of small dataset (apn_8bit.txt, 8157 functions):
  python src/analyze_dataset.py --file data/apn_8bit.txt --task full

Input format: one s-box per line as [v0,v1,...,v255]  (Beierle / sboxU format)

Output: data/analysis/<basename>_membership.json  (and .txt summary)
"""

import argparse
import json
import sys
import time
from pathlib import Path
from collections import defaultdict

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
_REPO = Path(__file__).resolve().parent.parent
DATA_DIR    = _REPO / "data"
FOUND_JSON  = DATA_DIR / "apn_found.json"
VERIF_DIR   = DATA_DIR / "verification"
ANALYSIS_DIR = DATA_DIR / "analysis"

N = 256
n = 8

# ---------------------------------------------------------------------------
# Build Class-22 matrix A and precompute lookup Ax[x] = A*x
# ---------------------------------------------------------------------------

def make_companion(poly, size):
    M = [[0]*size for _ in range(size)]
    for i in range(size):
        M[i][size-1] = poly[i]
    for i in range(1, size):
        M[i][i-1] = 1
    return M

def block_diag(B1, B2):
    n1, n2 = len(B1), len(B2)
    sz = n1 + n2
    M = [[0]*sz for _ in range(sz)]
    for i in range(n1):
        for j in range(n1):
            M[i][j] = B1[i][j]
    for i in range(n2):
        for j in range(n2):
            M[n1+i][n1+j] = B2[i][j]
    return M

def build_linear_map_lookup(M, n=8):
    """Precompute lookup[x] = M*x (GF(2)^n -> GF(2)^n integer representation)."""
    N = 1 << n
    lookup = [0] * N
    for x in range(N):
        bits = [(x >> i) & 1 for i in range(n)]
        val = 0
        for i in range(n):
            b = 0
            for j in range(n):
                b ^= M[i][j] * bits[j]
            val |= (b & 1) << i
        lookup[x] = val
    return lookup

# Build once at module load
_C4   = make_companion([1, 1, 1, 1], 4)
_A22  = block_diag(_C4, _C4)
_Ax22 = build_linear_map_lookup(_A22)   # Ax22[x] = A22 * x

# Also build Class-30 (for future use)
_C2   = make_companion([1, 0], 2)
_I2   = [[1, 0], [0, 1]]
_A30  = block_diag(block_diag(_I2, _C2), block_diag(_C2, _C2))
_Ax30 = build_linear_map_lookup(_A30)

CLASS_LOOKUPS = {
    "class22": _Ax22,
    "class30": _Ax30,
}

# ---------------------------------------------------------------------------
# Core checks
# ---------------------------------------------------------------------------

def check_self_equiv(sbox, Ax_lookup):
    """Return True if F(A*x) == A*F(x) for all x. O(256) operations."""
    for x in range(N):
        if sbox[Ax_lookup[x]] != Ax_lookup[sbox[x]]:
            return False
    return True

def check_apn(sbox):
    """Return True if differential uniformity = 2."""
    for a in range(1, N):
        counts = [0] * N
        for x in range(N):
            counts[sbox[x] ^ sbox[x ^ a]] += 1
        if max(counts) > 2:
            return False
    return True

def algebraic_degree(sbox):
    """Compute algebraic degree via Mobius transform."""
    f_all = []
    for bit in range(n):
        f = [(sbox[x] >> bit) & 1 for x in range(N)]
        for i in range(n):
            step = 1 << i
            for mask in range(N):
                if mask & step:
                    f[mask] ^= f[mask ^ step]
        d = max((bin(mask).count('1') for mask in range(N) if f[mask]), default=0)
        f_all.append(d)
    return max(f_all)

def image_size(sbox):
    return len(set(sbox))

# ---------------------------------------------------------------------------
# File parsing
# ---------------------------------------------------------------------------

def parse_line(line):
    """Parse a line like [0,1,2,...,255] into a list of 256 ints. Returns None on error."""
    line = line.strip()
    if not line or not line.startswith('['):
        return None
    try:
        vals = json.loads(line)
        if isinstance(vals, list) and len(vals) == N and all(0 <= v < N for v in vals):
            return vals
    except Exception:
        pass
    return None

def iter_sboxes(path, max_lines=None):
    """Yield (line_number, sbox) from a file."""
    with open(path, 'r', encoding='utf-8', errors='replace') as f:
        for i, line in enumerate(f):
            if max_lines and i >= max_lines:
                break
            sb = parse_line(line)
            if sb is not None:
                yield i + 1, sb

# ---------------------------------------------------------------------------
# Tasks
# ---------------------------------------------------------------------------

def task_membership(path, max_lines=None, classes=("class22",), verbose=True):
    """
    For each function in the file, check which self-equivalence subspaces
    it belongs to. Fast: ~9 microseconds per function per class.
    """
    counters = {cls: 0 for cls in classes}
    total = 0
    members = {cls: [] for cls in classes}
    t0 = time.time()
    report_every = 100_000

    for line_no, sbox in iter_sboxes(path, max_lines):
        total += 1
        for cls in classes:
            if check_self_equiv(sbox, CLASS_LOOKUPS[cls]):
                counters[cls] += 1
                members[cls].append({'line': line_no, 'sbox_first8': sbox[:8]})

        if verbose and total % report_every == 0:
            elapsed = time.time() - t0
            rate = total / elapsed
            pct = {cls: counters[cls]/total*100 for cls in classes}
            print(f"  {total:>9,}  elapsed={elapsed:.0f}s  rate={rate:.0f}/s  "
                  + "  ".join(f"{cls}={counters[cls]}({pct[cls]:.4f}%)" for cls in classes))

    elapsed = time.time() - t0
    result = {
        'file': str(path),
        'total_lines': total,
        'elapsed_sec': round(elapsed, 2),
        'rate_per_sec': round(total / elapsed, 0),
        'classes': {}
    }
    for cls in classes:
        pct = counters[cls] / total * 100 if total > 0 else 0.0
        result['classes'][cls] = {
            'count': counters[cls],
            'percent': round(pct, 6),
            'members_preview': members[cls][:20]  # first 20 matches
        }
    return result, members


def task_full(path, max_lines=None, classes=("class22",), our_classes=None):
    """
    Full analysis for small datasets:
    - membership check
    - algebraic degree and image_size for members
    """
    result, members = task_membership(path, max_lines, classes, verbose=False)

    # For each member, compute extra invariants
    for cls in classes:
        detailed = []
        if members[cls]:
            print(f"\nComputing invariants for {len(members[cls])} {cls} members...")
            for i, item in enumerate(members[cls][:100]):  # cap at 100
                line_no = item['line']
                # Re-read that line (small files only)
                sbox = None
                with open(path, 'r') as f:
                    for j, line in enumerate(f, 1):
                        if j == line_no:
                            sbox = parse_line(line)
                            break
                if sbox is None:
                    continue
                inv = {
                    'line': line_no,
                    'sbox_first8': sbox[:8],
                    'degree': algebraic_degree(sbox),
                    'image_size': image_size(sbox),
                    'perm': len(set(sbox)) == N,
                }
                detailed.append(inv)
                if (i+1) % 20 == 0:
                    print(f"  {i+1}/{min(len(members[cls]), 100)} done")
            result['classes'][cls]['detailed'] = detailed
    return result


# ---------------------------------------------------------------------------
# Load our known ortho-signatures for comparison
# ---------------------------------------------------------------------------

def load_our_signatures():
    """Load ortho-derivative signatures of our 6 CCZ classes."""
    path = VERIF_DIR / "found_apn_classes_report.json"
    if not path.exists():
        return {}
    with open(path) as f:
        d = json.load(f)
    sigs = {}
    for cls in d.get('classes', []):
        cls_id = cls['class_id']
        sig = cls.get('ortho_sig_hash', '')
        sigs[cls_id] = sig
    return sigs


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Analyze APN dataset against our search subspace.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Tasks:
  membership  Fast: count how many functions lie in each self-equiv subspace
  full        Membership + invariants (for small files, <50K functions)

Examples:
  python src/analyze_dataset.py --file data/new_apns.txt --task membership
  python src/analyze_dataset.py --file data/new_apns.txt --task membership --max-lines 50000
  python src/analyze_dataset.py --file data/apn_8bit.txt --task full
  python src/analyze_dataset.py --file data/new_apns.txt --task membership --classes class22 class30
""")
    parser.add_argument('--file', required=True, type=Path,
                        help='Input file (one sbox per line as [v0,...,v255])')
    parser.add_argument('--task', choices=['membership', 'full'], default='membership',
                        help='Analysis task (default: membership)')
    parser.add_argument('--classes', nargs='+', default=['class22'],
                        choices=list(CLASS_LOOKUPS.keys()),
                        help='Which self-equivalence subspaces to check (default: class22)')
    parser.add_argument('--max-lines', type=int, default=None,
                        help='Limit to first N lines (for testing)')
    parser.add_argument('--out-dir', type=Path, default=ANALYSIS_DIR,
                        help=f'Output directory (default: {ANALYSIS_DIR})')
    args = parser.parse_args()

    if not args.file.exists():
        print(f"ERROR: file not found: {args.file}")
        print(f"Put datasets in {DATA_DIR}/  (e.g. data/new_apns.txt)")
        sys.exit(1)

    ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    fname = args.file.stem
    suffix = f"_first{args.max_lines}" if args.max_lines else ""
    out_json = args.out_dir / f"{fname}{suffix}_{args.task}.json"
    out_txt  = args.out_dir / f"{fname}{suffix}_{args.task}.txt"

    print(f"File    : {args.file}")
    print(f"Task    : {args.task}")
    print(f"Classes : {args.classes}")
    if args.max_lines:
        print(f"Max lines: {args.max_lines:,}")
    print(f"Output  : {out_json}")
    print()

    t0 = time.time()

    if args.task == 'membership':
        result, _ = task_membership(args.file, args.max_lines, args.classes, verbose=True)
    else:
        result = task_full(args.file, args.max_lines, args.classes)

    elapsed = time.time() - t0

    # Save JSON
    with open(out_json, 'w', encoding='utf-8') as f:
        json.dump(result, f, indent=2)

    # Save text summary
    our_sigs = load_our_signatures()
    lines = []
    lines.append("=" * 65)
    lines.append(f"DATASET ANALYSIS: {args.file.name}")
    lines.append(f"Task: {args.task}  |  Total time: {elapsed:.1f}s")
    lines.append("=" * 65)
    lines.append(f"Total functions read : {result['total_lines']:,}")
    lines.append(f"Processing rate      : {result['rate_per_sec']:,.0f} functions/sec")
    lines.append("")
    for cls, cdata in result['classes'].items():
        lines.append(f"--- {cls.upper()} self-equivalence subspace ---")
        lines.append(f"  Members found : {cdata['count']:,}")
        lines.append(f"  Percentage    : {cdata['percent']:.6f}%")
        lines.append(f"  Expected dim  : {'40' if cls == 'class22' else '116'} (subspace size 2^dim)")
        if cdata.get('detailed'):
            lines.append(f"  Detailed invariants (first {len(cdata['detailed'])} members):")
            img_sizes = [d['image_size'] for d in cdata['detailed']]
            degs      = [d['degree']     for d in cdata['detailed']]
            lines.append(f"    degree range     : {min(degs)} - {max(degs)}")
            lines.append(f"    image_size range : {min(img_sizes)} - {max(img_sizes)}")
            lines.append(f"    unique image_sizes: {sorted(set(img_sizes))}")
        if cdata.get('members_preview'):
            lines.append(f"  First matches (line, sbox[0:8]):")
            for m in cdata['members_preview'][:5]:
                lines.append(f"    line {m['line']:>7}: {m['sbox_first8']}")
    lines.append("")
    if our_sigs:
        lines.append("Our CCZ classes (for reference):")
        for cls_id, sig in sorted(our_sigs.items()):
            lines.append(f"  {cls_id}: ortho_sig={sig}")
    lines.append("=" * 65)

    summary = "\n".join(lines)
    with open(out_txt, 'w', encoding='utf-8') as f:
        f.write(summary + "\n")

    print()
    print(summary)
    print(f"\nSaved: {out_json}")
    print(f"Saved: {out_txt}")


if __name__ == '__main__':
    main()
