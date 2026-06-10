#!/usr/bin/env python3
"""
make_gold_test.py
=================
Creates gold_test.json with canonical Gold x^3 function
for testing check_vs_beierle_db.py.

Run in WSL:
  conda activate sage
  python make_gold_test.py
  python src/check_vs_beierle_db.py --found gold_test.json --known data/new_apns.txt --max-known 100000 --progress 10000

If Gold x^3 gives HIT -> pipeline correct, our new functions are genuinely new.
If Gold x^3 gives NO_MATCH -> problem with database format or pipeline.
"""

import json, hashlib

N, n = 256, 8
POLY = 285  # x^8 + x^4 + x^3 + x^2 + 1

def gf_mul(a, b):
    r = 0
    while b:
        if b & 1: r ^= a
        a <<= 1
        if a & 256: a ^= POLY
        b >>= 1
    return r

def gf_pow(a, k):
    r = 1
    for _ in range(k): r = gf_mul(r, a)
    return r

# Gold x^3: F(x) = x^3 in GF(2^8)
gold3  = [gf_pow(x, 3)  for x in range(N)]
gold9  = [gf_pow(x, 9)  for x in range(N)]
gold33 = [gf_pow(x, 33) for x in range(N)]

print(f"Gold x^3:  F(0..7) = {gold3[:8]}")
print(f"Gold x^9:  F(0..7) = {gold9[:8]}")
print(f"Gold x^33: F(0..7) = {gold33[:8]}")
print(f"All start with 0 (F(0)=0): {gold3[0]==0 and gold9[0]==0 and gold33[0]==0}")

def sh(sb): return hashlib.sha256(bytes(sb)).hexdigest()[:16]

# Check APN
def is_apn(sb):
    from collections import Counter
    for a in range(1, N):
        if max(Counter(sb[x]^sb[x^a] for x in range(N)).values()) > 2:
            return False
    return True

print(f"Gold x^3  APN: {is_apn(gold3)}")
print(f"Gold x^9  APN: {is_apn(gold9)}")

# Write test JSON   format compatible with check_vs_beierle_db.py
data = {
    "source": "gold_test",
    "functions": {
        "GOLD-x3": {
            "apn_id": "GOLD-x3",
            "sbox": gold3,
            "sbox_hash": sh(gold3),
            "perm": False,
            "ccz_status": "known_gold",
        },
        "GOLD-x9": {
            "apn_id": "GOLD-x9",
            "sbox": gold9,
            "sbox_hash": sh(gold9),
            "perm": False,
            "ccz_status": "known_gold",
        },
        "GOLD-x33": {
            "apn_id": "GOLD-x33",
            "sbox": gold33,
            "sbox_hash": sh(gold33),
            "perm": False,
            "ccz_status": "known_gold",
        },
    }
}

with open("gold_test.json", "w") as f:
    json.dump(data, f, indent=2)

print()
print("Written: gold_test.json")
print()
print("Now run:")
print("  python src/check_vs_beierle_db.py \\")
print("    --found gold_test.json \\")
print("    --known data/new_apns.txt \\")
print("    --max-known 100000 \\")
print("    --progress 10000 \\")
print("    --skip-pairwise")
print()
print("Expected result:")
print("  GOLD-x3 => HIT (Gold x^3 must be in Beierle 2025)")
print("  GOLD-x9 => HIT (Gold x^9 must be in Beierle 2025)")
print()
print("If NO_MATCH -> check database format:")
print("  head -3 data/new_apns.txt")
