/**********************************************************************
  ccz_test_found_sboxes_full.m

  Complete CCZ-equivalence test for the APN S-boxes found in:
    - seed03: 4 solutions
    - seed06: 2 solutions

  The reference is Gold x^3 over GF(2^8) with polynomial modulus 0x11D
  (decimal 285), using the same little-endian integer/bit convention as
  the APN scripts.

  Method: graph-code criterion.
  For F : F_2^8 -> F_2^8, build the binary linear code with parity-check
  matrix columns (1, x, F(x)). CCZ-equivalence corresponds to equivalence
  of these graph codes.

  Usage:
      load "ccz_test_found_sboxes_full.m";
**********************************************************************/

SetVerbose("Code", 0);

n := 8;
N := 2^n;
MODPOLY := 285; // x^8 + x^4 + x^3 + x^2 + 1 = 0x11D

// ------------------------------------------------------------
// GF(2^8) arithmetic in integer polynomial basis, little-endian
// ------------------------------------------------------------

function GF2MulInt(a, b, modpoly)
    aa := Integers()!a;
    bb := Integers()!b;
    res := 0;

    while bb gt 0 do
        if (bb mod 2) eq 1 then
            res := BitwiseXor(res, aa);
        end if;
        bb := bb div 2;
        aa := 2*aa;
        if aa ge 256 then
            aa := BitwiseXor(aa, modpoly);
        end if;
    end while;

    return res mod 256;
end function;

function GF2PowInt(a, e, modpoly)
    result := 1;
    base := Integers()!a;
    ee := Integers()!e;

    while ee gt 0 do
        if (ee mod 2) eq 1 then
            result := GF2MulInt(result, base, modpoly);
        end if;
        base := GF2MulInt(base, base, modpoly);
        ee := ee div 2;
    end while;

    return result;
end function;

function GoldSboxPower3()
    return [ GF2PowInt(x, 3, MODPOLY) : x in [0..255] ];
end function;

// ------------------------------------------------------------
// Basic utilities
// ------------------------------------------------------------

function IntToBitsLE(a, n)
    bits := [];
    x := Integers()!a;
    for i in [1..n] do
        Append(~bits, x mod 2);
        x := x div 2;
    end for;
    return bits;
end function;

function CheckSbox(S)
    if #S ne N then
        print "Bad S-box length:", #S, "expected", N;
        return false;
    end if;
    for v in S do
        if (v lt 0) or (v ge N) then
            print "Bad S-box value:", v;
            return false;
        end if;
    end for;
    return true;
end function;

function IsPermutationSbox(S)
    return #Seqset(S) eq #S;
end function;

function DifferentialUniformity(S)
    maxdu := 0;
    for a in [1..N-1] do
        counts := [0 : i in [1..N]];
        for x in [0..N-1] do
            b := BitwiseXor(S[x+1], S[BitwiseXor(x,a)+1]);
            counts[b+1] +:= 1;
        end for;
        ma := Maximum(counts);
        if ma gt maxdu then maxdu := ma; end if;
    end for;
    return maxdu;
end function;

function WalshAbsSpectrum(S)
    // Returns sorted pairs <abs_value, multiplicity> for vectorial Walsh transform.
    spec := AssociativeArray(Integers());
    for u in [0..N-1] do
        for v in [1..N-1] do
            sum := 0;
            for x in [0..N-1] do
                ux := BitwiseAnd(u, x);
                vy := BitwiseAnd(v, S[x+1]);
                parity := 0;
                z := BitwiseXor(ux, vy);
                while z gt 0 do
                    parity := BitwiseXor(parity, z mod 2);
                    z := z div 2;
                end while;
                if parity eq 0 then
                    sum +:= 1;
                else
                    sum -:= 1;
                end if;
            end for;
            av := Abs(sum);
            if IsDefined(spec, av) then
                spec[av] +:= 1;
            else
                spec[av] := 1;
            end if;
        end for;
    end for;
    keys := Sort(Setseq(Keys(spec)));
    return [ <k, spec[k]> : k in keys ];
end function;

function GraphParityCheckMatrix(S)
    F2 := GF(2);
    H := ZeroMatrix(F2, 1 + 2*n, N);

    for x in [0..N-1] do
        col := x + 1;
        H[1, col] := F2!1;

        xb := IntToBitsLE(x, n);
        yb := IntToBitsLE(S[col], n);

        for i in [1..n] do
            H[1+i, col] := F2!xb[i];
            H[1+n+i, col] := F2!yb[i];
        end for;
    end for;

    return H;
end function;

function GraphKernelCode(S)
    H := GraphParityCheckMatrix(S);
    return Dual(LinearCode(H));
end function;

function CCZEquivalentByCode(S, T)
    if not CheckSbox(S) then return false; end if;
    if not CheckSbox(T) then return false; end if;

    C1 := GraphKernelCode(S);
    C2 := GraphKernelCode(T);

    print "  Code parameters C1:", Length(C1), Dimension(C1);
    print "  Code parameters C2:", Length(C2), Dimension(C2);

    ok := false;
    try
        ok := IsIsomorphic(C1, C2);
    catch e
        try
            ok, iso := IsIsomorphic(C1, C2);
        catch ee
            print "  ERROR: IsIsomorphic failed.";
            print ee`Object;
            return false;
        end try;
    end try;

    return ok;
end function;

procedure TestOne(label, S, GOLD)
    print "=======================================================";
    print "Testing", label;
    print "  length:", #S;
    print "  permutation:", IsPermutationSbox(S);
    print "  DU:", DifferentialUniformity(S);
    print "  Walsh_abs_spectrum:", WalshAbsSpectrum(S);
    ok := CCZEquivalentByCode(S, GOLD);
    print "  CCZ-equivalent to Gold x^3:", ok;
end procedure;

// ------------------------------------------------------------
// Reference Gold x^3 and found S-boxes
// ------------------------------------------------------------

GOLD_SBOX := GoldSboxPower3();

FOUND_SBOXES := [
    <"seed03_sol1", [
        0, 0, 0, 63, 0, 42, 97, 116, 0, 253, 75, 137, 194, 21, 232, 0,
        0, 24, 207, 232, 165, 151, 11, 6, 183, 82, 51, 233, 208, 31, 53, 197,
        0, 122, 47, 106, 113, 33, 63, 80, 186, 61, 222, 102, 9, 164, 12, 158,
        19, 113, 243, 174, 199, 143, 70, 49, 30, 129, 181, 21, 8, 189, 194, 72,
        0, 211, 235, 7, 65, 184, 203, 13, 226, 204, 66, 83, 97, 101, 160, 155,
        178, 121, 150, 98, 86, 183, 19, 205, 231, 209, 136, 129, 193, 221, 207, 236,
        144, 57, 84, 194, 160, 35, 5, 185, 200, 156, 71, 44, 58, 68, 212, 149,
        49, 128, 58, 180, 164, 63, 206, 106, 222, 146, 158, 237, 137, 239, 168, 241,
        0, 217, 86, 176, 57, 202, 14, 194, 130, 166, 159, 132, 121, 119, 5, 52,
        95, 158, 198, 56, 195, 40, 59, 239, 106, 86, 184, 187, 52, 34, 135, 174,
        34, 129, 91, 199, 106, 227, 114, 196, 26, 68, 40, 73, 144, 228, 195, 136,
        110, 213, 216, 92, 131, 18, 84, 250, 225, 167, 28, 101, 206, 162, 82, 1,
        76, 70, 241, 196, 52, 20, 232, 247, 44, 219, 218, 18, 150, 75, 1, 227,
        161, 179, 211, 254, 124, 68, 111, 104, 118, 153, 79, 159, 105, 172, 49, 203,
        254, 142, 108, 35, 247, 173, 4, 97, 36, 169, 253, 79, 239, 72, 87, 207,
        0, 104, 93, 10, 172, 238, 144, 237, 109, 248, 123, 209, 3, 188, 116, 244
    ]>,
    <"seed03_sol2", [
        0, 0, 0, 63, 0, 42, 97, 116, 0, 253, 75, 137, 194, 21, 232, 0,
        0, 24, 207, 232, 165, 151, 11, 6, 183, 82, 51, 233, 208, 31, 53, 197,
        0, 122, 47, 106, 113, 33, 63, 80, 186, 61, 222, 102, 9, 164, 12, 158,
        19, 113, 243, 174, 199, 143, 70, 49, 30, 129, 181, 21, 8, 189, 194, 72,
        0, 211, 235, 7, 65, 184, 203, 13, 226, 204, 66, 83, 97, 101, 160, 155,
        178, 121, 150, 98, 86, 183, 19, 205, 231, 209, 136, 129, 193, 221, 207, 236,
        38, 143, 226, 116, 22, 149, 179, 15, 126, 42, 241, 154, 140, 242, 98, 35,
        135, 54, 140, 2, 18, 137, 120, 220, 104, 36, 40, 91, 63, 89, 30, 71,
        0, 217, 86, 176, 57, 202, 14, 194, 130, 166, 159, 132, 121, 119, 5, 52,
        95, 158, 198, 56, 195, 40, 59, 239, 106, 86, 184, 187, 52, 34, 135, 174,
        148, 55, 237, 113, 220, 85, 196, 114, 172, 242, 158, 255, 38, 82, 117, 62,
        216, 99, 110, 234, 53, 164, 226, 76, 87, 17, 170, 211, 120, 20, 228, 183,
        76, 70, 241, 196, 52, 20, 232, 247, 44, 219, 218, 18, 150, 75, 1, 227,
        161, 179, 211, 254, 124, 68, 111, 104, 118, 153, 79, 159, 105, 172, 49, 203,
        254, 142, 108, 35, 247, 173, 4, 97, 36, 169, 253, 79, 239, 72, 87, 207,
        0, 104, 93, 10, 172, 238, 144, 237, 109, 248, 123, 209, 3, 188, 116, 244
    ]>,
    <"seed03_sol3", [
        0, 0, 0, 63, 0, 42, 97, 116, 0, 253, 75, 137, 194, 21, 232, 0,
        0, 24, 207, 232, 165, 151, 11, 6, 183, 82, 51, 233, 208, 31, 53, 197,
        0, 122, 47, 106, 113, 33, 63, 80, 186, 61, 222, 102, 9, 164, 12, 158,
        19, 113, 243, 174, 199, 143, 70, 49, 30, 129, 181, 21, 8, 189, 194, 72,
        0, 211, 235, 7, 65, 184, 203, 13, 226, 204, 66, 83, 97, 101, 160, 155,
        178, 121, 150, 98, 86, 183, 19, 205, 231, 209, 136, 129, 193, 221, 207, 236,
        25, 176, 221, 75, 41, 170, 140, 48, 65, 21, 206, 165, 179, 205, 93, 28,
        184, 9, 179, 61, 45, 182, 71, 227, 87, 27, 23, 100, 0, 102, 33, 120,
        0, 217, 86, 176, 57, 202, 14, 194, 130, 166, 159, 132, 121, 119, 5, 52,
        95, 158, 198, 56, 195, 40, 59, 239, 106, 86, 184, 187, 52, 34, 135, 174,
        171, 8, 210, 78, 227, 106, 251, 77, 147, 205, 161, 192, 25, 109, 74, 1,
        231, 92, 81, 213, 10, 155, 221, 115, 104, 46, 149, 236, 71, 43, 219, 136,
        76, 70, 241, 196, 52, 20, 232, 247, 44, 219, 218, 18, 150, 75, 1, 227,
        161, 179, 211, 254, 124, 68, 111, 104, 118, 153, 79, 159, 105, 172, 49, 203,
        254, 142, 108, 35, 247, 173, 4, 97, 36, 169, 253, 79, 239, 72, 87, 207,
        0, 104, 93, 10, 172, 238, 144, 237, 109, 248, 123, 209, 3, 188, 116, 244
    ]>,
    <"seed03_sol4", [
        0, 0, 0, 63, 0, 42, 97, 116, 0, 253, 75, 137, 194, 21, 232, 0,
        0, 24, 207, 232, 165, 151, 11, 6, 183, 82, 51, 233, 208, 31, 53, 197,
        0, 122, 47, 106, 113, 33, 63, 80, 186, 61, 222, 102, 9, 164, 12, 158,
        19, 113, 243, 174, 199, 143, 70, 49, 30, 129, 181, 21, 8, 189, 194, 72,
        0, 211, 235, 7, 65, 184, 203, 13, 226, 204, 66, 83, 97, 101, 160, 155,
        178, 121, 150, 98, 86, 183, 19, 205, 231, 209, 136, 129, 193, 221, 207, 236,
        175, 6, 107, 253, 159, 28, 58, 134, 247, 163, 120, 19, 5, 123, 235, 170,
        14, 191, 5, 139, 155, 0, 241, 85, 225, 173, 161, 210, 182, 208, 151, 206,
        0, 217, 86, 176, 57, 202, 14, 194, 130, 166, 159, 132, 121, 119, 5, 52,
        95, 158, 198, 56, 195, 40, 59, 239, 106, 86, 184, 187, 52, 34, 135, 174,
        29, 190, 100, 248, 85, 220, 77, 251, 37, 123, 23, 118, 175, 219, 252, 183,
        81, 234, 231, 99, 188, 45, 107, 197, 222, 152, 35, 90, 241, 157, 109, 62,
        76, 70, 241, 196, 52, 20, 232, 247, 44, 219, 218, 18, 150, 75, 1, 227,
        161, 179, 211, 254, 124, 68, 111, 104, 118, 153, 79, 159, 105, 172, 49, 203,
        254, 142, 108, 35, 247, 173, 4, 97, 36, 169, 253, 79, 239, 72, 87, 207,
        0, 104, 93, 10, 172, 238, 144, 237, 109, 248, 123, 209, 3, 188, 116, 244
    ]>,
    <"seed06_sol1", [
        0, 0, 0, 247, 0, 175, 30, 70, 0, 212, 177, 146, 35, 88, 140, 0,
        0, 105, 168, 54, 8, 206, 190, 143, 186, 7, 163, 233, 145, 131, 150, 115,
        0, 230, 205, 220, 191, 246, 108, 210, 15, 61, 115, 182, 147, 14, 241, 155,
        57, 182, 92, 36, 142, 174, 245, 34, 140, 215, 88, 244, 24, 236, 210, 209,
        0, 217, 60, 18, 117, 3, 87, 214, 145, 156, 28, 230, 199, 101, 84, 1,
        46, 158, 186, 253, 83, 76, 217, 49, 5, 97, 32, 179, 91, 144, 96, 92,
        8, 55, 249, 49, 194, 82, 45, 74, 150, 125, 214, 202, 127, 59, 33, 146,
        31, 73, 70, 231, 221, 36, 154, 148, 59, 185, 211, 166, 218, 247, 44, 246,
        0, 2, 93, 168, 103, 202, 36, 126, 234, 60, 6, 39, 174, 215, 92, 210,
        252, 151, 9, 149, 147, 87, 120, 75, 172, 19, 232, 160, 224, 240, 186, 93,
        38, 194, 182, 165, 254, 181, 112, 204, 195, 243, 226, 37, 56, 167, 7, 111,
        227, 110, 219, 161, 51, 17, 21, 192, 188, 229, 53, 155, 79, 185, 216, 217,
        197, 30, 164, 136, 215, 163, 168, 43, 190, 177, 110, 150, 143, 47, 65, 22,
        23, 165, 222, 155, 13, 16, 218, 48, 214, 176, 174, 63, 239, 38, 137, 183,
        235, 214, 71, 141, 70, 212, 244, 145, 159, 118, 130, 156, 17, 87, 18, 163,
        0, 84, 4, 167, 165, 94, 191, 179, 206, 78, 123, 12, 72, 103, 227, 59
    ]>,
    <"seed06_sol2", [
        0, 0, 0, 247, 0, 175, 30, 70, 0, 212, 177, 146, 35, 88, 140, 0,
        0, 105, 168, 54, 8, 206, 190, 143, 186, 7, 163, 233, 145, 131, 150, 115,
        0, 230, 205, 220, 191, 246, 108, 210, 15, 61, 115, 182, 147, 14, 241, 155,
        57, 182, 92, 36, 142, 174, 245, 34, 140, 215, 88, 244, 24, 236, 210, 209,
        0, 217, 60, 18, 117, 3, 87, 214, 145, 156, 28, 230, 199, 101, 84, 1,
        46, 158, 186, 253, 83, 76, 217, 49, 5, 97, 32, 179, 91, 144, 96, 92,
        109, 82, 156, 84, 167, 55, 72, 47, 243, 24, 179, 175, 26, 94, 68, 247,
        122, 44, 35, 130, 184, 65, 255, 241, 94, 220, 182, 195, 191, 146, 73, 147,
        0, 2, 93, 168, 103, 202, 36, 126, 234, 60, 6, 39, 174, 215, 92, 210,
        252, 151, 9, 149, 147, 87, 120, 75, 172, 19, 232, 160, 224, 240, 186, 93,
        67, 167, 211, 192, 155, 208, 21, 169, 166, 150, 135, 64, 93, 194, 98, 10,
        134, 11, 190, 196, 86, 116, 112, 165, 217, 128, 80, 254, 42, 220, 189, 188,
        197, 30, 164, 136, 215, 163, 168, 43, 190, 177, 110, 150, 143, 47, 65, 22,
        23, 165, 222, 155, 13, 16, 218, 48, 214, 176, 174, 63, 239, 38, 137, 183,
        235, 214, 71, 141, 70, 212, 244, 145, 159, 118, 130, 156, 17, 87, 18, 163,
        0, 84, 4, 167, 165, 94, 191, 179, 206, 78, 123, 12, 72, 103, 227, 59
    ]>
];



// ------------------------------------------------------------
// Optimized CCZ test via the DUAL graph code
//
// Original graph-code criterion often compares C_F = Nullspace(H_F),
// a [256,239] code.  This is heavy.
// Since permutation equivalence preserves duals, C_F ~ C_G iff
// C_F^perp ~ C_G^perp.  Here C_F^perp is just the row space of H_F,
// a [256,17] code.  This is usually much faster.
// ------------------------------------------------------------

function GraphDualCode(S)
    H := GraphParityCheckMatrix(S);
    return LinearCode(H);
end function;

function CodeWeightDistributionCompact(C)
    WD := WeightDistribution(C);

    // In Magma V2.20 WeightDistribution(C) usually returns
    // a sequence of tuples <weight, multiplicity>.
    // Keep it in this compact canonical form.
    if #WD eq 0 then
        return [];
    end if;

    return [ <Integers()!t[1], Integers()!t[2]> : t in WD | Integers()!t[2] ne 0 ];
end function;

function TryIsIsomorphicSmall(C1, C2)
    ok := false;
    iso_ok := false;

    try
        ok := IsIsomorphic(C1, C2);
        iso_ok := true;
    catch e
        try
            ok, iso := IsIsomorphic(C1, C2);
            iso_ok := true;
        catch ee
            print "  ERROR in IsIsomorphic:";
            print ee`Object;
            return false, false;
        end try;
    end try;

    return ok, iso_ok;
end function;

procedure PrintDualCodeBasic(label, C)
    print label;
    print "  length:", Length(C), "dimension:", Dimension(C);
    print "  min distance:", MinimumDistance(C);
    print "  weight distribution:", CodeWeightDistributionCompact(C);
end procedure;

function FastCCZEquivalentByDualCode(S, T)
    if not CheckSbox(S) then return false; end if;
    if not CheckSbox(T) then return false; end if;

    C1 := GraphDualCode(S);
    C2 := GraphDualCode(T);

    print "  Dual graph code parameters:";
    print "    C1:", Length(C1), Dimension(C1), "d=", MinimumDistance(C1);
    print "    C2:", Length(C2), Dimension(C2), "d=", MinimumDistance(C2);

    WD1 := CodeWeightDistributionCompact(C1);
    WD2 := CodeWeightDistributionCompact(C2);

    if WD1 ne WD2 then
        print "  Dual weight distributions differ -> NOT CCZ-equivalent.";
        return false;
    end if;

    print "  Dual weight distributions match.";
    print "  Running IsIsomorphic on [256,17] dual graph codes...";
    t := Cputime();
    ok, iso_ok := TryIsIsomorphicSmall(C1, C2);
    print "  IsIsomorphic done in", Cputime(t), "sec";
    if not iso_ok then
        print "  IsIsomorphic failed.";
        return false;
    end if;

    return ok;
end function;

procedure TestOneFast(label, S, GOLD)
    print "=======================================================";
    print "Testing", label, "vs Gold x^3";
    print "  length:", #S;
    print "  permutation:", IsPermutationSbox(S);
    print "  DU:", DifferentialUniformity(S);
    ok := FastCCZEquivalentByDualCode(S, GOLD);
    print "  FAST CCZ-equivalent to Gold x^3:", ok;
end procedure;

print "=======================================================";
print "FAST CCZ TEST USING DUAL GRAPH CODES";
print "Reference: Gold x^3, modulus 285 = 0x11D";
print "Gold permutation:", IsPermutationSbox(GOLD_SBOX);
print "Gold DU:", DifferentialUniformity(GOLD_SBOX);
print "Number of found S-boxes:", #FOUND_SBOXES;
print "=======================================================";

CG := GraphDualCode(GOLD_SBOX);
PrintDualCodeBasic("Gold dual graph code:", CG);

print "=======================================================";
print "TEST FOUND S-BOXES AGAINST GOLD";
print "=======================================================";

for rec in FOUND_SBOXES do
    TestOneFast(rec[1], rec[2], GOLD_SBOX);
end for;

print "=======================================================";
print "PAIRWISE FAST CCZ TESTS AMONG FOUND S-BOXES";
print "=======================================================";

for i in [1..#FOUND_SBOXES] do
    for j in [i+1..#FOUND_SBOXES] do
        print "-------------------------------------------------------";
        print "Pair", FOUND_SBOXES[i][1], "vs", FOUND_SBOXES[j][1];
        tpair := Cputime();
        ok := FastCCZEquivalentByDualCode(FOUND_SBOXES[i][2], FOUND_SBOXES[j][2]);
        print "  FAST CCZ-equivalent:", ok, "time:", Cputime(tpair), "sec";
    end for;
end for;

print "DONE.";
