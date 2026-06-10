/**********************************************************************
  invariant_test_found_sboxes.m

  Lightweight invariant test for the six found APN S-boxes.
  No CCZ code-equivalence call is used here.

  It computes:
    - permutation flag
    - image size
    - differential uniformity
    - DDT entry spectrum
    - vectorial Walsh absolute spectrum
    - vectorial autocorrelation absolute spectrum (optional)
    - algebraic degree of coordinate functions
    - exact equality / Hamming distances

  Usage:
      load "invariant_test_found_sboxes.m";
**********************************************************************/

n := 8;
N := 2^n;
MODPOLY := 285; // 0x11D
COMPUTE_AUTOCORR := true;  // set false if you want a faster run

// ------------------------------------------------------------
// GF(2^8) arithmetic, little-endian integer basis
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

function ParityInt(z)
    p := 0;
    zz := Integers()!z;
    while zz gt 0 do
        p := BitwiseXor(p, zz mod 2);
        zz := zz div 2;
    end while;
    return p;
end function;

function HammingWeightInt(z)
    w := 0;
    zz := Integers()!z;
    while zz gt 0 do
        w +:= zz mod 2;
        zz := zz div 2;
    end while;
    return w;
end function;

function IsPermutationSbox(S)
    return #Seqset(S) eq #S;
end function;

function ImageSize(S)
    return #Seqset(S);
end function;

function HammingDistanceSboxes(A, B)
    return #[ i : i in [1..#A] | A[i] ne B[i] ];
end function;

function SpectrumFromList(L)
    spec := AssociativeArray(Integers());
    for v in L do
        if IsDefined(spec, v) then
            spec[v] +:= 1;
        else
            spec[v] := 1;
        end if;
    end for;
    keys := Sort(Setseq(Keys(spec)));
    return [ <k, spec[k]> : k in keys ];
end function;

// ------------------------------------------------------------
// Cryptographic/vectorial invariants
// ------------------------------------------------------------

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

function DifferentialSpectrum(S)
    vals := [];
    for a in [1..N-1] do
        counts := [0 : i in [1..N]];
        for x in [0..N-1] do
            b := BitwiseXor(S[x+1], S[BitwiseXor(x,a)+1]);
            counts[b+1] +:= 1;
        end for;
        vals cat:= counts;
    end for;
    return SpectrumFromList(vals);
end function;

function WalshAbsSpectrum(S)
    vals := [];
    for u in [0..N-1] do
        for v in [1..N-1] do
            sum := 0;
            for x in [0..N-1] do
                parity := ParityInt(BitwiseXor(BitwiseAnd(u,x), BitwiseAnd(v,S[x+1])));
                if parity eq 0 then
                    sum +:= 1;
                else
                    sum -:= 1;
                end if;
            end for;
            Append(~vals, Abs(sum));
        end for;
    end for;
    return SpectrumFromList(vals);
end function;

function AutocorrelationAbsSpectrum(S)
    vals := [];
    for a in [1..N-1] do
        for v in [1..N-1] do
            sum := 0;
            for x in [0..N-1] do
                fx := S[x+1];
                fy := S[BitwiseXor(x,a)+1];
                parity := ParityInt(BitwiseAnd(v, BitwiseXor(fx, fy)));
                if parity eq 0 then
                    sum +:= 1;
                else
                    sum -:= 1;
                end if;
            end for;
            Append(~vals, Abs(sum));
        end for;
    end for;
    return SpectrumFromList(vals);
end function;

function ComponentAlgebraicDegrees(S)
    degs := [];
    for bit in [0..n-1] do
        f := [ (S[x+1] div 2^bit) mod 2 : x in [0..N-1] ];
        // Mobius transform truth table -> ANF coefficients.
        for i in [0..n-1] do
            step := 2^i;
            for mask in [0..N-1] do
                if BitwiseAnd(mask, step) ne 0 then
                    f[mask+1] := (f[mask+1] + f[mask-step+1]) mod 2;
                end if;
            end for;
        end for;
        d := -1;
        for mask in [0..N-1] do
            if f[mask+1] ne 0 then
                w := HammingWeightInt(mask);
                if w gt d then d := w; end if;
            end if;
        end for;
        Append(~degs, d);
    end for;
    return degs;
end function;

function AlgebraicDegree(S)
    return Maximum(ComponentAlgebraicDegrees(S));
end function;

function FullInvariantRecord(S)
    rec := AssociativeArray();
    rec["perm"] := IsPermutationSbox(S);
    rec["image"] := ImageSize(S);
    rec["DU"] := DifferentialUniformity(S);
    rec["DDT"] := DifferentialSpectrum(S);
    rec["WalshAbs"] := WalshAbsSpectrum(S);
    rec["coord_degs"] := ComponentAlgebraicDegrees(S);
    rec["alg_deg"] := Maximum(rec["coord_degs"]);
    if COMPUTE_AUTOCORR then
        rec["ACAbs"] := AutocorrelationAbsSpectrum(S);
    end if;
    return rec;
end function;

procedure PrintRecord(label, rec, gold_rec)
    print "-------------------------------------------------------";
    print label;
    print "  permutation:", rec["perm"], " same as Gold:", rec["perm"] eq gold_rec["perm"];
    print "  image_size :", rec["image"], " same as Gold:", rec["image"] eq gold_rec["image"];
    print "  DU         :", rec["DU"], " same as Gold:", rec["DU"] eq gold_rec["DU"];
    print "  alg_degree :", rec["alg_deg"], " same as Gold:", rec["alg_deg"] eq gold_rec["alg_deg"];
    print "  coord_degs :", rec["coord_degs"], " same as Gold:", rec["coord_degs"] eq gold_rec["coord_degs"];
    print "  DDT spectrum same as Gold:", rec["DDT"] eq gold_rec["DDT"];
    print "    DDT:", rec["DDT"];
    print "  Walsh abs spectrum same as Gold:", rec["WalshAbs"] eq gold_rec["WalshAbs"];
    print "    WalshAbs:", rec["WalshAbs"];
    if COMPUTE_AUTOCORR then
        print "  Autocorr abs spectrum same as Gold:", rec["ACAbs"] eq gold_rec["ACAbs"];
        print "    ACAbs:", rec["ACAbs"];
    end if;
end procedure;

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
print "=======================================================";
print "LIGHTWEIGHT INVARIANT TEST: FOUND APN S-BOXES";
print "Reference: Gold x^3, modulus 285 = 0x11D";
print "Autocorrelation enabled:", COMPUTE_AUTOCORR;
print "Number of found S-boxes:", #FOUND_SBOXES;
print "=======================================================";

t0 := Cputime();
print "Computing Gold invariants...";
gold_rec := FullInvariantRecord(GOLD_SBOX);
PrintRecord("Gold x^3", gold_rec, gold_rec);
print "Gold invariants time:", Cputime(t0), "sec";

records := [];
for item in FOUND_SBOXES do
    label := item[1];
    S := item[2];
    t := Cputime();
    rec := FullInvariantRecord(S);
    Append(~records, <label, rec>);
    PrintRecord(label, rec, gold_rec);
    print "  time:", Cputime(t), "sec";
end for;

print "=======================================================";
print "PAIRWISE EXACT COMPARISON";
print "=======================================================";
for i in [1..#FOUND_SBOXES] do
    for j in [i+1..#FOUND_SBOXES] do
        A := FOUND_SBOXES[i][2];
        B := FOUND_SBOXES[j][2];
        print FOUND_SBOXES[i][1], "vs", FOUND_SBOXES[j][1],
              " exact_equal=", A eq B,
              " hamming_distance=", HammingDistanceSboxes(A,B);
    end for;
end for;

print "=======================================================";
print "PAIRWISE INVARIANT EQUALITY";
print "=======================================================";
for i in [1..#records] do
    for j in [i+1..#records] do
        ri := records[i][2];
        rj := records[j][2];
        same := (ri["perm"] eq rj["perm"]) and
                (ri["image"] eq rj["image"]) and
                (ri["DU"] eq rj["DU"]) and
                (ri["DDT"] eq rj["DDT"]) and
                (ri["WalshAbs"] eq rj["WalshAbs"]) and
                (ri["alg_deg"] eq rj["alg_deg"]) and
                (ri["coord_degs"] eq rj["coord_degs"]);
        if COMPUTE_AUTOCORR then
            same := same and (ri["ACAbs"] eq rj["ACAbs"]);
        end if;
        print records[i][1], "vs", records[j][1], " all_invariants_equal=", same;
    end for;
end for;

print "=======================================================";
print "DONE. Total time:", Cputime(t0), "sec";
print "Note: matching invariants do NOT prove CCZ-equivalence; they are filters.";
print "=======================================================";
