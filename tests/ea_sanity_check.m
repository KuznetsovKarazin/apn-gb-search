// =======================================================
// EA sanity check for FAST dual graph-code isomorphism test
// Magma V2.20-safe
//
// Builds Gold F(x)=x^3 over GF(2^8) with modulus 0x11D = 285,
// constructs several deliberately EA-equivalent functions, and checks
// whether the fast dual graph-code test returns true.
//
// Run:
//   load "ea_sanity_check_dual_fast.m";
// =======================================================

SetVerbose("Code", 0);

n := 8;
 N := 2^n;
MODPOLY := 285; // 0x11D = x^8+x^4+x^3+x^2+1

// -------------------------------
// GF(2^8) arithmetic on integers
// -------------------------------
function GFAdd(a,b)
    return BitwiseXor(a,b);
end function;

function GFMul(a,b)
    res := 0;
    aa := a;
    bb := b;
    while bb gt 0 do
        if (bb mod 2) eq 1 then
            res := BitwiseXor(res, aa);
        end if;
        bb := bb div 2;
        aa := 2*aa;
        if aa ge 256 then
            aa := BitwiseXor(aa, MODPOLY);
        end if;
    end while;
    return res mod 256;
end function;

function GFPow(a,e)
    r := 1;
    x := a;
    ee := e;
    while ee gt 0 do
        if (ee mod 2) eq 1 then
            r := GFMul(r,x);
        end if;
        x := GFMul(x,x);
        ee := ee div 2;
    end while;
    return r;
end function;

function GoldSbox()
    return [ GFPow(x,3) : x in [0..255] ];
end function;

// -------------------------------
// bit/vector conversions
// -------------------------------
function IntToBits(x, m)
    return [ ((x div 2^i) mod 2) : i in [0..m-1] ];
end function;

function BitsToInt(bits)
    s := 0;
    for i in [1..#bits] do
        if bits[i] mod 2 eq 1 then
            s +:= 2^(i-1);
        end if;
    end for;
    return s;
end function;

function Parity(x)
    p := 0;
    y := x;
    while y gt 0 do
        p := (p + (y mod 2)) mod 2;
        y := y div 2;
    end while;
    return p;
end function;

// -------------------------------
// simple invertible linear maps on F_2^8
// -------------------------------
function FrobeniusMapInt(x, k)
    // x -> x^(2^k), F_2-linear over GF(2^8)
    y := x;
    for i in [1..k] do
        y := GFMul(y,y);
    end for;
    return y;
end function;

function MulConstMapInt(x, c)
    // x -> c*x, invertible if c != 0
    return GFMul(c,x);
end function;

function RotateLeft8(x, r)
    bits := IntToBits(x,8);
    nb := [ bits[((i-r-1) mod 8)+1] : i in [1..8] ];
    return BitsToInt(nb);
end function;

function LinearMaskMapInt(x, masks)
    // output bit j is parity(x & masks[j])
    bits := [];
    for j in [1..8] do
        Append(~bits, Parity(BitwiseAnd(x, masks[j])));
    end for;
    return BitsToInt(bits);
end function;

// -------------------------------
// EA-equivalent constructors
// -------------------------------
function AddLinearOutput(S, masks)
    // G(x) = S(x) + L(x)
    return [ BitwiseXor(S[x+1], LinearMaskMapInt(x, masks)) : x in [0..255] ];
end function;

function OutputLinear(S, c)
    // G(x) = c*S(x), c != 0
    return [ MulConstMapInt(S[x+1], c) : x in [0..255] ];
end function;

function InputLinear(S, c)
    // G(x) = S(c*x), c != 0
    return [ S[MulConstMapInt(x,c)+1] : x in [0..255] ];
end function;

function InputFrobenius(S, k)
    // G(x) = S(x^(2^k))
    return [ S[FrobeniusMapInt(x,k)+1] : x in [0..255] ];
end function;

function OutputFrobenius(S, k)
    // G(x) = S(x)^(2^k)
    return [ FrobeniusMapInt(S[x+1],k) : x in [0..255] ];
end function;

function AddAffineConstant(S, b)
    // G(x)=S(x)+b
    return [ BitwiseXor(S[x+1], b) : x in [0..255] ];
end function;

function FullEAExample(S)
    // G(x)= A(S(Bx + a)) + L(x) + b
    // Here B: multiplication by 5; input shift a=37;
    // A: Frobenius^3 followed by multiplication by 17;
    // L: simple invertible-ish linear mask map; b=91.
    masks := [1,2,4,8,16,32,64,128];
    G := [];
    for x in [0..255] do
        bx := BitwiseXor(MulConstMapInt(x,5), 37);
        y := S[bx+1];
        y := FrobeniusMapInt(y,3);
        y := MulConstMapInt(y,17);
        y := BitwiseXor(y, LinearMaskMapInt(x,masks));
        y := BitwiseXor(y, 91);
        Append(~G, y);
    end for;
    return G;
end function;

// -------------------------------
// basic cryptographic checks
// -------------------------------
function IsPermutation(S)
    return #Seqset(S) eq #S;
end function;

function DifferentialUniformity(S)
    maxv := 0;
    for a in [1..255] do
        counts := [0 : i in [0..255]];
        for x in [0..255] do
            b := BitwiseXor(S[x+1], S[BitwiseXor(x,a)+1]);
            counts[b+1] +:= 1;
        end for;
        m := Max(counts);
        if m gt maxv then maxv := m; end if;
    end for;
    return maxv;
end function;

function ImageSize(S)
    return #Seqset(S);
end function;

// -------------------------------
// dual graph-code construction
// H has rows: 1, input bits, output bits. C = row span(H), dim <= 17.
// Isomorphism of these [256,17] codes is used as the fast test.
// -------------------------------
function DualGraphCode(S)
    F2 := GF(2);
    rows := [];

    Append(~rows, [ F2!1 : x in [0..255] ]);

    for i in [0..7] do
        Append(~rows, [ F2!(((x div 2^i) mod 2)) : x in [0..255] ]);
    end for;

    for j in [0..7] do
        Append(~rows, [ F2!(((S[x+1] div 2^j) mod 2)) : x in [0..255] ]);
    end for;

    M := Matrix(F2, rows);
    return LinearCode(M);
end function;

function WDCompact(C)
    WD := WeightDistribution(C);
    return [ <Integers()!t[1], Integers()!t[2]> : t in WD | Integers()!t[2] ne 0 ];
end function;

function FastDualCodeEquivalent(S,T)
    C1 := DualGraphCode(S);
    C2 := DualGraphCode(T);

    if Length(C1) ne Length(C2) then return false, 0.0; end if;
    if Dimension(C1) ne Dimension(C2) then return false, 0.0; end if;
    if MinimumDistance(C1) ne MinimumDistance(C2) then return false, 0.0; end if;
    if WDCompact(C1) ne WDCompact(C2) then return false, 0.0; end if;

    t := Cputime();
    ok := IsIsomorphic(C1, C2);
    return ok, Cputime(t);
end function;

procedure PrintSboxBasic(label, S, Gold)
    print "-------------------------------------------------------";
    print label;
    print "  permutation:", IsPermutation(S), " same as Gold:", IsPermutation(S) eq IsPermutation(Gold);
    print "  image size :", ImageSize(S), " same as Gold:", ImageSize(S) eq ImageSize(Gold);
    print "  DU         :", DifferentialUniformity(S), " same as Gold:", DifferentialUniformity(S) eq DifferentialUniformity(Gold);
end procedure;

// =======================================================
// Main
// =======================================================
Gold := GoldSbox();

// masks for L(x). Identity masks: L(x)=x. Nontrivial but invertible.
IDMasks := [1,2,4,8,16,32,64,128];
MixMasks := [3,6,12,24,48,96,192,129];

EA_LIST := [
    <"Gold", Gold>,
    <"Gold_plus_const_91", AddAffineConstant(Gold,91)>,
    <"Gold_plus_identity_linear", AddLinearOutput(Gold, IDMasks)>,
    <"Gold_plus_mixed_linear", AddLinearOutput(Gold, MixMasks)>,
    <"Output_mul_17", OutputLinear(Gold,17)>,
    <"Input_mul_5", InputLinear(Gold,5)>,
    <"Output_Frobenius_3", OutputFrobenius(Gold,3)>,
    <"Input_Frobenius_2", InputFrobenius(Gold,2)>,
    <"Full_EA_example", FullEAExample(Gold)>
];

print "=======================================================";
print "EA SANITY CHECK FOR FAST DUAL GRAPH-CODE TEST";
print "Reference: Gold x^3, modulus 285 = 0x11D";
print "Expected result: all deliberately EA-equivalent examples should test TRUE.";
print "=======================================================";

Cgold := DualGraphCode(Gold);
print "Gold dual code: length=", Length(Cgold), " dim=", Dimension(Cgold), " d=", MinimumDistance(Cgold);
print "Gold dual WD:", WDCompact(Cgold);

for item in EA_LIST do
    label := item[1];
    S := item[2];
    PrintSboxBasic(label, S, Gold);
    ok, tiso := FastDualCodeEquivalent(S, Gold);
    print "  fast dual-code isomorphic to Gold:", ok, " time:", tiso, "sec";
end for;

print "=======================================================";
print "PAIRWISE AMONG EA EXAMPLES";
print "=======================================================";
for i in [1..#EA_LIST] do
    for j in [i+1..#EA_LIST] do
        ok, tiso := FastDualCodeEquivalent(EA_LIST[i][2], EA_LIST[j][2]);
        print EA_LIST[i][1], "vs", EA_LIST[j][1], " -> ", ok, " time:", tiso;
    end for;
end for;

print "=======================================================";
print "DONE.";
print "If some EA examples return false, the fast dual-code test is not a sufficient positive CCZ test in this implementation.";
print "If all return true, it passes this sanity check.";
print "=======================================================";
