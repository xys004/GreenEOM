(* ::Package:: *)

(* GreenEOM -- Green-function equations of motion from a second-quantised
   Hamiltonian.

   Given H written with creation/annihilation operators, the package
     1. computes [A,H] using the fundamental (anti)commutators,
     2. writes the Fourier-transformed equation of motion for <<A;B>>_w,
     3. closes the hierarchy automatically (exactly, if H is quadratic),
     4. solves the resulting linear algebraic system,
     5. builds G^r, G^a and G^< in equilibrium.

   Source is pure ASCII: every special character uses a \[Name] escape so the
   file survives transport through shells and non-UTF8 code pages.

   Conventions
     G(t,t') = -i <T_c A(t) B(t')>
     i d_t G = delta(t,t') <[A,B]_-+> + <<[A,H];B>>
     Fourier: w G(w) = norm <[A,B]_-+> + <<[A,H];B>>_w
     with norm = 1/Sqrt[2 Pi] by default (reference convention); use 1/(2 Pi) for
     Zubarev's, or 1 for the bare resolvent.
*)

BeginPackage["GreenEOM`"];

Ann::usage =
  "Ann[sp, idx] is an annihilation operator of species sp with index idx. \
idx may be any expression, e.g. a site number or a list {site, spin}.";
Cre::usage =
  "Cre[sp, idx] is the corresponding creation operator.";
DeclareSpecies::usage =
  "DeclareSpecies[sp, \"Fermion\"] or DeclareSpecies[sp, \"Boson\"] sets the \
statistics of species sp. Species default to \"Fermion\".";
SpeciesStatistics::usage =
  "SpeciesStatistics[sp] returns \"Fermion\" or \"Boson\".";

Mono::usage =
  "Mono[{o1, o2, ...}] is an ordered operator monomial. Mono[{}] is the identity.";
NCTimes::usage =
  "NCTimes[a, b, ...] is the non-commutative product, linear in every slot. \
Operator expressions may also be written with ** .";
NOrder::usage =
  "NOrder[expr] normal-orders expr (creation operators to the left), producing \
contraction terms from the fundamental (anti)commutators.";
Comm::usage = "Comm[a, b] is the commutator a b - b a, normal-ordered.";
AComm::usage = "AComm[a, b] is the anticommutator a b + b a, normal-ordered.";
OperatorParity::usage =
  "OperatorParity[expr] is 0 for an even (bosonic) operator expression, 1 for \
an odd (fermionic) one, and $Failed if expr mixes parities.";

GF::usage =
  "GF[A, B] denotes the Fourier-transformed Green function <<A;B>>_w. It is \
linear in both slots.";
GFOrder::usage = "GFOrder[GF[A,B]] is the number of operators in the left slot.";
Expect::usage =
  "Expect[expr] denotes the thermal expectation value <expr> of an operator \
expression that survived normal ordering.";

EOMEquation::usage =
  "EOMEquation[GF[A,B], H] returns the Fourier-space equation of motion for \
<<A;B>>_w as an Equal.";
CloseEOM::usage =
  "CloseEOM[seeds, H] repeatedly applies the equation of motion starting from \
the seed Green functions until no new ones appear, and returns an Association \
with keys \"Equations\", \"Unknowns\" and \"Frequency\".";
SolveGF::usage =
  "SolveGF[closure] solves the closed linear system and returns a list of \
replacement rules for the Green functions.";
GFMatrix::usage =
  "GFMatrix[closure] returns {M, v} such that M . G == v, where G is the vector \
of unknowns in the order given by closure[\"Unknowns\"]. M is (w I - H) in \
disguise.";

SpectralDecomposition::usage =
  "SpectralDecomposition[g, w] writes the rational function g of w as a sum of \
simple poles and returns {{pole, residue}, ...}.";
GLesser::usage =
  "GLesser[g, w, beta, eta] builds the lesser Green function in thermal \
equilibrium from the retarded/advanced pair obtained by w -> w +- I eta.";

\[Omega]::usage = "\\[Omega] is the default frequency symbol used by the equations of motion.";
FrequencySymbol::usage = "Option: the symbol used for the frequency. Default \\[Omega].";
SourceNormalization::usage =
  "Option: prefactor of the inhomogeneous term. Default 1/Sqrt[2 Pi].";
SourceTime::usage = "Option: the time t' in the source phase Exp[-I w t']. Default 0.";
MaxOperators::usage =
  "Option: highest number of operators allowed in the left slot of a Green \
function. Default 1, which is exact for quadratic Hamiltonians.";
DecouplingRule::usage =
  "Option: a function applied to any GF exceeding MaxOperators, e.g. a \
mean-field factorisation. Default None, which makes CloseEOM stop with an error.";
MaxEquations::usage = "Option: safety cap on the number of generated equations. Default 2000.";
ParticleStatistics::usage = "Option: \"Fermion\" or \"Boson\", for GLesser.";

CloseEOM::open =
  "The equation of motion generated `1`, which has `2` operators in its left \
slot (MaxOperators -> `3`). The hierarchy does not close: H is not quadratic. \
Supply a DecouplingRule to truncate it.";
CloseEOM::runaway =
  "More than `1` equations were generated. The index set is probably infinite: \
build H over an explicit finite set of sites.";

Begin["`Private`"];

(* ------------------------------------------------------------------ *)
(* statistics                                                          *)
(* ------------------------------------------------------------------ *)

$statistics = <||>;

DeclareSpecies[sp_, s_String] /; MemberQ[{"Fermion", "Boson"}, s] :=
  ($statistics[sp] = s; s);

SpeciesStatistics[sp_] := Lookup[$statistics, sp, "Fermion"];

fermionQ[Ann[sp_, _]] := SpeciesStatistics[sp] === "Fermion";
fermionQ[Cre[sp_, _]] := SpeciesStatistics[sp] === "Fermion";
fermionQ[_] := False;

(* ------------------------------------------------------------------ *)
(* Kronecker delta on possibly composite, possibly symbolic indices     *)
(* ------------------------------------------------------------------ *)

kron[i_List, j_List] /; Length[i] === Length[j] := Times @@ MapThread[kron, {i, j}];
kron[i_List, j_List] := 0;
kron[i_, j_] := KroneckerDelta[i, j];

(* ------------------------------------------------------------------ *)
(* monomials                                                           *)
(* ------------------------------------------------------------------ *)

Mono[{}] := 1;

scalarQ[e_] := FreeQ[e, _Ann | _Cre | _Mono];

ncd[x_?scalarQ, y_] := x y;
ncd[x_, y_?scalarQ] := x y;
ncd[c1_. Mono[l1_], c2_. Mono[l2_]] := c1 c2 Mono[Join[l1, l2]];

ncPair[a_, b_] := Distribute[ncd[Expand[a], Expand[b]]];

(* toMono must be idempotent: operators already sitting inside a Mono must not
   be lifted a second time, or we build Mono[{Mono[{...}], ...}]. Existing
   monomials are therefore parked behind an opaque placeholder while the bare
   operators are lifted. *)
toMono[e_] :=
  Module[{x, existing, park, unpark},
    x = Expand[e];
    existing = DeleteDuplicates[Cases[x, _Mono, {0, Infinity}]];
    park = MapIndexed[#1 -> monoSlot[First[#2]] &, existing];
    unpark = MapIndexed[monoSlot[First[#2]] -> #1 &, existing];
    x = x /. park;
    x = x /. op : (_Ann | _Cre) :> Mono[{op}];
    x = x /. unpark;
    x = x //. {
       NonCommutativeMultiply[a_] :> a,
       NonCommutativeMultiply[a_, b__] :> ncPair[a, NonCommutativeMultiply[b]]
    };
    Expand[x]
  ];

NCTimes[a_] := toMono[a];
NCTimes[a_, b__] := Fold[ncPair, toMono[a], toMono /@ {b}];

(* ------------------------------------------------------------------ *)
(* normal ordering                                                     *)
(*   x y = swapSign[x,y] y x + kernel[x,y]                             *)
(* ------------------------------------------------------------------ *)

opKey[Cre[sp_, i_]] := {0, sp, i};
opKey[Ann[sp_, i_]] := {1, sp, i};

swapSign[x_, y_] := If[fermionQ[x] && fermionQ[y], -1, 1];

kernel[Ann[sp_, i_], Cre[sp_, j_]] := kron[i, j];
kernel[Cre[sp_, i_], Ann[sp_, j_]] :=
  If[SpeciesStatistics[sp] === "Fermion", kron[i, j], -kron[i, j]];
kernel[_, _] := 0;

firstInversion[l_List] :=
  Catch[
    Do[If[! OrderedQ[{opKey[l[[i]]], opKey[l[[i + 1]]]}], Throw[i]],
      {i, Length[l] - 1}];
    None
  ];

noMono[l_List] := noMono[l] =
  Module[{p, x, y},
    p = firstInversion[l];
    If[p === None,
      If[AnyTrue[Range[Length[l] - 1],
           (l[[#]] === l[[# + 1]] && fermionQ[l[[#]]]) &],
        0,
        Mono[l]],
      x = l[[p]]; y = l[[p + 1]];
      Expand[
        swapSign[x, y] noMono[ReplacePart[l, {p -> y, p + 1 -> x}]] +
        kernel[x, y] noMono[Delete[l, {{p}, {p + 1}}]]
      ]
    ]
  ];

NOrder[e_] := Expand[Expand[toMono[e]] /. Mono[l_] :> noMono[l]];

Comm[a_, b_] := NOrder[NCTimes[a, b] - NCTimes[b, a]];
AComm[a_, b_] := NOrder[NCTimes[a, b] + NCTimes[b, a]];

OperatorParity[e_] :=
  Module[{x, ps},
    x = Expand[toMono[e]];
    If[scalarQ[x], Return[0]];
    ps = Union[Cases[x, Mono[l_] :> Mod[Count[l, _?fermionQ], 2], {0, Infinity}]];
    If[Length[ps] === 1, First[ps], $Failed]
  ];

(* ------------------------------------------------------------------ *)
(* Green functions: linear in both slots                                *)
(* ------------------------------------------------------------------ *)

GF[0, _] := 0;
GF[_, 0] := 0;
GF[a_Plus, b_] := Total[GF[#, b] & /@ (List @@ a)];
GF[a_, b_Plus] := Total[GF[a, #] & /@ (List @@ b)];
GF[c_ x_Mono, b_] /; scalarQ[c] := c GF[x, b];
GF[a_, c_ y_Mono] /; scalarQ[c] := c GF[a, y];
GF[a : (_Ann | _Cre), b_] := GF[Mono[{a}], b];
GF[a_, b : (_Ann | _Cre)] := GF[a, Mono[{b}]];

GFOrder[GF[Mono[l_], _]] := Length[l];
GFOrder[GF[a_, _]] := If[scalarQ[a], 0, $Failed];

Expect[e_Plus] := Total[Expect /@ (List @@ e)];
Expect[c_ x_Mono] /; scalarQ[c] := c Expect[x];
Expect[c_?scalarQ] := c;

expectify[e_] := Expand[e] /. m_Mono :> Expect[m];

(* ------------------------------------------------------------------ *)
(* equation of motion                                                   *)
(* ------------------------------------------------------------------ *)

Options[EOMEquation] = {
  FrequencySymbol -> \[Omega],
  SourceNormalization -> 1/Sqrt[2 Pi],
  SourceTime -> 0
};

eomRHS[g : GF[a_, b_], h_, norm_, t0_, w_] :=
  Module[{pa, pb, brk, src},
    pa = OperatorParity[a];
    pb = OperatorParity[b];
    brk = If[pa === 1 && pb === 1, AComm[a, b], Comm[a, b]];
    src = norm Exp[-I w t0] expectify[brk];
    Expand[src + GF[Comm[a, h], b]]
  ];

EOMEquation[g : GF[_, _], h_, opts : OptionsPattern[]] :=
  Module[{w, norm, t0},
    w = OptionValue[FrequencySymbol];
    norm = OptionValue[SourceNormalization];
    t0 = OptionValue[SourceTime];
    w g == eomRHS[g, h, norm, t0, w]
  ];

(* ------------------------------------------------------------------ *)
(* hierarchy closure                                                    *)
(* ------------------------------------------------------------------ *)

Options[CloseEOM] = Join[Options[EOMEquation], {
  MaxOperators -> 1,
  DecouplingRule -> None,
  MaxEquations -> 2000
}];

CloseEOM[seeds_List, h_, opts : OptionsPattern[]] :=
  Module[{w, norm, t0, maxOps, dec, cap, pending, done, eqs, g, eq, news, bad},
    w = OptionValue[FrequencySymbol];
    norm = OptionValue[SourceNormalization];
    t0 = OptionValue[SourceTime];
    maxOps = OptionValue[MaxOperators];
    dec = OptionValue[DecouplingRule];
    cap = OptionValue[MaxEquations];

    pending = DeleteDuplicates[seeds];
    done = {};
    eqs = {};

    While[pending =!= {},
      g = First[pending];
      pending = Rest[pending];
      If[MemberQ[done, g], Continue[]];
      AppendTo[done, g];

      eq = w g == eomRHS[g, h, norm, t0, w];

      If[dec =!= None,
        bad = Select[DeleteDuplicates[Cases[eq, _GF, {0, Infinity}]],
                TrueQ[GFOrder[#] > maxOps] &];
        If[bad =!= {}, eq = Expand[eq /. (Rule[#, dec[#]] & /@ bad)]]
      ];

      news = DeleteDuplicates[Cases[eq, _GF, {0, Infinity}]];
      bad = Select[news, TrueQ[GFOrder[#] > maxOps] &];
      If[bad =!= {},
        Message[CloseEOM::open, First[bad], GFOrder[First[bad]], maxOps];
        Return[$Failed]
      ];

      AppendTo[eqs, eq];
      pending = Join[pending, Complement[news, done]];

      If[Length[eqs] > cap,
        Message[CloseEOM::runaway, cap];
        Return[$Failed]
      ]
    ];

    <|"Equations" -> eqs, "Unknowns" -> done, "Frequency" -> w|>
  ];

SolveGF[c_Association] := Solve[c["Equations"], c["Unknowns"]];

GFMatrix[c_Association] :=
  Module[{ca},
    ca = CoefficientArrays[
       (#[[1]] - #[[2]]) & /@ c["Equations"], c["Unknowns"]];
    {Normal[ca[[2]]], -Normal[ca[[1]]]}
  ];

(* ------------------------------------------------------------------ *)
(* equilibrium post-processing                                          *)
(* ------------------------------------------------------------------ *)

SpectralDecomposition[g_, w_] :=
  Module[{ap, terms},
    ap = Apart[Together[g], w];
    terms = If[Head[ap] === Plus, List @@ ap, {ap}];
    DeleteCases[
      Map[
        Function[term,
          Module[{den, pole},
            den = Denominator[Together[term]];
            If[Exponent[den, w] =!= 1, Nothing,
              pole = w /. First[Solve[den == 0, w]];
              {pole, Simplify[Residue[term, {w, pole}]]}
            ]
          ]
        ], terms],
      Nothing]
  ];

Options[GLesser] = {ParticleStatistics -> "Fermion"};

GLesser[g_, w_, beta_, eta_, OptionsPattern[]] :=
  Module[{gr, ga},
    gr = g /. w -> w + I eta;
    ga = g /. w -> w - I eta;
    Switch[OptionValue[ParticleStatistics],
      "Fermion", -(gr - ga)/(Exp[beta w] + 1),
      "Boson", (gr - ga)/(Exp[beta w] - 1)
    ]
  ];

End[];
EndPackage[];
