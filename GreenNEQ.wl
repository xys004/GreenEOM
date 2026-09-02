(* ::Package:: *)

(* GreenNEQ -- steady-state non-equilibrium (Keldysh) layer on top of GreenEOM.
   Pure ASCII source.

   What changes with respect to equilibrium
     Equilibrium ties G^< to G^r through one Fermi function. With several
     reservoirs at different chemical potentials that link is cut and G^<
     becomes an independent object,  G^< = G^r Sigma^< G^a.
     Reservoirs must also be continua: discrete levels only give poles, never
     a finite level width and never a current.

   Scope
     Quadratic H, steady state, no initial correlations. That case is exactly
     solvable and is what this module does. Interacting Keldysh (vertex
     corrections, Kadanoff-Baym time propagation) is NOT here.

   Conventions
     H_ij defined by  [c_i, H] = Sum_j H_ij c_j
     G^r(w) = (w I - H - Sigma^r)^-1,  G^a = (G^r)^dagger for real w
     Gamma_alpha = i (Sigma^r_alpha - Sigma^a_alpha)
     Sigma^< = i Sum_alpha Gamma_alpha f_alpha
     All symbolic parameters are assumed real.
*)

BeginPackage["GreenNEQ`", {"GreenEOM`"}];

HamiltonianMatrix::usage =
  "HamiltonianMatrix[H, ops] returns the single-particle matrix H_ij defined by \
[c_i, H] = Sum_j H_ij c_j, where ops is the ordered list of annihilation \
operators. It fails if H is not quadratic.";
Downfold::usage =
  "Downfold[h, keep, w] eliminates the degrees of freedom outside `keep` by a \
Schur complement and returns {h_keep, Sigma[w]}. The self-energy is the Schur \
complement -- it is derived, not postulated.";

CProd::usage = "CProd[a, b, ...] is a product of contour-ordered quantities.";
ContourRetarded::usage = "ContourRetarded[CProd[...]] applies the Langreth rule (AB)^r = A^r B^r.";
ContourAdvanced::usage = "ContourAdvanced[CProd[...]] applies (AB)^a = A^a B^a.";
ContourLesser::usage = "ContourLesser[CProd[...]] applies (AB)^< = A^r B^< + A^< B^a.";
ContourGreater::usage = "ContourGreater[CProd[...]] applies (AB)^> = A^r B^> + A^> B^a.";
MDot::usage = "MDot[a, b, ...] is the (non-commutative) product of real-time components.";
Rt::usage = "Rt[x] is the retarded component of x.";
Ad::usage = "Ad[x] is the advanced component of x.";
Ls::usage = "Ls[x] is the lesser component of x.";
Gt::usage = "Gt[x] is the greater component of x.";

Dagger::usage =
  "Dagger[m] is the Hermitian conjugate of m, assuming every symbol in m is real.";
WideBandSelfEnergy::usage =
  "WideBandSelfEnergy[gamma] returns the retarded self-energy -I gamma/2 of a \
wide-band reservoir with coupling matrix gamma.";
CouplingMatrix::usage =
  "CouplingMatrix[sigmaR] returns Gamma = I (Sigma^r - Sigma^a).";

RetardedG::usage = "RetardedG[h, sigmaR, w] = Inverse[w I - h - sigmaR].";
AdvancedG::usage = "AdvancedG[gr] = Dagger[gr].";
LesserSigma::usage =
  "LesserSigma[{gamma1, gamma2, ...}, {f1, f2, ...}] = I Sum_a Gamma_a f_a.";
LesserG::usage = "LesserG[gr, sigmaLess] = gr . sigmaLess . Dagger[gr].";
GreaterG::usage = "GreaterG[gr, sigmaGreater] = gr . sigmaGreater . Dagger[gr].";

Transmission::usage =
  "Transmission[gr, gammaL, gammaR] = Tr[gammaL . gr . gammaR . Dagger[gr]].";
LandauerCurrent::usage =
  "LandauerCurrent[tfun, muL, muR, beta, {wmin, wmax}] integrates \
T(w)(f_L - f_R) dw/(2 Pi) numerically.";
MeirWingreenCurrent::usage =
  "MeirWingreenCurrent[gr, glss, gamma, f, w] is the integrand \
(1/2 Pi) Tr[Gamma (G^< + f (G^r - G^a))] of the Meir-Wingreen current into \
the corresponding lead.";
Occupation::usage =
  "Occupation[glss, j] = -I G^<_jj / (2 Pi), the integrand of the occupation \
of site j.";
FermiF::usage = "FermiF[w, mu, beta] = 1/(Exp[beta (w - mu)] + 1).";

HamiltonianMatrix::notquad =
  "[c_`1`, H] contains `2`, which is not a single operator: H is not quadratic, \
so it has no single-particle matrix.";

Begin["`Private`"];

(* ------------------------------------------------------------------ *)
(* single-particle matrix, straight out of the verified commutator      *)
(* ------------------------------------------------------------------ *)

(* Catch/Throw, not Return: Return would only leave the inner Module of the
   Table and quietly hand back a row of $Failed. *)
HamiltonianMatrix[h_, ops_List] :=
  Catch[
    Table[
      Module[{cm, terms, bad},
        cm = Expand[Comm[ops[[i]], h]];
        terms = If[Head[cm] === Plus, List @@ cm, {cm}];
        bad = Select[terms, ! MatchQ[#, _. Mono[{_Ann}]] && # =!= 0 &];
        If[bad =!= {},
          Message[HamiltonianMatrix::notquad, i, First[bad]];
          Throw[$Failed, "GreenNEQnotquad"]];
        Table[Coefficient[cm, Mono[{ops[[j]]}]], {j, Length[ops]}]
      ], {i, Length[ops]}],
    "GreenNEQnotquad"];

(* ------------------------------------------------------------------ *)
(* lead elimination = Schur complement                                  *)
(* ------------------------------------------------------------------ *)

Downfold[h_, keep_List, w_] :=
  Module[{all, drop, hkk, hke, hek, hee, sigma},
    all = Range[Length[h]];
    drop = Complement[all, keep];
    hkk = h[[keep, keep]];
    If[drop === {}, Return[{hkk, ConstantArray[0, {Length[keep], Length[keep]}]}]];
    hke = h[[keep, drop]];
    hek = h[[drop, keep]];
    hee = h[[drop, drop]];
    sigma = hke . Inverse[w IdentityMatrix[Length[drop]] - hee] . hek;
    {hkk, sigma}
  ];

(* ------------------------------------------------------------------ *)
(* Langreth rules                                                       *)
(* ------------------------------------------------------------------ *)

(* MDot must be linear as well as flat, or the Langreth recursion leaves
   MDot[x, MDot[..] + MDot[..]] unexpanded. Distribute first, then flatten. *)
MDot[a___, x_Plus, b___] := Total[MDot[a, #, b] & /@ (List @@ x)];
MDot[a___, MDot[b___], c___] := MDot[a, b, c];

ContourRetarded[CProd[a_]] := Rt[a];
ContourAdvanced[CProd[a_]] := Ad[a];
ContourLesser[CProd[a_]] := Ls[a];
ContourGreater[CProd[a_]] := Gt[a];

ContourRetarded[CProd[a_, b__]] := MDot[Rt[a], ContourRetarded[CProd[b]]];
ContourAdvanced[CProd[a_, b__]] := MDot[Ad[a], ContourAdvanced[CProd[b]]];
ContourLesser[CProd[a_, b__]] :=
  MDot[Rt[a], ContourLesser[CProd[b]]] + MDot[Ls[a], ContourAdvanced[CProd[b]]];
ContourGreater[CProd[a_, b__]] :=
  MDot[Rt[a], ContourGreater[CProd[b]]] + MDot[Gt[a], ContourAdvanced[CProd[b]]];

ContourRetarded[a_Plus] := ContourRetarded /@ a;
ContourAdvanced[a_Plus] := ContourAdvanced /@ a;
ContourLesser[a_Plus] := ContourLesser /@ a;
ContourGreater[a_Plus] := ContourGreater /@ a;

(* ------------------------------------------------------------------ *)
(* Keldysh assembly                                                     *)
(* ------------------------------------------------------------------ *)

conjR[e_] := e /. Complex[re_, im_] :> Complex[re, -im];
Dagger[m_] := Transpose[conjR[m]];

WideBandSelfEnergy[gamma_] := -I gamma/2;
CouplingMatrix[sigmaR_] := I (sigmaR - Dagger[sigmaR]);

RetardedG[h_, sigmaR_, w_] :=
  Inverse[w IdentityMatrix[Length[h]] - h - sigmaR];
AdvancedG[gr_] := Dagger[gr];

LesserSigma[gammas_List, fs_List] := I Total[MapThread[#1 #2 &, {gammas, fs}]];

LesserG[gr_, sigmaLess_] := gr . sigmaLess . Dagger[gr];
GreaterG[gr_, sigmaGreater_] := gr . sigmaGreater . Dagger[gr];

FermiF[w_, mu_, beta_] := 1/(Exp[beta (w - mu)] + 1);

Transmission[gr_, gammaL_, gammaR_] := Tr[gammaL . gr . gammaR . Dagger[gr]];

LandauerCurrent[tfun_, muL_, muR_, beta_, {wmin_, wmax_}] :=
  NIntegrate[
    tfun[w] (FermiF[w, muL, beta] - FermiF[w, muR, beta])/(2 Pi),
    {w, wmin, wmax}];

MeirWingreenCurrent[gr_, glss_, gamma_, f_] :=
  Tr[gamma . (glss + f (gr - Dagger[gr]))]/(2 Pi);

Occupation[glss_, j_] := -I glss[[j, j]]/(2 Pi);

End[];
EndPackage[];
