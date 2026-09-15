(* The three UCV rings, same method, different Hamiltonians:

     Rashba-Dresselhaus     NN hopping + Rashba(alpha) + Dresselhaus(beta)
     graphene               + Kane-Mele intrinsic SOC at SECOND neighbours
                              (lambda_EO, with the chirality sign nu_n)
                            + Rashba lambda_R1 at first neighbours
     silicene               + intrinsic Rashba lambda_R2 at SECOND neighbours

   All three are quadratic, so the engine handles all three unchanged: the only
   thing that differs is the H you type. This script builds the graphene and
   silicene rings, checks hermiticity, and compares the machine-generated equation
   of motion against the reference form. Pure ASCII. *)

$pass = 0; $fail = 0;
checkTrue[name_String, v_] :=
  If[TrueQ[v], $pass++; Print["  PASS  ", name],
    $fail++; Print["  FAIL  ", name, "  -> ",
      StringTake[ToString[InputForm[v]], UpTo[150]]]];

DeclareSpecies[c, "Fermion"];

ns = 8;                          (* ring sites; must be even for the A/B sign *)
nxt[n_] := Mod[n, ns] + 1;
nx2[n_] := Mod[n + 1, ns] + 1;
(* THE THESES' CHIRALITY SIGN, deliberately: alternating, inherited from the
   honeycomb lattice. Sec. IV C of the paper argues that on a convex ring every
   two-step path turns the same way, so nu is CONSTANT there and is absorbed
   into lambda_EO -- which is what scans.wl, figure4.wl, phsym.wl and
   open_family.wl do, and where every number the paper quotes comes from. This
   suite reproduces the construction as written in the theses, so it keeps
   their convention; the two are not equivalent (spectra differ by up to
   0.28 t) and that difference is the point of Sec. IV C. *)
nu[n_] := (-1)^(n - 1);
th[n_] := 2 Pi (n - 1)/ns;       (* bond angle, reference convention *)
sz[s_] := s/2;                   (* S_z eigenvalue, hbar = 1 *)

(* Hermitian conjugate of a quadratic operator expression. Hand-signing the
   h.c. is where the sign errors live -- a term with an overall i needs a minus
   that a term without one does not -- so derive it instead. All named
   parameters are assumed real. *)
hcOf[e_] := Expand[e] /. cc_. Mono[{Cre[sp_, i_], Ann[sp2_, j_]}] :>
   (cc /. Complex[re_, im_] :> Complex[re, -im]) Mono[{Cre[sp2, j], Ann[sp, i]}];
herm[e_] := Expand[NCTimes[e] + hcOf[NCTimes[e]]];

(* --- the four building blocks, each written Hermitian --------------- *)

(* A: on-site *)
hA[t0_] := t0 Sum[NCTimes[Cre[c, {n, s}], Ann[c, {n, s}]], {n, ns}, {s, {1, -1}}];

(* B: nearest-neighbour hopping with the Aharonov-Bohm phase *)
hB[t_, ph_] := t Sum[
   Exp[I ph] NCTimes[Cre[c, {n, s}], Ann[c, {nxt[n], s}]] +
   Exp[-I ph] NCTimes[Cre[c, {nxt[n], s}], Ann[c, {n, s}]],
   {n, ns}, {s, {1, -1}}];

(* C: Kane-Mele intrinsic SOC, second neighbours, spin conserving.
      i lam nu_n c^dag_n S_z c_{n+2}  + h.c.  -- the h.c. is what makes it
      Hermitian; nu_ij = -nu_ji does the same job in the standard KM form. *)
(* C without the h.c., exactly as printed in the theses *)
hClit[lam_] := I lam Sum[
   nu[n] sz[s] NCTimes[Cre[c, {n, s}], Ann[c, {nx2[n], s}]],
   {n, ns}, {s, {1, -1}}];
hC[lam_] := herm[hClit[lam]];

(* D: Rashba at first neighbours, spin flipping.
      i lamR c^dag_n (S x d_n)_z c_{n+1} + h.c.
      (S x d)_z = (1/2)(sigma_x sin th - sigma_y cos th), whose only nonzero
      entries are (up,down) = i e^{-i th}/2 and (down,up) = -i e^{i th}/2. *)
hDlit[lr_] := -(lr/2) Sum[
   s Exp[-I s th[n]] NCTimes[Cre[c, {n, s}], Ann[c, {nxt[n], -s}]],
   {n, ns}, {s, {1, -1}}];
hD[lr_] := herm[hDlit[lr]];

(* E: intrinsic Rashba of silicene, second neighbours, spin flipping *)
hElit[lr2_] := -(lr2/2) Sum[
   nu[n] s Exp[-I s th[n]] NCTimes[Cre[c, {n, s}], Ann[c, {nx2[n], -s}]],
   {n, ns}, {s, {1, -1}}];
hE[lr2_] := herm[hElit[lr2]];

ops = Flatten[Table[Ann[c, {n, s}], {n, ns}, {s, {1, -1}}]];

hGraphene = hA[t0] + hB[tt, ph] + hC[lEO] + hD[lR1];          (* graphene *)
hSilicene = hA[t0] + hB[tt, ph] + hC[lEO] + hD[lR1] + hE[lR2];   (* silicene *)

Print["=============================================="];
Print[" 1. hermiticity"];
Print["=============================================="];

Module[{m},
  m = HamiltonianMatrix[hGraphene, ops];
  checkTrue["graphene ring H is Hermitian (with h.c. restored)",
    Simplify[Max[Abs[N[(m - ConjugateTranspose[m]) /.
      {t0 -> 1/5, tt -> 1, ph -> 1/4, lEO -> 3/10, lR1 -> 1/5}]]]] < 10^-12];
  m = HamiltonianMatrix[hSilicene, ops];
  checkTrue["silicene ring H is Hermitian (with h.c. restored)",
    Simplify[Max[Abs[N[(m - ConjugateTranspose[m]) /.
      {t0 -> 1/5, tt -> 1, ph -> 1/4, lEO -> 3/10, lR1 -> 1/5, lR2 -> 1/10}]]]] < 10^-12];
  (* and the literal transcription, without the h.c., is NOT Hermitian *)
  m = HamiltonianMatrix[hA[t0] + hB[tt, ph] + hClit[lEO] + hD[lR1], ops];
  checkTrue["the intrinsic-SOC term as literally printed is NOT Hermitian",
    Max[Abs[N[(m - ConjugateTranspose[m]) /.
      {t0 -> 1/5, tt -> 1, ph -> 1/4, lEO -> 3/10, lR1 -> 1/5}]]] > 10^-3];
];

Print[""];
Print["=============================================="];
Print[" 2. the equation of motion, machine vs reference"];
Print["=============================================="];
Print["  The reference calculation writes the EOM for"];
Print["  G_{n+1 s, n s} with these Green functions on the right:"];
Print["      c_{n+1,s}   c_{n+2,s}   c_{n,s}   c_{n+3,s}"];
Print["      c_{n,-s}    c_{n+2,-s}"];
Print[""];

Module[{n0 = 3, sg = 1, cm, appear, expected, extra, missing, fmt},
  cm = Comm[Ann[c, {nxt[n0], sg}], hGraphene];
  appear = Sort[DeleteDuplicates[
     Cases[cm, Ann[c, idx_] :> idx, {0, Infinity}]]];
  fmt[{k_, s_}] := "c_{" <> ToString[k] <> "," <> If[s === 1, "+", "-"] <> "}";
  Print["  machine, for n = ", n0, " (so n+1 = ", nxt[n0], "), s = +:"];
  Print["      ", StringRiffle[fmt /@ appear, "  "]];
  expected = Sort[{{nxt[n0], sg}, {nx2[n0], sg}, {n0, sg},
                   {nxt[nx2[n0]], sg}, {n0, -sg}, {nx2[n0], -sg}}];
  Print["  reference list, same labels:"];
  Print["      ", StringRiffle[fmt /@ expected, "  "]];
  extra = Complement[appear, expected];
  missing = Complement[expected, appear];
  Print["  in machine but not in reference: ", If[extra === {}, "none", fmt /@ extra]];
  Print["  in reference but not in machine: ", If[missing === {}, "none", fmt /@ missing]];
  checkTrue["every Green function of the reference is produced by the engine",
    missing === {}];
  checkTrue["the engine finds ONE extra: the backward second-neighbour hop",
    Length[extra] === 1 && First[extra] === {Mod[nxt[n0] - 2, ns, 1], sg}];
];

Print[""];
Print["  Same test on the literal (non-Hermitian) transcription:"];

Module[{n0 = 3, sg = 1, cm, appear, expected, extra},
  cm = Comm[Ann[c, {nxt[n0], sg}], hA[t0] + hB[tt, ph] + hClit[lEO] + hD[lR1]];
  appear = Sort[DeleteDuplicates[Cases[cm, Ann[c, idx_] :> idx, {0, Infinity}]]];
  expected = Sort[{{nxt[n0], sg}, {nx2[n0], sg}, {n0, sg},
                   {nxt[nx2[n0]], sg}, {n0, -sg}, {nx2[n0], -sg}}];
  extra = Complement[appear, expected];
  Print["    extra terms: ", If[extra === {}, "none", extra]];
  checkTrue["without the h.c. the engine reproduces the reference EOM exactly",
    Sort[appear] === expected];
];

Print[""];
Print["=============================================="];
Print[" 3. the hierarchy still closes"];
Print["=============================================="];

Module[{cl, cls},
  cl = CloseEOM[{GF[Ann[c, {2, 1}], Cre[c, {1, 1}]]}, hGraphene];
  Print["  graphene ring, ", ns, " sites: ", Length[cl["Unknowns"]], " unknowns"];
  checkTrue["graphene closes at 2N", Length[cl["Unknowns"]] === 2 ns];
  cls = CloseEOM[{GF[Ann[c, {2, 1}], Cre[c, {1, 1}]]}, hSilicene];
  Print["  silicene ring, ", ns, " sites: ", Length[cls["Unknowns"]], " unknowns"];
  checkTrue["silicene closes at 2N", Length[cls["Unknowns"]] === 2 ns];
  Print["  -> second- and third-neighbour couplings extend the REACH of the"];
  Print["     equation of motion, not the SIZE of the closed system: periodicity"];
  Print["     still caps it at one Green function per (site, spin)."];
];

Print[""];
Print["=============================================="];
Print[" 4. spectra are real and respond to each coupling"];
Print["=============================================="];

Module[{subs, ev, spread},
  subs = {t0 -> 0, tt -> 1, ph -> 1/4};
  ev[h_, extra_] := Sort[N[Eigenvalues[
     N[HamiltonianMatrix[h, ops] /. Join[subs, extra]]]]];
  spread[a_, b_] := Max[Abs[a - b]];
  checkTrue["graphene spectrum real",
    Max[Abs[Im[ev[hGraphene, {lEO -> 3/10, lR1 -> 1/5}]]]] < 10^-10];
  checkTrue["silicene spectrum real",
    Max[Abs[Im[ev[hSilicene, {lEO -> 3/10, lR1 -> 1/5, lR2 -> 1/10}]]]] < 10^-10];
  checkTrue["intrinsic SOC moves the spectrum",
    spread[ev[hGraphene, {lEO -> 0, lR1 -> 0}],
           ev[hGraphene, {lEO -> 3/10, lR1 -> 0}]] > 10^-3];
  checkTrue["Rashba moves it further",
    spread[ev[hGraphene, {lEO -> 3/10, lR1 -> 0}],
           ev[hGraphene, {lEO -> 3/10, lR1 -> 1/5}]] > 10^-3];
  checkTrue["silicene's extra term is not redundant with the others",
    spread[ev[hSilicene, {lEO -> 3/10, lR1 -> 1/5, lR2 -> 0}],
           ev[hSilicene, {lEO -> 3/10, lR1 -> 1/5, lR2 -> 1/10}]] > 10^-3];
];

Print[""];
Print["=============================================="];
Print["  PASSED: ", $pass, "     FAILED: ", $fail];
Print["=============================================="];
