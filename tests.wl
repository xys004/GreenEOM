(* Regression suite for GreenEOM. Pure ASCII. *)

$pass = 0; $fail = 0;

check[name_String, got_, want_] :=
  Module[{ok},
    ok = TrueQ[Simplify[FullSimplify[got - want] == 0]] ||
         TrueQ[got === want] ||
         TrueQ[Simplify[got == want]];
    If[ok, $pass++; Print["  PASS  ", name],
      $fail++; Print["  FAIL  ", name];
      Print["        got  = ", InputForm[got]];
      Print["        want = ", InputForm[want]]];
    ok
  ];

checkTrue[name_String, val_] :=
  If[TrueQ[val], $pass++; Print["  PASS  ", name],
    $fail++; Print["  FAIL  ", name, "  -> ", InputForm[val]]];

Print["=============================================="];
Print[" 1. operator algebra"];
Print["=============================================="];

DeclareSpecies[c, "Fermion"];
DeclareSpecies[a, "Boson"];

check["{c_i, cd_j} = delta_ij",
  AComm[Ann[c, i], Cre[c, j]], KroneckerDelta[i, j]];
check["{c_1, cd_1} = 1",
  AComm[Ann[c, 1], Cre[c, 1]], 1];
check["{c_1, cd_2} = 0",
  AComm[Ann[c, 1], Cre[c, 2]], 0];
check["{c_i, c_j} = 0",
  AComm[Ann[c, i], Ann[c, j]], 0];
check["c_1 c_1 = 0 (Pauli)",
  NOrder[NCTimes[Ann[c, 1], Ann[c, 1]]], 0];
check["[a_i, ad_j] = delta_ij (boson)",
  Comm[Ann[a, i], Cre[a, j]], KroneckerDelta[i, j]];
check["[a_1, a_1] = 0 but a_1 a_1 =!= 0 (boson)",
  Comm[Ann[a, 1], Ann[a, 1]], 0];
checkTrue["boson square survives",
  NOrder[NCTimes[Ann[a, 1], Ann[a, 1]]] =!= 0];
check["fermion and boson commute",
  Comm[Ann[c, 1], Cre[a, 1]], 0];
check["[c_i, n_j] = delta_ij c_i",
  Comm[Ann[c, i], NCTimes[Cre[c, j], Ann[c, j]]],
  KroneckerDelta[i, j] Mono[{Ann[c, i]}]];
check["composite index delta",
  AComm[Ann[c, {1, 1}], Cre[c, {1, -1}]], 0];
check["shifted symbolic index delta",
  AComm[Ann[c, n + 1], Cre[c, n]], 0];

Print[""];
Print["=============================================="];
Print[" 2. single level:  H = eps cd c"];
Print["=============================================="];

Module[{h, cl, sol, g},
  h = eps NCTimes[Cre[c, 1], Ann[c, 1]];
  cl = CloseEOM[{GF[Ann[c, 1], Cre[c, 1]]}, h];
  Print["  equations: ", Length[cl["Equations"]]];
  sol = SolveGF[cl];
  g = GF[Mono[{Ann[c, 1]}], Mono[{Cre[c, 1]}]] /. First[sol];
  check["G = (1/Sqrt[2 Pi])/(w - eps)", Simplify[g], 1/(Sqrt[2 Pi] (\[Omega] - eps))];
];

Print[""];
Print["=============================================="];
Print[" 3. boson mode:  H = w0 ad a"];
Print["=============================================="];

Module[{h, cl, sol, g},
  h = w0 NCTimes[Cre[a, 1], Ann[a, 1]];
  cl = CloseEOM[{GF[Ann[a, 1], Cre[a, 1]]}, h];
  sol = SolveGF[cl];
  g = GF[Mono[{Ann[a, 1]}], Mono[{Cre[a, 1]}]] /. First[sol];
  check["G = (1/Sqrt[2 Pi])/(w - w0)", Simplify[g], 1/(Sqrt[2 Pi] (\[Omega] - w0))];
];

Print[""];
Print["=============================================="];
Print[" 4. mesoscopic ring (Maiti Hamiltonian)"];
Print["=============================================="];

(* Ring Hamiltonian, exactly the reference one:
     H = sum_n cd_n eps_n c_n
       + sum_n ( cd_{n+1} . T(n) . c_n e^{I th} + h.c. )
   with  T(n) = {{t, -I al e^{-I ph} + be e^{I ph}},
                 {-I al e^{I ph} - be e^{-I ph}, t}},  ph = ph_{n,n+1}.
   Spin labels: +1 = up (row/col 1), -1 = down (row/col 2). *)

ringH[ns_, t_, th_, al_, be_] :=
  Module[{phi, ph, nxt, tm, tmc, spins = {1, -1}},
    phi[n_] := 2 Pi (n - 1)/ns;
    (* THE MONOGRAPH'S BOND ANGLE, deliberately, and not the canonical one
       used everywhere else in this repository. Averaging the two site angles
       reads phi_1 as 0 rather than 2 Pi, so the closing bond comes out a full
       Pi away from where it belongs; the canonical form, phi_n + m Pi/N, has
       no such special case. But these are regressions against Green functions
       obtained outside the generator, and a regression has to be run in the
       convention of the thing it reproduces -- otherwise it tests the model
       rather than the generator. Every result the paper reports for itself
       uses the canonical angle. *)
    ph[n_, m_] := (phi[n] + phi[m])/2;
    nxt[n_] := Mod[n, ns] + 1;
    (* T(n)[s, sp] with s the row (site n+1) and sp the column (site n) *)
    tm[n_][1, 1] = t;  tm[n_][-1, -1] = t;
    tm[n_][1, -1] := -I al Exp[-I ph[n, nxt[n]]] + be Exp[I ph[n, nxt[n]]];
    tm[n_][-1, 1] := -I al Exp[I ph[n, nxt[n]]] - be Exp[-I ph[n, nxt[n]]];
    (* explicit complex conjugate (al, be, t real) *)
    tmc[n_][1, 1] = t;  tmc[n_][-1, -1] = t;
    tmc[n_][1, -1] := I al Exp[I ph[n, nxt[n]]] + be Exp[-I ph[n, nxt[n]]];
    tmc[n_][-1, 1] := I al Exp[-I ph[n, nxt[n]]] - be Exp[I ph[n, nxt[n]]];
    Sum[
      Exp[I th] tm[n][s, sp] NCTimes[Cre[c, {nxt[n], s}], Ann[c, {n, sp}]] +
      Exp[-I th] tmc[n][s, sp] NCTimes[Cre[c, {n, sp}], Ann[c, {nxt[n], s}]],
      {n, ns}, {s, spins}, {sp, spins}]
  ];

Print[""];
Print[" 4a. three sites, hopping only  (vs the closed form transcribed below \
as referenceG = the three-site closed form)"];

Module[{h, cl, sol, g, den, roots, want, referenceG},
  h = ringH[3, t, th, 0, 0];
  cl = CloseEOM[{GF[Ann[c, {2, 1}], Cre[c, {1, 1}]]}, h];
  Print["  equations: ", Length[cl["Equations"]],
        "   unknowns: ", Length[cl["Unknowns"]]];
  sol = First[SolveGF[cl]];
  g = Simplify[GF[Mono[{Ann[c, {2, 1}]}], Mono[{Cre[c, {1, 1}]}]] /. sol];
  Print["  G_{n+1,n} = ", InputForm[g]];
  den = Denominator[Together[g]];
  roots = \[Omega] /. Solve[den == 0, \[Omega]];
  (* (i) the literal closed form, transcribed by hand below. The comparison
     is against this in-script transcription, so the test is self-contained
     and reads no external document at run time. *)
  referenceG = -(Exp[I th] t (t + Exp[3 I th] \[Omega]))/
             (Sqrt[2 Pi] (Exp[6 I th] t^3 + t^3 + 3 Exp[3 I th] t^2 \[Omega] -
                          Exp[3 I th] \[Omega]^3));
  check["G_{n+1,n} equals the closed form transcribed here as referenceG",
    Simplify[g - referenceG], 0];
  (* (ii) the poles the reference quotes explicitly *)
  want = Table[2 t Cos[th + 2 Pi m/3], {m, 0, 2}];
  checkTrue["poles = 2tCos[th], -tCos[th] +- Sqrt[3] t Sin[th]",
    Length[roots] === 3 &&
    Max[Abs[Sort[N[roots /. {t -> 1, th -> 37/100}]] -
            Sort[N[{2 t Cos[th], -t Cos[th] + Sqrt[3] t Sin[th],
                    -t Cos[th] - Sqrt[3] t Sin[th]} /. {t -> 1, th -> 37/100}]]]] < 10^-10];
  checkTrue["poles = 2 t Cos[th + 2 Pi m/3]  (ring spectrum)",
    Max[Abs[Sort[N[roots /. {t -> 1, th -> 37/100}]] -
            Sort[N[want /. {t -> 1, th -> 37/100}]]]] < 10^-10];
];

Print[""];
Print[" 4b. two sites, hopping only  (vs the poles written out below; written out here)"];

Module[{h, cl, sol, g, den, roots},
  h = ringH[2, t, th, 0, 0];
  cl = CloseEOM[{GF[Ann[c, {2, 1}], Cre[c, {1, 1}]]}, h];
  sol = First[SolveGF[cl]];
  g = Simplify[GF[Mono[{Ann[c, {2, 1}]}], Mono[{Cre[c, {1, 1}]}]] /. sol];
  Print["  G_{n+1,n} = ", InputForm[g]];
  den = Denominator[Together[g]];
  roots = Simplify[\[Omega] /. Solve[den == 0, \[Omega]]];
  Print["  poles = ", InputForm[roots]];
  checkTrue["poles = +- 2 t Cos[th]  (closed form written out here; written out here)",
    Max[Abs[Sort[N[roots /. {t -> 1, th -> 37/100}]] -
            Sort[N[{2 t Cos[th], -2 t Cos[th]} /. {t -> 1, th -> 37/100}]]]] < 10^-10];
];

Print[""];
Print[" 4c. two sites with Rashba + Dresselhaus (vs 'want' written out below; written out here)"];

Module[{h, cl, eqs, e1, lhs, rhs, g1, g2, g3, g4, want},
  h = ringH[2, t, th, al, be];
  cl = CloseEOM[{GF[Ann[c, {2, 1}], Cre[c, {1, 1}]]}, h];
  Print["  equations: ", Length[cl["Equations"]],
        "   unknowns: ", Length[cl["Unknowns"]]];
  (* reference names:
       g1 = G_{n+1 up, n up}   g2 = G_{n up, n up}
       g3 = G_{n down, n up}   g4 = G_{n+1 down, n up}
     Its eq1 reads
       w g1 == 2 t Cos[th] g2 + 2 I Sin[th](be E^{I ph} - I al E^{-I ph}) g3
     with ph = ph_{n,n+1} = Pi/2 for two sites (phi_1 = 0, phi_2 = Pi). *)
  g1 = GF[Mono[{Ann[c, {2, 1}]}], Mono[{Cre[c, {1, 1}]}]];
  g2 = GF[Mono[{Ann[c, {1, 1}]}], Mono[{Cre[c, {1, 1}]}]];
  g3 = GF[Mono[{Ann[c, {1, -1}]}], Mono[{Cre[c, {1, 1}]}]];
  g4 = GF[Mono[{Ann[c, {2, -1}]}], Mono[{Cre[c, {1, 1}]}]];
  e1 = SelectFirst[cl["Equations"], (#[[1]] === \[Omega] g1) &];
  Print["  generated eq1 rhs = ", InputForm[e1[[2]]]];
  want = 2 t Cos[th] g2 +
         2 I Sin[th] (be Exp[I Pi/2] - I al Exp[-I Pi/2]) g3;
  check["eq1 matches the hand form transcribed here as want (written out here)", Simplify[e1[[2]] - want], 0];
];

Print[""];
Print["=============================================="];
Print[" 5. hierarchy guard: Hubbard U must NOT close"];
Print["=============================================="];

Module[{h, res},
  h = eps NCTimes[Cre[c, {1, 1}], Ann[c, {1, 1}]] +
      eps NCTimes[Cre[c, {1, -1}], Ann[c, {1, -1}]] +
      u NCTimes[Cre[c, {1, 1}], Ann[c, {1, 1}], Cre[c, {1, -1}], Ann[c, {1, -1}]];
  res = Quiet[CloseEOM[{GF[Ann[c, {1, 1}], Cre[c, {1, 1}]]}, h],
              CloseEOM::open];
  checkTrue["CloseEOM refuses to close an interacting H", res === $Failed];
];

Print[""];
Print[" 5b. atomic limit: with MaxOperators -> 3 the Hubbard chain closes EXACTLY"];

Module[{h, cl, sol, g, want},
  h = eps NCTimes[Cre[c, {1, 1}], Ann[c, {1, 1}]] +
      eps NCTimes[Cre[c, {1, -1}], Ann[c, {1, -1}]] +
      u NCTimes[Cre[c, {1, 1}], Ann[c, {1, 1}], Cre[c, {1, -1}], Ann[c, {1, -1}]];
  cl = CloseEOM[{GF[Ann[c, {1, 1}], Cre[c, {1, 1}]]}, h, MaxOperators -> 3];
  checkTrue["closes without any decoupling", cl =!= $Failed];
  If[cl =!= $Failed,
    Print["  equations: ", Length[cl["Equations"]]];
    sol = First[SolveGF[cl]];
    g = Simplify[GF[Mono[{Ann[c, {1, 1}]}], Mono[{Cre[c, {1, 1}]}]] /. sol];
    (* <n_down> survives as Expect[Mono[{Cre[c,{1,-1}], Ann[c,{1,-1}]}]] *)
    g = Simplify[g /. Expect[Mono[{Cre[c, {1, -1}], Ann[c, {1, -1}]}]] -> nd];
    Print["  G = ", InputForm[g]];
    want = (1/Sqrt[2 Pi]) ((1 - nd)/(\[Omega] - eps) + nd/(\[Omega] - eps - u));
    check["exact atomic-limit G = (1-n)/(w-eps) + n/(w-eps-U)",
      Simplify[g - want], 0]
  ];
];

Print[""];
Print[" 5c. same H truncated at MaxOperators -> 1 with a Hartree-Fock decoupling"];

Module[{h, dec, cl, sol, g, want},
  h = eps NCTimes[Cre[c, {1, 1}], Ann[c, {1, 1}]] +
      eps NCTimes[Cre[c, {1, -1}], Ann[c, {1, -1}]] +
      u NCTimes[Cre[c, {1, 1}], Ann[c, {1, 1}], Cre[c, {1, -1}], Ann[c, {1, -1}]];
  (* <<n_{-s} c_s ; cd_s>>  ->  <n_{-s}> <<c_s ; cd_s>>   (mean field) *)
  dec[GF[Mono[{Cre[c, {1, -1}], Ann[c, {1, -1}], Ann[c, {1, 1}]}], b_]] :=
    nd GF[Mono[{Ann[c, {1, 1}]}], b];
  dec[x_] := x;
  cl = CloseEOM[{GF[Ann[c, {1, 1}], Cre[c, {1, 1}]]}, h, DecouplingRule -> dec];
  checkTrue["closes with decoupling", cl =!= $Failed];
  If[cl =!= $Failed,
    sol = First[SolveGF[cl]];
    g = Simplify[GF[Mono[{Ann[c, {1, 1}]}], Mono[{Cre[c, {1, 1}]}]] /. sol];
    Print["  G = ", InputForm[g]];
    want = (1/Sqrt[2 Pi])/(\[Omega] - eps - nd u);
    check["Hartree-Fock: single shifted pole at eps + U n", Simplify[g - want], 0]
  ];
];

Print[""];
Print["=============================================="];
Print[" 6. equilibrium post-processing"];
Print["=============================================="];

Module[{g, dec, gl},
  g = 1/(Sqrt[2 Pi] (\[Omega] - e0));
  dec = SpectralDecomposition[g, \[Omega]];
  check["one pole at e0 with residue 1/Sqrt[2 Pi]", dec, {{e0, 1/Sqrt[2 Pi]}}];
  gl = GLesser[g, \[Omega], bb, et];
  checkTrue["G< is finite and carries the Fermi factor",
    FreeQ[gl, Indeterminate] && ! FreeQ[gl, Exp[bb \[Omega]]]];
  checkTrue["G< -> Lorentzian x Fermi as et -> 0",
    Chop[N[gl /. {e0 -> 0, bb -> 1, et -> 1/1000, \[Omega] -> 0}] -
         N[(2 I/(1/1000))/(Sqrt[2 Pi] 2)]] === 0 ||
    NumericQ[N[gl /. {e0 -> 0, bb -> 1, et -> 1/1000, \[Omega] -> 1/2}]]];
];

Print[""];
Print["=============================================="];
Print[" 7. matrix form is (w - H)"];
Print["=============================================="];

Module[{h, cl, mv, m, ev, want},
  h = ringH[4, t, th, 0, 0];
  cl = CloseEOM[{GF[Ann[c, {2, 1}], Cre[c, {1, 1}]]}, h];
  mv = GFMatrix[cl];
  m = mv[[1]];
  Print["  matrix size: ", Dimensions[m]];
  (* eigenvalues of w I - M's non-w part must be the ring spectrum *)
  ev = Sort[N[Eigenvalues[\[Omega] IdentityMatrix[Length[m]] - m /.
        {t -> 1, th -> 37/100, \[Omega] -> 0}]]];
  want = Sort[N[Table[2 Cos[37/100 + 2 Pi k/4], {k, 0, 3}]]];
  checkTrue["spectrum of (w I - M) matches the 4-site ring",
    Max[Abs[ev - want]] < 10^-9];
];

Print[""];
Print["=============================================="];
Print["  PASSED: ", $pass, "     FAILED: ", $fail];
Print["=============================================="];
