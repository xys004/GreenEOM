(* Worked examples for GreenEOM. Doubles as the user manual. Pure ASCII. *)

$pass = 0; $fail = 0;
check[name_String, got_, want_] :=
  If[TrueQ[Simplify[got - want == 0]] || TrueQ[got === want],
    $pass++; Print["  PASS  ", name],
    $fail++; Print["  FAIL  ", name];
    Print["        got  = ", InputForm[got]]];
checkTrue[name_String, v_] :=
  If[TrueQ[v], $pass++; Print["  PASS  ", name], $fail++; Print["  FAIL  ", name]];

DeclareSpecies[c, "Fermion"];
DeclareSpecies[d, "Fermion"];
DeclareSpecies[b, "Boson"];

Print["##############################################"];
Print["# A. Resonant level: a self-energy, for free  "];
Print["##############################################"];
Print[""];
Print["  H = e0 dd d + sum_k ek cdk ck + sum_k V (dd ck + cdk d)"];
Print[""];

Module[{h, cl, sol, g, sigma, want},
  h = e0 NCTimes[Cre[d, 0], Ann[d, 0]] +
      Sum[ek[k] NCTimes[Cre[c, k], Ann[c, k]], {k, 2}] +
      Sum[v[k] (NCTimes[Cre[d, 0], Ann[c, k]] + NCTimes[Cre[c, k], Ann[d, 0]]), {k, 2}];
  cl = CloseEOM[{GF[Ann[d, 0], Cre[d, 0]]}, h];
  Print["  equations generated: ", Length[cl["Equations"]]];
  sol = First[SolveGF[cl]];
  g = Simplify[GF[Mono[{Ann[d, 0]}], Mono[{Cre[d, 0]}]] /. sol];
  Print["  G_dd = ", InputForm[g]];
  sigma = Sum[v[k]^2/(\[Omega] - ek[k]), {k, 2}];
  want = 1/(Sqrt[2 Pi] (\[Omega] - e0 - sigma));
  check["G_dd = norm / (w - e0 - Sum_k V_k^2/(w - ek))", Simplify[g - want], 0];
  Print["  -> the hybridisation self-energy was never typed in; it is what"];
  Print["     eliminating the bath Green functions produces."];
];

Print[""];
Print["##############################################"];
Print["# B. Mesoscopic ring (thesis Hamiltonian)     "];
Print["##############################################"];

(* Reusable ring builder.  Spin labels: +1 = up, -1 = down.
     H = sum_n ( cd_{n+1} . T(n) . c_n e^{I th} + h.c. )
     T(n) = {{t,                             -I al e^{-I ph} + be e^{I ph}},
             {-I al e^{I ph} - be e^{-I ph},  t                          }}
   with ph = ph_{n,n+1} = (phi_n + phi_{n+1})/2,  phi_n = 2 Pi (n-1)/N. *)

ringH[ns_, t_, th_, al_, be_] :=
  Module[{phi, ph, nxt, tm, tmc, spins = {1, -1}},
    phi[n_] := 2 Pi (n - 1)/ns;
    (* midpoint of the arc the hop spans, as an offset from phi_n:
       averaging the two site angles puts the closing bond, where
       phi_1 reads as 0 rather than 2 Pi, a full Pi away from the rest *)
    ph[n_, m_] := phi[n] + Mod[m - n, ns] Pi/ns;
    nxt[n_] := Mod[n, ns] + 1;
    tm[n_][1, 1] = t; tm[n_][-1, -1] = t;
    tm[n_][1, -1] := -I al Exp[-I ph[n, nxt[n]]] + be Exp[I ph[n, nxt[n]]];
    tm[n_][-1, 1] := -I al Exp[I ph[n, nxt[n]]] - be Exp[-I ph[n, nxt[n]]];
    tmc[n_][1, 1] = t; tmc[n_][-1, -1] = t;
    tmc[n_][1, -1] := I al Exp[I ph[n, nxt[n]]] + be Exp[-I ph[n, nxt[n]]];
    tmc[n_][-1, 1] := I al Exp[-I ph[n, nxt[n]]] - be Exp[I ph[n, nxt[n]]];
    Sum[
      Exp[I th] tm[n][s, sp] NCTimes[Cre[c, {nxt[n], s}], Ann[c, {n, sp}]] +
      Exp[-I th] tmc[n][s, sp] NCTimes[Cre[c, {n, sp}], Ann[c, {nxt[n], s}]],
      {n, ns}, {s, spins}, {sp, spins}]
  ];

Print[""];
Print[" B1. ten sites, hopping only -- against the thesis result (eq. green1)"];

Module[{h, cl, sol, g1, g2, w1, w2, den},
  h = ringH[10, t, th, 0, 0];
  cl = CloseEOM[{GF[Ann[c, {2, 1}], Cre[c, {1, 1}]],
                 GF[Ann[c, {1, 1}], Cre[c, {2, 1}]]}, h];
  Print["  equations: ", Length[cl["Equations"]],
        "   unknowns: ", Length[cl["Unknowns"]]];
  sol = First[SolveGF[cl]];
  g1 = FullSimplify[GF[Mono[{Ann[c, {2, 1}]}], Mono[{Cre[c, {1, 1}]}]] /. sol];
  g2 = FullSimplify[GF[Mono[{Ann[c, {1, 1}]}], Mono[{Cre[c, {2, 1}]}]] /. sol];

  den = Sqrt[2 Pi] (2 t^10 Cos[10 th] +
        (2 t^2 - \[Omega]^2) (t^8 - 12 t^6 \[Omega]^2 + 19 t^4 \[Omega]^4 -
                              8 t^2 \[Omega]^6 + \[Omega]^8));
  w1 = Exp[-9 I th] (-t^9 - Exp[10 I th] t (t^8 - 10 t^6 \[Omega]^2 +
        15 t^4 \[Omega]^4 - 7 t^2 \[Omega]^6 + \[Omega]^8))/den;
  w2 = -Exp[-I th] t ((1 + Exp[10 I th]) t^8 - 10 t^6 \[Omega]^2 +
        15 t^4 \[Omega]^4 - 7 t^2 \[Omega]^6 + \[Omega]^8)/den;

  check["G_{n+1,n} matches thesis eq. (green1)", FullSimplify[g1 - w1], 0];
  check["G_{n,n+1} matches thesis eq. (610)", FullSimplify[g2 - w2], 0];
];

Print[""];
Print[" B2. spin-orbit switched on: how the system size grows"];

Module[{sizes},
  sizes = Table[
    Module[{cl},
      cl = CloseEOM[{GF[Ann[c, {2, 1}], Cre[c, {1, 1}]]}, ringH[n, t, th, al, be]];
      {n, Length[cl["Unknowns"]]}], {n, {2, 3, 4, 6}}];
  Print["  {sites, unknowns} = ", InputForm[sizes]];
  checkTrue["with Rashba/Dresselhaus the system is 2N, not N",
    And @@ (#[[2]] === 2 #[[1]] & /@ sizes)];
];

Print[""];
Print["##############################################"];
Print["# C. From G to the persistent charge current  "];
Print["##############################################"];
Print[""];
Print["  J_c(w) = (1/N) sum_n ( e^{-I th} t G^<_{n+1,n} - e^{I th} t G^<_{n,n+1} )"];
Print[""];

Module[{ns = 6, h, cl, sol, gl, jc, num, poles},
  h = ringH[ns, t, th, 0, 0];
  cl = CloseEOM[
     Flatten[Table[{GF[Ann[c, {Mod[n, ns] + 1, 1}], Cre[c, {n, 1}]],
                    GF[Ann[c, {n, 1}], Cre[c, {Mod[n, ns] + 1, 1}]]}, {n, ns}]], h];
  Print["  equations: ", Length[cl["Equations"]]];
  sol = First[SolveGF[cl]];
  gl[aa_, bb_] := GLesser[
     GF[Mono[{Ann[c, aa]}], Mono[{Cre[c, bb]}]] /. sol, \[Omega], beta, eta];
  jc = (1/ns) Sum[
     Exp[-I th] t gl[{Mod[n, ns] + 1, 1}, {n, 1}] -
     Exp[I th] t gl[{n, 1}, {Mod[n, ns] + 1, 1}], {n, ns}];
  (* numeric sanity: the current must be real, odd in the flux, and vanish at th=0 *)
  num[thv_, wv_] := Chop[N[jc /. {t -> 1, th -> thv, beta -> 5, eta -> 1/50,
                                  \[Omega] -> wv}], 10^-12];
  Print["  J_c(th=0.0,  w=0.5) = ", InputForm[num[0, 1/2]]];
  Print["  J_c(th=0.4,  w=0.5) = ", InputForm[num[2/5, 1/2]]];
  Print["  J_c(th=-0.4, w=0.5) = ", InputForm[num[-2/5, 1/2]]];
  checkTrue["J_c vanishes at zero flux", Abs[num[0, 1/2]] < 10^-10];
  checkTrue["J_c is odd in the flux",
    Abs[num[2/5, 1/2] + num[-2/5, 1/2]] < 10^-8];
  checkTrue["J_c is real", Abs[Im[num[2/5, 1/2]]] < 10^-8];
];

Print[""];
Print["##############################################"];
Print["# D. Bosons: two coupled modes                "];
Print["##############################################"];

Module[{h, cl, sol, g, want},
  h = w1 NCTimes[Cre[b, 1], Ann[b, 1]] + w2 NCTimes[Cre[b, 2], Ann[b, 2]] +
      jj (NCTimes[Cre[b, 1], Ann[b, 2]] + NCTimes[Cre[b, 2], Ann[b, 1]]);
  cl = CloseEOM[{GF[Ann[b, 1], Cre[b, 1]]}, h];
  sol = First[SolveGF[cl]];
  g = Simplify[GF[Mono[{Ann[b, 1]}], Mono[{Cre[b, 1]}]] /. sol];
  Print["  G_11 = ", InputForm[g]];
  want = (\[Omega] - w2)/(Sqrt[2 Pi] ((\[Omega] - w1) (\[Omega] - w2) - jj^2));
  check["two-mode boson propagator", Simplify[g - want], 0];
  checkTrue["poles are the normal modes",
    Simplify[(\[Omega] /. Solve[Denominator[Together[g]] == 0, \[Omega]])] =!= {}];
];

Print[""];
Print["##############################################"];
Print["  PASSED: ", $pass, "     FAILED: ", $fail];
Print["##############################################"];
