(* Does  P_z(alpha,beta) = -P_z(beta,alpha)  survive impurities?

   Stated prediction, before any number is computed:

     The identity comes from U = spin rotation by pi about [110], which swaps
     sigma_x <-> sigma_y and flips sigma_z. U acts on spin space only and is
     site independent, so it commutes with ANY term proportional to the
     identity in spin.

       (1) charge disorder  eps_n * 1   is such a term
           => the identity survives exactly, however strong and however random,
              and alpha = beta still kills P_z exactly.
       (2) a term touching sigma_z (Zeeman impurity eps_n*sigma_z, or
           spin-polarised contacts) is NOT invariant: U flips its sign
           => the identity breaks, deforming into
              P_z(alpha,beta,h) = -P_z(beta,alpha,-h),
              which at alpha = beta no longer forces zero.

   Everything below is an attempt to falsify that. Pure ASCII. *)

$pass = 0; $fail = 0;
checkTrue[name_String, v_] :=
  If[TrueQ[v], $pass++; Print["  PASS  ", name],
    $fail++; Print["  FAIL  ", name, "  -> ",
      StringTake[ToString[InputForm[v]], UpTo[140]]]];

DeclareSpecies[c, "Fermion"];

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

idx[n_, s_] := 2 (n - 1) + If[s === 1, 1, 2];
ringOps[ns_] := Flatten[Table[Ann[c, {n, s}], {n, ns}, {s, {1, -1}}]];

hsym[ns_] := hsym[ns] = HamiltonianMatrix[ringH[ns, 1, thS, alS, beS], ringOps[ns]];
hclean[ns_, alv_, bev_, thv_] :=
  N[hsym[ns] /. {thS -> thv, alS -> alv, beS -> bev}];

(* impurities are diagonal: eps[[n]] is the charge part, zee[[n]] the sigma_z part *)
hwith[ns_, alv_, bev_, thv_, eps_List, zee_List] :=
  hclean[ns, alv, bev, thv] +
  DiagonalMatrix[N[Flatten[Table[{eps[[n]] + zee[[n]], eps[[n]] - zee[[n]]},
     {n, ns}]]]];

$eta = 10^-8;

(* leads on sites nl and nr; gu/gd allow spin-polarised contacts *)
projm[ns_, n_, s_] := projm[ns, n, s] =
  Normal[SparseArray[{{idx[n, s], idx[n, s]} -> 1}, {2 ns, 2 ns}]];

polz[h_, ns_, nl_, nr_, gu_, gd_, w_] :=
  Module[{gl, gr2, sig, gr, ga, m},
    gl = gu projm[ns, nl, 1] + gd projm[ns, nl, -1];
    gr2 = gu projm[ns, nr, 1] + gd projm[ns, nr, -1];
    sig = WideBandSelfEnergy[gl + gr2];
    gr = RetardedG[h, sig, N[w] + I $eta];
    ga = Dagger[gr];
    m = Table[
      Re[Tr[(gu^(Boole[sout === 1]) gd^(Boole[sout === -1]) projm[ns, nr, sout]) .
            gr . (gu^(Boole[sin === 1]) gd^(Boole[sin === -1]) projm[ns, nl, sin]) . ga]],
      {sout, {1, -1}}, {sin, {1, -1}}];
    (m[[1, 1]] + m[[1, 2]] - m[[2, 1]] - m[[2, 2]])/Total[m, 2]
  ];

zeros[ns_] := ConstantArray[0, ns];

Print["=============================================="];
Print[" 0. the diagonal shortcut is the real thing"];
Print["=============================================="];

Module[{ns = 6, eps, hop, hdiag, hfull},
  eps = {3/10, -1/5, 0, 2/5, 1/10, -7/10};
  hop = ringH[6, 1, 1/5, 3/10, 1/10] +
        Sum[eps[[n]] NCTimes[Cre[c, {n, s}], Ann[c, {n, s}]], {n, 6}, {s, {1, -1}}];
  hfull = N[HamiltonianMatrix[hop, ringOps[6]]];
  hdiag = hwith[6, 3/10, 1/10, 1/5, eps, zeros[6]];
  checkTrue["adding a diagonal == putting eps_n in the operator H",
    Max[Abs[hfull - hdiag]] < 10^-12];
];

Print[""];
Print["=============================================="];
Print[" 1. control: clean ring (the identity as found)"];
Print["=============================================="];

Module[{ns = 6, worst},
  worst = Max[Flatten[Table[
     Abs[polz[hwith[ns, av, bv, 2 Pi p/ns, zeros[ns], zeros[ns]], ns, 1, 4, 3/10, 3/10, wv] +
         polz[hwith[ns, bv, av, 2 Pi p/ns, zeros[ns], zeros[ns]], ns, 1, 4, 3/10, 3/10, wv]],
     {av, {1/5, 1/2}}, {bv, {0, 7/10}}, {p, {1/8, 1/3}}, {wv, {-4/5, 3/5}}]]];
  Print["  max |P_z(a,b) + P_z(b,a)|, clean : ", worst];
  checkTrue["control reproduces the identity", worst < 10^-12];
];

Print[""];
Print["=============================================="];
Print[" 2. PREDICTION 1: charge disorder cannot break it"];
Print["=============================================="];

Module[{ns = 6, eps1, worst1, worst2, worst3, worstBig, worstEq, reals},
  (* 2a. one strong impurity *)
  eps1 = {9/10, 0, 0, 0, 0, 0};
  worst1 = Max[Flatten[Table[
     Abs[polz[hwith[ns, av, bv, 2 Pi p/ns, eps1, zeros[ns]], ns, 1, 4, 3/10, 3/10, wv] +
         polz[hwith[ns, bv, av, 2 Pi p/ns, eps1, zeros[ns]], ns, 1, 4, 3/10, 3/10, wv]],
     {av, {1/5, 1/2}}, {bv, {0, 7/10}}, {p, {1/8, 1/3}}, {wv, {-4/5, 3/5}}]]];
  Print["  2a. one impurity eps_1 = 0.9      : ", worst1];
  checkTrue["single charge impurity: identity survives", worst1 < 10^-12];

  (* 2b. full random disorder, several realisations *)
  SeedRandom[20260816];
  reals = Table[RandomReal[{-1, 1}, ns], {8}];
  worst2 = Max[Flatten[Table[
     Abs[polz[hwith[ns, av, bv, 2 Pi p/ns, r, zeros[ns]], ns, 1, 4, 3/10, 3/10, wv] +
         polz[hwith[ns, bv, av, 2 Pi p/ns, r, zeros[ns]], ns, 1, 4, 3/10, 3/10, wv]],
     {r, reals}, {av, {1/5, 1/2}}, {bv, {0, 7/10}}, {p, {1/8, 1/3}}, {wv, {-4/5, 3/5}}]]];
  Print["  2b. 8 random realisations, W = 1  : ", worst2];
  checkTrue["random charge disorder: identity survives", worst2 < 10^-12];

  (* 2c. very strong disorder, well past the clean bandwidth *)
  SeedRandom[7];
  worstBig = Max[Flatten[Table[
     Abs[polz[hwith[ns, av, bv, 2 Pi p/ns, r, zeros[ns]], ns, 1, 4, 3/10, 3/10, wv] +
         polz[hwith[ns, bv, av, 2 Pi p/ns, r, zeros[ns]], ns, 1, 4, 3/10, 3/10, wv]],
     {r, Table[RandomReal[{-5, 5}, ns], {5}]},
     {av, {3/10, 9/10}}, {bv, {0, 1/2}}, {p, {1/4}}, {wv, {-3/2, 1/2}}]]];
  Print["  2c. strong disorder, W = 10       : ", worstBig];
  checkTrue["strong charge disorder: identity survives", worstBig < 10^-12];

  (* 2d. break the lead symmetry too *)
  SeedRandom[99];
  worst3 = Max[Flatten[Table[
     Abs[polz[hwith[ns, av, bv, 2 Pi p/ns, r, zeros[ns]], ns, 1, 3, 3/10, 3/10, wv] +
         polz[hwith[ns, bv, av, 2 Pi p/ns, r, zeros[ns]], ns, 1, 3, 3/10, 3/10, wv]],
     {r, Table[RandomReal[{-1, 1}, ns], {4}]},
     {av, {1/5, 3/5}}, {bv, {0, 2/5}}, {p, {1/4}}, {wv, {-1/2, 4/5}}]]];
  Print["  2d. disorder + asymmetric leads   : ", worst3];
  checkTrue["no spatial symmetry left: identity STILL survives", worst3 < 10^-12];

  (* 2e. and therefore alpha = beta must still give exactly zero *)
  SeedRandom[555];
  worstEq = Max[Flatten[Table[
     Abs[polz[hwith[ns, av, av, 2 Pi p/ns, r, zeros[ns]], ns, 1, 3, 3/10, 3/10, wv]],
     {r, Table[RandomReal[{-2, 2}, ns], {6}]},
     {av, {1/5, 1/2, 1}}, {p, {1/8, 1/4}}, {wv, {-1, 1/2}}]]];
  Print["  2e. |P_z| at alpha = beta, dirty  : ", worstEq];
  checkTrue["alpha = beta still kills P_z exactly in a dirty ring",
    worstEq < 10^-12];
];

Print[""];
Print["=============================================="];
Print[" 3. PREDICTION 2: a sigma_z term must break it"];
Print["=============================================="];

Module[{ns = 6, zee, viol, defo, tab},
  zee = {3/10, 0, 0, 0, 0, 0};
  viol = Max[Flatten[Table[
     Abs[polz[hwith[ns, av, bv, 2 Pi p/ns, zeros[ns], zee], ns, 1, 4, 3/10, 3/10, wv] +
         polz[hwith[ns, bv, av, 2 Pi p/ns, zeros[ns], zee], ns, 1, 4, 3/10, 3/10, wv]],
     {av, {1/5, 1/2}}, {bv, {0, 7/10}}, {p, {1/8, 1/3}}, {wv, {-4/5, 3/5}}]]];
  Print["  3a. Zeeman impurity h_1 = 0.3, max violation : ", viol];
  checkTrue["a sigma_z impurity DOES break the identity", viol > 10^-3];

  (* the deformed identity should hold instead *)
  defo = Max[Flatten[Table[
     Abs[polz[hwith[ns, av, bv, 2 Pi p/ns, zeros[ns], zee], ns, 1, 4, 3/10, 3/10, wv] +
         polz[hwith[ns, bv, av, 2 Pi p/ns, zeros[ns], -zee], ns, 1, 4, 3/10, 3/10, wv]],
     {av, {1/5, 1/2}}, {bv, {0, 7/10}}, {p, {1/8, 1/3}}, {wv, {-4/5, 3/5}}]]];
  Print["  3b. deformed identity P(a,b,h) = -P(b,a,-h)  : ", defo];
  checkTrue["the identity is deformed, not destroyed", defo < 10^-12];

  Print[""];
  Print["  3c. how P_z at alpha = beta grows with the Zeeman strength"];
  Print["      (alpha = beta = 0.3, phi = 0.25, w = 0.5, one site)"];
  tab = Table[{N[h],
     polz[hwith[ns, 3/10, 3/10, 2 Pi (1/4)/ns, zeros[ns],
        ReplacePart[zeros[ns], 1 -> h]], ns, 1, 4, 3/10, 3/10, 1/2]},
     {h, {0, 1/100, 1/20, 1/10, 3/10, 1/2, 1}}];
  Do[Print["      h = ", PaddedForm[N[r[[1]]], {4, 2}],
           "   P_z = ", PaddedForm[r[[2]], {8, 5}]], {r, tab}];
  checkTrue["P_z at alpha=beta grows monotonically from zero with h",
    tab[[1, 2]] < 10^-12 && Abs[tab[[-1, 2]]] > 10^-3];
];

Print[""];
Print["=============================================="];
Print[" 4. spin-polarised contacts also break it"];
Print["=============================================="];

Module[{ns = 6, viol, viol2},
  viol = Max[Flatten[Table[
     Abs[polz[hwith[ns, av, bv, 2 Pi p/ns, zeros[ns], zeros[ns]], ns, 1, 4, 4/10, 2/10, wv] +
         polz[hwith[ns, bv, av, 2 Pi p/ns, zeros[ns], zeros[ns]], ns, 1, 4, 4/10, 2/10, wv]],
     {av, {1/5, 1/2}}, {bv, {0, 7/10}}, {p, {1/8, 1/3}}, {wv, {-4/5, 3/5}}]]];
  Print["  Gamma_up = 0.4, Gamma_down = 0.2, max violation : ", viol];
  checkTrue["spin-polarised leads break the identity", viol > 10^-3];
  (* control: equal Gammas of a different value must not break it *)
  viol2 = Max[Flatten[Table[
     Abs[polz[hwith[ns, av, bv, 2 Pi p/ns, zeros[ns], zeros[ns]], ns, 1, 4, 3/5, 3/5, wv] +
         polz[hwith[ns, bv, av, 2 Pi p/ns, zeros[ns], zeros[ns]], ns, 1, 4, 3/5, 3/5, wv]],
     {av, {1/5, 1/2}}, {bv, {0, 7/10}}, {p, {1/8, 1/3}}, {wv, {-4/5, 3/5}}]]];
  Print["  control, Gamma_up = Gamma_down = 0.6            : ", viol2];
  checkTrue["it is the spin asymmetry that matters, not the coupling size",
    viol2 < 10^-12];
];

Print[""];
Print["=============================================="];
Print[" 5. does disorder change the SIZE of P_z?"];
Print["=============================================="];
Print["  The identity is untouched, but the polarisation itself is not."];
Print[""];

Module[{ns = 6, clean, dirty, ws, mk},
  mk[r_] := Max[Table[
     Abs[polz[hwith[ns, 3/10, 0, 2 Pi p/ns, r, zeros[ns]], ns, 1, 4, 3/10, 3/10, 1/2]],
     {p, Range[0, 50]/100}]];
  clean = mk[zeros[ns]];
  SeedRandom[2026];
  dirty = Table[mk[RandomReal[{-1, 1}, ns]], {6}];
  Print["  max |P_z| over flux, clean ring        : ", clean];
  Print["  same, 6 disorder realisations (W = 1)  : ", Sort[dirty]];
  Print["  median of the dirty realisations       : ", Median[dirty]];
  checkTrue["disorder shifts the magnitude while the identity holds",
    Max[Abs[dirty - clean]] > 10^-3];
];

Print[""];
Print["=============================================="];
Print[" 6. export: residual vs perturbation strength"];
Print["=============================================="];

Module[{ns = 6, grid, resid, strengths, d1, d2, hs},
  grid = Flatten[Table[{av, bv, p, wv},
     {av, {1/5, 1/2}}, {bv, {0, 7/10}}, {p, {1/8, 1/3}}, {wv, {-4/5, 3/5}}], 3];
  (* residual of the identity under a perturbation of amplitude s *)
  resid[eps_, zee_, gu_, gd_] := Max[Table[
     Abs[polz[hwith[ns, q[[1]], q[[2]], 2 Pi q[[3]]/ns, eps, zee], ns, 1, 4, gu, gd, q[[4]]] +
         polz[hwith[ns, q[[2]], q[[1]], 2 Pi q[[3]]/ns, eps, zee], ns, 1, 4, gu, gd, q[[4]]]],
     {q, grid}]];
  SeedRandom[31415];
  strengths = Join[{0}, Table[10^e, {e, -3, 1/2, 1/8}]];
  d1 = Table[
     Module[{r = RandomReal[{-1, 1}, ns]},
       {N[s],
        Max[10^-17, resid[s r, zeros[ns], 3/10, 3/10]],          (* charge   *)
        Max[10^-17, resid[zeros[ns], s r, 3/10, 3/10]],          (* Zeeman   *)
        (* lead polarisation p = s/(1+s) stays in [0,1), so Gamma_down > 0:
           a negative coupling is not a stronger perturbation, it is nonsense *)
        Max[10^-17, resid[zeros[ns], zeros[ns],
              3/10 (1 + s/(1 + s)), 3/10 (1 - s/(1 + s))]]}],
     {s, strengths}];
  Print["  s = 0      : ", Rest[d1[[1]]]];
  Print["  s = 0.1    : ", Rest[SelectFirst[d1, Abs[#[[1]] - 1/10] < 10^-3 &]]];
  Print["  s = 1.0    : ", Rest[SelectFirst[d1, Abs[#[[1]] - 1] < 10^-3 &]]];
  Export["/home/astrum/imp_residual.csv",
     Prepend[d1, {"s", "charge", "zeeman", "polarised_leads"}]];

  hs = Join[Range[0, 20]/100, Range[25, 100, 5]/100];
  d2 = Table[{N[h],
     polz[hwith[ns, 3/10, 3/10, 2 Pi (1/4)/ns, zeros[ns],
        ReplacePart[zeros[ns], 1 -> h]], ns, 1, 4, 3/10, 3/10, 1/2]}, {h, hs}];
  Export["/home/astrum/imp_pz.csv", Prepend[d2, {"h", "Pz_at_alpha_eq_beta"}]];
  Print["  exported 2 CSV files"];
];

Print[""];
Print["=============================================="];
Print["  PASSED: ", $pass, "     FAILED: ", $fail];
Print["=============================================="];
