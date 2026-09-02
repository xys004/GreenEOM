(* Regression suite for GreenNEQ (steady-state Keldysh). Pure ASCII. *)

$pass = 0; $fail = 0;
check[name_String, got_, want_] :=
  If[TrueQ[Simplify[got - want == 0]] || TrueQ[got === want],
    $pass++; Print["  PASS  ", name],
    $fail++; Print["  FAIL  ", name]; Print["        got  = ", InputForm[got]];
    Print["        want = ", InputForm[want]]];
checkTrue[name_String, v_] :=
  If[TrueQ[v], $pass++; Print["  PASS  ", name],
    $fail++; Print["  FAIL  ", name, "  -> ", InputForm[v]]];
checkNum[name_String, got_, want_, tol_: 10^-9] :=
  checkTrue[name, Abs[N[got] - N[want]] < tol];

DeclareSpecies[c, "Fermion"];
DeclareSpecies[dd, "Fermion"];

Print["=============================================="];
Print[" 1. single-particle matrix from the commutator"];
Print["=============================================="];

Module[{h, m},
  h = Sum[e[n] NCTimes[Cre[c, n], Ann[c, n]], {n, 3}] +
      Sum[tt (NCTimes[Cre[c, n + 1], Ann[c, n]] +
              NCTimes[Cre[c, n], Ann[c, n + 1]]), {n, 2}];
  m = HamiltonianMatrix[h, Table[Ann[c, n], {n, 3}]];
  check["open chain gives a tridiagonal matrix", m,
    {{e[1], tt, 0}, {tt, e[2], tt}, {0, tt, e[3]}}];
  checkTrue["matrix is Hermitian", m === Transpose[m]];
];

Module[{h, m},
  h = u NCTimes[Cre[c, 1], Ann[c, 1], Cre[c, 2], Ann[c, 2]];
  m = Quiet[HamiltonianMatrix[h, {Ann[c, 1], Ann[c, 2]}],
        HamiltonianMatrix::notquad];
  checkTrue["refuses a quartic H", m === $Failed];
];

Print[""];
Print["=============================================="];
Print[" 2. the self-energy IS the Schur complement"];
Print["=============================================="];

Module[{h, m, dfl, sigma, want},
  (* device level 0 hybridised with two discrete bath levels *)
  h = e0 NCTimes[Cre[dd, 0], Ann[dd, 0]] +
      Sum[ek[k] NCTimes[Cre[c, k], Ann[c, k]], {k, 2}] +
      Sum[v[k] (NCTimes[Cre[dd, 0], Ann[c, k]] +
                NCTimes[Cre[c, k], Ann[dd, 0]]), {k, 2}];
  m = HamiltonianMatrix[h, {Ann[dd, 0], Ann[c, 1], Ann[c, 2]}];
  dfl = Downfold[m, {1}, \[Omega]];
  sigma = Simplify[dfl[[2]][[1, 1]]];
  want = Sum[v[k]^2/(\[Omega] - ek[k]), {k, 2}];
  check["Sigma = Sum_k V_k^2/(w - ek), same as the EOM route",
    Simplify[sigma - want], 0];
];

Print[""];
Print["=============================================="];
Print[" 3. Langreth rules"];
Print["=============================================="];

check["(AB)^r = A^r B^r", ContourRetarded[CProd[aa, bb]], MDot[Rt[aa], Rt[bb]]];
check["(AB)^< = A^r B^< + A^< B^a",
  ContourLesser[CProd[aa, bb]], MDot[Rt[aa], Ls[bb]] + MDot[Ls[aa], Ad[bb]]];
check["(ABC)^< has the three standard terms",
  ContourLesser[CProd[aa, bb, cc]],
  MDot[Rt[aa], Rt[bb], Ls[cc]] + MDot[Rt[aa], Ls[bb], Ad[cc]] +
  MDot[Ls[aa], Ad[bb], Ad[cc]]];
checkTrue["(ABCD)^< has four terms",
  Length[List @@ ContourLesser[CProd[aa, bb, cc, ee]]] === 4];

Module[{keldysh},
  (* Dyson on the contour, G = g + g Sigma G, projected with Langreth.
     Dropping the boundary term g^< (steady state, no initial correlations)
     the surviving structure is exactly G^< = G^r Sigma^< G^a. *)
  keldysh = ContourLesser[CProd[gr, sig, gg]];
  Print["  (g Sigma G)^< = ", InputForm[keldysh]];
  checkTrue["Keldysh equation term G^r Sigma^< G^a appears",
    ! FreeQ[keldysh, MDot[Rt[gr], Ls[sig], Ad[gg]]]];
];

Print[""];
Print["=============================================="];
Print[" 4. resonant level: exact Breit-Wigner"];
Print["=============================================="];

Module[{h, gl, gr2, sig, gret, tt2, want},
  h = {{e0}};
  gl = {{gamL}}; gr2 = {{gamR}};
  sig = WideBandSelfEnergy[gl + gr2];
  gret = RetardedG[h, sig, \[Omega]];
  tt2 = Simplify[Transmission[gret, gl, gr2]];
  Print["  T(w) = ", InputForm[tt2]];
  want = gamL gamR/((\[Omega] - e0)^2 + (gamL + gamR)^2/4);
  check["T(w) = GamL GamR / ((w-e0)^2 + ((GamL+GamR)/2)^2)",
    Simplify[tt2 - want], 0];
  checkNum["T = 1 on resonance with symmetric coupling",
    tt2 /. {gamL -> 1/5, gamR -> 1/5, \[Omega] -> 0, e0 -> 0}, 1];
  checkTrue["T <= 1 for asymmetric coupling",
    And @@ Table[N[tt2 /. {gamL -> 3/10, gamR -> 1/10, e0 -> 0,
        \[Omega] -> wv}] <= 1 + 10^-12, {wv, -2, 2, 1/4}]];
  checkNum["spectral sum rule  Int A dw = 1",
    NIntegrate[
      Re[-2 Im[gret[[1, 1]]] /. {e0 -> 0, gamL -> 1/5, gamR -> 3/10}] /(2 Pi),
      {\[Omega], -Infinity, Infinity}], 1, 10^-6];
];

Print[""];
Print["=============================================="];
Print[" 5. equilibrium must be RECOVERED, not imposed"];
Print["=============================================="];

Module[{h, gl, grr, sig, gret, gadv, slss, glss, fdt},
  h = {{e0}};
  gl = {{gamL}}; grr = {{gamR}};
  sig = WideBandSelfEnergy[gl + grr];
  gret = RetardedG[h, sig, \[Omega]];
  gadv = AdvancedG[gret];
  (* same chemical potential in both leads = equilibrium *)
  slss = LesserSigma[{gl, grr}, {ff, ff}];
  glss = LesserG[gret, slss];
  fdt = Simplify[glss + ff (gret - gadv)];
  check["G^< = -(G^r - G^a) f  emerges at equal mu", Simplify[fdt[[1, 1]]], 0];
  Print["  -> the fluctuation-dissipation relation is an OUTPUT here,"];
  Print["     not an assumption: it is what G^r Sigma^< G^a collapses to."];
];

Print[""];
Print["=============================================="];
Print[" 6. Meir-Wingreen: conservation and Landauer"];
Print["=============================================="];

Module[{h, gl, grr, sig, gret, gadv, slss, glss, il, ir, tr, subs},
  h = {{e0, tc}, {tc, e1}};                        (* double dot *)
  gl = {{gamL, 0}, {0, 0}}; grr = {{0, 0}, {0, gamR}};
  sig = WideBandSelfEnergy[gl + grr];
  gret = RetardedG[h, sig, \[Omega]];
  gadv = AdvancedG[gret];
  slss = LesserSigma[{gl, grr}, {fL, fR}];
  glss = LesserG[gret, slss];
  il = MeirWingreenCurrent[gret, glss, gl, fL];
  ir = MeirWingreenCurrent[gret, glss, grr, fR];
  check["current conservation  I_L + I_R = 0 pointwise in w",
    Simplify[il + ir], 0];
  tr = Transmission[gret, gl, grr];
  check["Meir-Wingreen reduces to Landauer for a quadratic H",
    Simplify[I il - tr (fL - fR)/(2 Pi)], 0];
  subs = {e0 -> 0, e1 -> 0, tc -> 1/2, gamL -> 1/5, gamR -> 1/5};
  checkNum["zero bias gives zero current",
    I il /. subs /. {fL -> ff, fR -> ff, \[Omega] -> 1/3}, 0];
  Print["  T(w) double dot = ", InputForm[Simplify[tr /. {gamL -> g0, gamR -> g0}]]];
];

Print[""];
Print["=============================================="];
Print[" 7. Aharonov-Bohm ring wired to two leads"];
Print["=============================================="];

ringHmat[ns_, t_, th_] :=
  Table[Which[
     Mod[j - i, ns] == 1, t Exp[I th],
     Mod[i - j, ns] == 1, t Exp[-I th],
     True, 0], {i, ns}, {j, ns}];

Module[{ns = 6, hm, gl, grr, sig, gret, tfun, t1, t2, t3, per},
  hm = ringHmat[ns, tt, th];
  gl = SparseArray[{{1, 1} -> gam}, {ns, ns}] // Normal;
  grr = SparseArray[{{4, 4} -> gam}, {ns, ns}] // Normal;
  sig = WideBandSelfEnergy[gl + grr];
  tfun[thv_, wv_] :=
    Module[{s, g},
      s = sig /. {gam -> 3/10};
      g = RetardedG[hm /. {tt -> 1, th -> thv}, s, wv];
      Re[N[Transmission[g, gl /. gam -> 3/10, grr /. gam -> 3/10]]]];
  t1 = tfun[2/5, 1/2];
  t2 = tfun[-2/5, 1/2];
  Print["  T(th=+0.4) = ", t1];
  Print["  T(th=-0.4) = ", t2];
  checkNum["two-terminal reciprocity: T(phi) = T(-phi)", t1, t2, 10^-9];
  per = tfun[2/5 + 2 Pi/ns, 1/2];
  checkNum["periodic in one flux quantum: th -> th + 2 Pi/N", t1, per, 10^-8];
  checkTrue["transmission actually depends on the flux (AB oscillation)",
    Abs[tfun[0, 1/2] - tfun[Pi/6, 1/2]] > 10^-3];
  checkTrue["0 <= T <= number of channels",
    And @@ Table[0 <= tfun[thv, 1/2] <= 1 + 10^-9, {thv, 0, Pi/3, Pi/30}]];
];

Print[""];
Print["=============================================="];
Print[" 8. occupation out of equilibrium"];
Print["=============================================="];

Module[{h, gl, grr, sig, gret, slss, glss, occ, nL, nR},
  h = {{0}};
  gl = {{gam}}; grr = {{gam}};
  sig = WideBandSelfEnergy[gl + grr];
  gret = RetardedG[h, sig, \[Omega]];
  slss = LesserSigma[{gl, grr}, {FermiF[\[Omega], muL, bt], FermiF[\[Omega], muR, bt]}];
  glss = LesserG[gret, slss];
  (* The Lorentzian tail falls off only as Gamma/(2 Pi w^2), so a window of
     +-20 loses Gamma/(2 Pi 20) = 3e-3 of the weight -- enough to break the
     sum rules. Integrate far out and hand NIntegrate the resonance. *)
  occ[muLv_, muRv_] :=
    NIntegrate[
      Re[N[Occupation[glss, 1] /. {gam -> 1/5, bt -> 20, muL -> muLv, muR -> muRv}]],
      {\[Omega], -10^5, -10, -1, 0, 1, 10, 10^5}];
  nL = occ[1, 1]; nR = occ[-1, -1];
  Print["  n at mu = +1 : ", nL];
  Print["  n at mu = -1 : ", nR];
  checkTrue["level fills as both reservoirs are raised", nL > 0.8 && nR < 0.2];
  checkNum["symmetric bias gives half filling", occ[1, -1], 1/2, 10^-4];
  checkNum["particle-hole symmetry: n(mu) + n(-mu) = 1", nL + nR, 1, 10^-4];
];

Print[""];
Print["=============================================="];
Print["  PASSED: ", $pass, "     FAILED: ", $fail];
Print["=============================================="];
