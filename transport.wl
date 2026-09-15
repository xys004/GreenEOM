(* The Rashba-Dresselhaus mesoscopic ring, opened up: two leads,
   Landauer-Buttiker transport.
   The ring Hamiltonian is the SAME operator-level H that reproduced the
   reference persistent-current results, so the transport inherits that
   validation. Pure ASCII. *)

$pass = 0; $fail = 0;
checkTrue[name_String, v_] :=
  If[TrueQ[v], $pass++; Print["  PASS  ", name],
    $fail++; Print["  FAIL  ", name, "  -> ",
      StringTake[ToString[InputForm[v]], UpTo[160]]]];
checkNum[name_String, got_, want_, tol_: 10^-8] :=
  checkTrue[name <> "  (" <> ToString[N[got]] <> " vs " <> ToString[N[want]] <> ")",
    Abs[N[got] - N[want]] < tol];

DeclareSpecies[c, "Fermion"];

(* ---------------------------------------------------------------- *)
(* the ring, exactly as in the reference                              *)
(* ---------------------------------------------------------------- *)

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

(* basis order: (site 1 up, site 1 down, site 2 up, ...) *)
idx[n_, s_] := 2 (n - 1) + If[s === 1, 1, 2];
ringOps[ns_] := Flatten[Table[Ann[c, {n, s}], {n, ns}, {s, {1, -1}}]];

(* The operator-level step (40 commutators over a 160-term H) is the expensive
   one, so do it ONCE per ring size, symbolically in th/al/be, and substitute
   numbers afterwards. Rebuilding it per sweep point is what made the first
   version hang. *)
Clear[hsym, hmat, proj, gamL, gamR, tmat];
hsym[ns_] := hsym[ns] =
  HamiltonianMatrix[ringH[ns, 1, thS, alS, beS], ringOps[ns]];
hmat[ns_, alv_, bev_, thv_] :=
  N[hsym[ns] /. {thS -> thv, alS -> alv, beS -> bev}];

(* lead projectors: lead L on site 1, lead R diametrically opposite *)
proj[ns_, n_, s_] := proj[ns, n, s] =
  Normal[SparseArray[{{idx[n, s], idx[n, s]} -> 1}, {2 ns, 2 ns}]];
gamL[ns_, g_] := gamL[ns, g] = g (proj[ns, 1, 1] + proj[ns, 1, -1]);
gamR[ns_, g_] := gamR[ns, g] = With[{nr = 1 + Quotient[ns, 2]},
   g (proj[ns, nr, 1] + proj[ns, nr, -1])];

(* Leads sit on diametrically opposite sites, so every eigenstate that is odd
   about the lead axis has a node on BOTH contacts: it is decoupled from the
   two terminals -- a dark state. At those energies (w = a ring eigenvalue)
   the matrix w - H - Sigma is genuinely singular. The transmission itself is
   still well defined, because a state with psi(1) = 0 contributes nothing to
   G^r_{1,nR}; only the full inverse is ill conditioned. A small eta
   regularises it without touching the observable -- and test 1 checks that
   the answer does not depend on eta. *)
$eta = 10^-8;

tmat[ns_, alv_, bev_, thv_, g_, w_] := tmat[ns, alv, bev, thv, g, w] =
  tmatCompute[ns, alv, bev, thv, g, w];

tmatCompute[ns_, alv_, bev_, thv_, g_, w_] :=
  Module[{h, sig, gr, ga, gl, gr2, nr},
    h = hmat[ns, alv, bev, thv];
    gl = gamL[ns, g]; gr2 = gamR[ns, g];
    nr = 1 + Quotient[ns, 2];
    sig = WideBandSelfEnergy[gl + gr2];
    gr = RetardedG[h, sig, N[w] + I $eta];
    ga = Dagger[gr];
    Table[Re[Tr[(g proj[ns, nr, sout]) . gr . (g proj[ns, 1, sin]) . ga]],
      {sout, {1, -1}}, {sin, {1, -1}}]
  ];

tcharge[ns_, alv_, bev_, thv_, g_, w_] := Total[tmat[ns, alv, bev, thv, g, w], 2];
tflip[ns_, alv_, bev_, thv_, g_, w_] :=
  With[{m = tmat[ns, alv, bev, thv, g, w]}, m[[1, 2]] + m[[2, 1]]];
(* polarisation of the outgoing current for an unpolarised incoming beam *)
pol[ns_, alv_, bev_, thv_, g_, w_] :=
  With[{m = tmat[ns, alv, bev, thv, g, w]},
    (m[[1, 1]] + m[[1, 2]] - m[[2, 1]] - m[[2, 2]])/Total[m, 2]];

Print["=============================================="];
Print[" Rashba-Dresselhaus ring, opened: setup"];
Print["=============================================="];
Print["  6 sites, leads on sites 1 and 4, wide band, Gamma = 0.3, t = 1"];
Print["  th = 2 Pi phi / N  is the AB phase per bond; phi in flux quanta"];
Print[""];

th6[phi_] := 2 Pi phi/6;

Print["=============================================="];
Print[" 1. consistency of the open ring"];
Print["=============================================="];

Module[{m0, m1, gg, hh, sg, grr},
  (* the single-entry Gamma shortcut T = GamL GamR |G^r_{out,in}|^2 must agree
     with the general Tr[Gam_R G^r Gam_L G^a] -- checks the projector bookkeeping *)
  gg = 3/10; hh = hmat[6, 3/10, 2/10, 1/5];
  sg = WideBandSelfEnergy[gamL[6, gg] + gamR[6, gg]];
  grr = RetardedG[hh, sg, N[1/2] + I $eta];
  checkNum["projector trace = GamL GamR |G^r_{out,in}|^2",
    tmat[6, 3/10, 2/10, 1/5, gg, 1/2][[1, 1]],
    N[gg^2 Abs[grr[[idx[4, 1], idx[1, 1]]]]^2], 10^-10];
  (* the regularisation must not be doing any physics *)
  checkTrue["result independent of eta (dark-state regularisation is inert)",
    Module[{a, b},
      a = Total[Block[{$eta = 10^-8},
            tmatCompute[6, 3/10, 3/10, th6[1/4], 3/10, 9/10]], 2];
      b = Total[Block[{$eta = 10^-12},
            tmatCompute[6, 3/10, 3/10, th6[1/4], 3/10, 9/10]], 2];
      (* G ~ 1/(w - E + I eta), so T shifts by O(eta/(w-E)^2): the tolerance
         has to be looser than the regulariser, not tighter *)
      Abs[a - b] < 10^-5]];
  (* and it survives sitting exactly on a dark state (w = 1 is a ring eigenvalue) *)
  checkTrue["finite T exactly on a dark-state energy",
    Module[{v}, v = Total[tmatCompute[6, 0, 0, 0, 3/10, 1], 2];
      NumericQ[v] && 0 <= v <= 2 + 10^-9]];
  checkTrue["H matrix is Hermitian with SOC on",
    Max[Abs[hmat[6, 3/10, 3/10, 2/10] -
            ConjugateTranspose[hmat[6, 3/10, 3/10, 2/10]]]] < 10^-12];
  (* Buttiker two-terminal reciprocity: charge conductance even in the flux,
     and it must survive spin-orbit coupling *)
  checkNum["charge T even in flux, WITH Rashba+Dresselhaus",
    tcharge[6, 3/10, 3/10, th6[1/4], 3/10, 1/2],
    tcharge[6, 3/10, 3/10, th6[-1/4], 3/10, 1/2], 10^-9];
  checkNum["periodic in one flux quantum",
    tcharge[6, 3/10, 2/10, th6[1/4], 3/10, 1/2],
    tcharge[6, 3/10, 2/10, th6[1/4 + 1], 3/10, 1/2], 10^-8];
  checkTrue["unitarity: 0 <= T <= 2 open channels",
    And @@ Table[0 <= tcharge[6, 3/10, 3/10, th6[p], 3/10, 1/2] <= 2 + 10^-9,
      {p, 0, 1, 1/20}]];
  (* without spin-orbit the two spin channels are degenerate and never mix *)
  m0 = tmat[6, 0, 0, th6[1/4], 3/10, 1/2];
  checkTrue["no SOC: zero spin-flip transmission",
    Abs[m0[[1, 2]]] + Abs[m0[[2, 1]]] < 10^-12];
  checkNum["no SOC: the two spin channels are degenerate",
    m0[[1, 1]], m0[[2, 2]], 10^-12];
  m1 = tmat[6, 3/10, 0, th6[1/4], 3/10, 1/2];
  checkTrue["Rashba switches spin-flip transmission on",
    m1[[1, 2]] + m1[[2, 1]] > 10^-4];
];

Print[""];
Print["=============================================="];
Print[" 2. time reversal forbids spin polarisation"];
Print["=============================================="];
Print["  Two-terminal + time-reversal symmetry => the transmitted current"];
Print["  cannot be spin polarised, no matter how strong the SOC."];
Print[""];

Module[{p0, p1},
  p0 = pol[6, 3/10, 3/10, 0, 3/10, 1/2];
  Print["  P_z at zero flux, alpha = beta = 0.3 : ", p0];
  checkTrue["P_z = 0 at zero flux (TRS intact)", Abs[p0] < 10^-9];
  checkTrue["P_z = 0 at zero flux for several energies and SOC strengths",
    And @@ Flatten[Table[
       Abs[pol[6, av, bv, 0, 3/10, wv]] < 10^-9,
       {av, {0, 1/5, 1/2}}, {bv, {0, 3/10}}, {wv, {-1, -1/3, 1/2, 1}}]]];
  p1 = pol[6, 3/10, 0, th6[1/4], 3/10, 1/2];
  Print["  P_z at phi = 0.25, alpha = 0.3, beta = 0 : ", p1];
  checkTrue["flux + SOC together do polarise the current", Abs[p1] > 10^-4];
];

Print[""];
Print["=============================================="];
Print[" 2b. Rashba <-> Dresselhaus antisymmetry of P_z"];
Print["=============================================="];
Print["  It is observed (for persistent currents) that Rashba and"];
Print["  Dresselhaus spin currents are opposite, so alpha = beta should give"];
Print["  'a density very close to zero'. In transport the statement sharpens:"];
Print["  the map alpha <-> beta is a spin rotation that flips sigma_z, so"];
Print["            P_z(alpha, beta) = -P_z(beta, alpha)"];
Print["  identically -- and alpha = beta then forces P_z = 0 exactly, not"];
Print["  approximately."];
Print[""];

Module[{worst, worstEq, flip},
  worst = Max[Flatten[Table[
     Abs[pol[6, av, bv, th6[p], 3/10, wv] + pol[6, bv, av, th6[p], 3/10, wv]],
     {av, {1/10, 3/10, 7/10}}, {bv, {0, 1/5, 1/2}},
     {p, {1/8, 1/4, 2/5}}, {wv, {-11/10, -3/10, 7/10, 13/10}}]]];
  Print["  max |P_z(a,b) + P_z(b,a)| over 108 (a,b,phi,w) points : ", worst];
  checkTrue["P_z(alpha,beta) = -P_z(beta,alpha) identically", worst < 10^-12];

  worstEq = Max[Flatten[Table[
     Abs[pol[ns2, av, av, 2 Pi p/ns2, 3/10, wv]],
     {ns2, {4, 6, 8}}, {av, {1/10, 3/10, 7/10, 12/10}},
     {p, {1/8, 1/4, 2/5}}, {wv, {-11/10, -3/10, 7/10}}]]];
  Print["  max |P_z| at alpha = beta over 108 points, 3 ring sizes : ", worstEq];
  checkTrue["alpha = beta kills P_z exactly, for every size/flux/energy",
    worstEq < 10^-12];

  (* but the spin flip itself does NOT vanish: it is only symmetrised *)
  flip = tflip[6, 3/10, 3/10, th6[1/4], 3/10, 1/2];
  Print["  spin-flip transmission at alpha = beta : ", flip];
  checkTrue["spin flip survives at alpha = beta (only the NET polarisation dies)",
    flip > 10^-3];
];

Print[""];
Print["=============================================="];
Print[" 3. RESULTS: transport through the ring"];
Print["=============================================="];

Module[{ns = 20, ws, phis, als, data1, data2, data3, data4},
  ns = 20;                                     (* same size as her figures *)
  th20[phi_] := 2 Pi phi/ns;

  Print[""];
  Print[" 3a. charge transmission vs energy (phi = 0.25), her four cases"];
  ws = Range[-24, 24]/10;
  data1 = Table[
     {N[w],
      tcharge[ns, 0, 0, th20[1/4], 3/10, w],
      tcharge[ns, 3/10, 0, th20[1/4], 3/10, w],
      tcharge[ns, 0, 3/10, th20[1/4], 3/10, w],
      tcharge[ns, 3/10, 3/10, th20[1/4], 3/10, w]}, {w, ws}];
  Print["  peak T_c (no SOC)      = ", Max[data1[[All, 2]]]];
  Print["  peak T_c (Rashba)      = ", Max[data1[[All, 3]]]];
  Print["  peak T_c (Dresselhaus) = ", Max[data1[[All, 4]]]];
  Print["  peak T_c (both)        = ", Max[data1[[All, 5]]]];

  Print[""];
  Print[" 3b. Aharonov-Bohm oscillation of the conductance (w = 0.5)"];
  phis = Range[-100, 100]/100;
  data2 = Table[
     {N[p],
      tcharge[ns, 0, 0, th20[p], 3/10, 1/2],
      tcharge[ns, 3/10, 0, th20[p], 3/10, 1/2],
      tcharge[ns, 3/10, 3/10, th20[p], 3/10, 1/2]}, {p, phis}];
  Print["  T_c range, no SOC : ", {Min[data2[[All, 2]]], Max[data2[[All, 2]]]}];
  Print["  T_c range, Rashba : ", {Min[data2[[All, 3]]], Max[data2[[All, 3]]]}];

  Print[""];
  Print[" 3c. Aharonov-Casher: transmission vs SOC strength at zero flux"];
  als = Range[0, 150]/100;
  data3 = Table[
     {N[a],
      tcharge[ns, a, 0, 0, 3/10, 1/2],
      tflip[ns, a, 0, 0, 3/10, 1/2]}, {a, als}];
  Print["  T_c at alpha=0    : ", data3[[1, 2]]];
  Print["  T_c at alpha=0.75 : ", data3[[76, 2]]];
  Print["  T_c at alpha=1.5  : ", data3[[151, 2]]];
  Print["  max spin-flip T over the sweep: ", Max[data3[[All, 3]]]];

  Print[""];
  Print[" 3d. spin polarisation vs flux, and the alpha = beta question"];
  data4 = Table[
     {N[p],
      pol[ns, 3/10, 0, th20[p], 3/10, 1/2],
      pol[ns, 0, 3/10, th20[p], 3/10, 1/2],
      pol[ns, 3/10, 3/10, th20[p], 3/10, 1/2],
      tflip[ns, 3/10, 0, th20[p], 3/10, 1/2],
      tflip[ns, 3/10, 3/10, th20[p], 3/10, 1/2]}, {p, phis}];
  Print["  max |P_z| Rashba only      : ", Max[Abs[data4[[All, 2]]]]];
  Print["  max |P_z| Dresselhaus only : ", Max[Abs[data4[[All, 3]]]]];
  Print["  max |P_z| alpha = beta     : ", Max[Abs[data4[[All, 4]]]]];
  Print["  max spin-flip T, Rashba only  : ", Max[data4[[All, 5]]]];
  Print["  max spin-flip T, alpha = beta : ", Max[data4[[All, 6]]]];
  Print[""];
  Print["  Rashba vs Dresselhaus polarisation, point by point (phi>0 half):"];
  Print["    max |P_R(phi) + P_D(phi)| = ",
    Max[Abs[data4[[All, 2]] + data4[[All, 3]]]]];
  Print["    max |P_R(phi)|            = ", Max[Abs[data4[[All, 2]]]]];

  Export["/home/astrum/ring_Tvsw.csv", Prepend[data1,
     {"w", "T_none", "T_rashba", "T_dressel", "T_both"}]];
  Export["/home/astrum/ring_AB.csv", Prepend[data2,
     {"phi", "T_none", "T_rashba", "T_both"}]];
  Export["/home/astrum/ring_AC.csv", Prepend[data3,
     {"alpha", "T_charge", "T_flip"}]];
  Export["/home/astrum/ring_pol.csv", Prepend[data4,
     {"phi", "P_rashba", "P_dressel", "P_both", "Tflip_rashba", "Tflip_both"}]];
  Print[""];
  Print["  exported 4 CSV files"];
];

Print[""];
Print["=============================================="];
Print["  PASSED: ", $pass, "     FAILED: ", $fail];
Print["=============================================="];
