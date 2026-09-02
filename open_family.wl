(* The three rings of Table II, opened with two reservoirs.

   Everything reported for graphene and silicene so far is for the CLOSED ring
   (persistent currents), while the Rashba-Dresselhaus identity was checked both
   closed and in transport. This closes that asymmetry and asks the question a
   measurement would ask: with the leads attached, does the omitted Peierls
   phase show up in the transmitted current?

   One builder for the whole family, so the three models are parameter choices
   of Eq. (13) exactly as the paper claims. Conventions follow the manuscript:
   theta_{n,n+1} = (phi_n + phi_{n+1})/2 is the bond angle, and the flux phase
   is m*ph on a hop of range m -- with pei = 0 restoring the published form,
   which carries the phase only on the range-one hops.

   Leads: wide band, on diametrically opposite sites, coupling Gamma.        *)

Get["GreenEOM.wl"]; Get["GreenNEQ.wl"];

$pass = 0; $fail = 0;
sf[x_] := ToString[CForm[N[x, 4]]];
nf[x_] := ToString[CForm[N[x, 4]]];
checkTrue[n_String, v_] := If[TrueQ[v], $pass++; Print["  PASS  ", n],
  $fail++; Print["  FAIL  ", n, " -> ",
                 StringTake[ToString[InputForm[v]], UpTo[110]]]];

(* ------------------------------------------------------------------ *)
(* the family Hamiltonian as a single-particle matrix                   *)
(* ------------------------------------------------------------------ *)
ix[n_, s_] := 2 (n - 1) + If[s === 1, 1, 2];

famH[ns_, t_, al_, be_, lEO_, lR1_, lR2_, ph_, pei_] :=
 Module[{m = ConstantArray[0, {2 ns, 2 ns}], add, nx, n2, site, bond, p2},
  site[n_] := 2 Pi (n - 1)/ns;
  nx[n_] := Mod[n, ns] + 1;
  n2[n_] := Mod[n + 1, ns] + 1;
  (* midpoint of the arc a range-k hop spans. Writing it as an offset
     from phi_n rather than as an average of the two site angles is what
     keeps the closing bond, where phi_1 would be read as 0 instead of
     2 Pi, from coming out Pi away from every other bond. *)
  bond[n_, k_] := site[n] + k Pi/ns;
  add[i_, j_, v_] := (m[[i, j]] += v; m[[j, i]] += Conjugate[v]);
  p2 := If[pei === 1, Exp[2 I ph], 1];                  (* range-two phase *)
  Do[
   (* range one: hopping, Rashba, Dresselhaus *)
   add[ix[n, s], ix[nx[n], s], t Exp[I ph]];
   add[ix[n, s], ix[nx[n], -s],
       (-I al s Exp[-I s bond[n, 1]] + be Exp[I s bond[n, 1]]) Exp[I ph]];
   add[ix[n, s], ix[nx[n], -s], -(lR1/2) s Exp[-I s bond[n, 1]] Exp[I ph]];
   (* range two: intrinsic spin-orbit, intrinsic Rashba *)
   add[ix[n, s], ix[n2[n], s], I lEO (s/2) p2];
   add[ix[n, s], ix[n2[n], -s], -(lR2/2) s Exp[-I s bond[n, 2]] p2],
   {n, ns}, {s, {1, -1}}];
  m];

(* ------------------------------------------------------------------ *)
(* two wide-band leads on diametrically opposite sites                  *)
(* ------------------------------------------------------------------ *)
$eta = 10^-8;
prj[ns_, n_, s_] := prj[ns, n, s] =
  Normal[SparseArray[{{ix[n, s], ix[n, s]} -> 1}, {2 ns, 2 ns}]];
gL[ns_, g_] := gL[ns, g] = g (prj[ns, 1, 1] + prj[ns, 1, -1]);
gR[ns_, g_] := gR[ns, g] = With[{nr = 1 + Quotient[ns, 2]},
   g (prj[ns, nr, 1] + prj[ns, nr, -1])];

tmat[ns_, pars__, ph_, pei_, g_, w_] :=
 Module[{h, sig, gr, ga, nr},
  h = N[famH[ns, pars, ph, pei]];
  nr = 1 + Quotient[ns, 2];
  sig = WideBandSelfEnergy[gL[ns, g] + gR[ns, g]];
  gr = RetardedG[h, sig, N[w] + I $eta];
  ga = Dagger[gr];
  Table[Re[Tr[(g prj[ns, nr, so]) . gr . (g prj[ns, 1, si]) . ga]],
   {so, {1, -1}}, {si, {1, -1}}]];

tchg[a__] := Total[tmat[a], 2];
pz[a__] := With[{m = tmat[a]},
  (m[[1, 1]] + m[[1, 2]] - m[[2, 1]] - m[[2, 2]])/Total[m, 2]];

(* parameter points: the three models of Table II *)
ns = 6; gam = 3/10; tt = 1;
rd[al_, be_, ph_] := {tt, al, be, 0, 0, 0, ph};          (* Rashba-Dresselhaus *)
gr[l1_, ph_] := {tt, 0, 0, 39/1000, l1, 0, ph};          (* graphene          *)
si[l1_, l2_, ph_] := {tt, 0, 0, 39/1000, l1, l2, ph};    (* silicene          *)

Print["=============================================="];
Print[" The three rings of Table II, opened"];
Print["=============================================="];
Print["  N = ", ns, ", leads on sites 1 and ", 1 + Quotient[ns, 2],
      ", wide band, Gamma = ", N[gam], ", t = 1"];
Print[""];

(* --- 1. the machinery still behaves ------------------------------- *)
Print[" 1. consistency"];
Module[{a, b},
  a = tchg[ns, Sequence @@ rd[0, 0, 2 Pi (1/4)/ns], 1, gam, 1/10];
  b = tchg[ns, Sequence @@ rd[0, 0, -2 Pi (1/4)/ns], 1, gam, 1/10];
  Print["    T(phi) - T(-phi), no SOC : ", sf[Abs[a - b]]];
  checkTrue["two-terminal reciprocity T(phi) = T(-phi)", Abs[a - b] < 10^-9];
  checkTrue["transmission is bounded by the number of channels",
            0 <= a <= 2 + 10^-9]];
Print[""];

(* --- 2. the alpha <-> beta mirror, in transport -------------------- *)
Print[" 2. Rashba <-> Dresselhaus, with leads"];
Module[{res},
  res = Table[
    Module[{p = 2 Pi q/ns, u, v},
      u = pz[ns, Sequence @@ rd[3/10, 1/10, p], 1, gam, w];
      v = pz[ns, Sequence @@ rd[1/10, 3/10, p], 1, gam, w];
      Abs[u + v]],
    {q, {1/8, 1/4, 2/5}}, {w, {-8/10, -2/10, 3/10, 9/10}}];
  Print["    max |P_z(al,be) + P_z(be,al)| : ",
        sf[Max[Flatten[res]]]];
  checkTrue["the mirror survives with reservoirs attached",
            Max[Flatten[res]] < 10^-9]];
Print[""];

(* --- 3. is there an analogue in the other two? -------------------- *)
Print[" 3. searching graphene/silicene for an analogue"];
Module[{res, mn},
  res = Table[
    Module[{p = 2 Pi q/ns, u, v},
      u = pz[ns, Sequence @@ si[l, m, p], 1, gam, w];
      v = pz[ns, Sequence @@ si[m, l, p], 1, gam, w];
      Abs[u + v]],
    {q, {1/8, 1/4, 2/5}}, {w, {-8/10, -2/10, 3/10, 9/10}},
    {l, {5/100, 15/100, 3/10}}, {m, {5/100, 15/100, 3/10}}];
  mn = Min[Flatten[res]];
  Print["    lambda_R1 <-> lambda_R2 exchange, ", Length[Flatten[res]],
        " points"];
  Print["    min |P_z(l1,l2) + P_z(l2,l1)| : ", sf[mn]];
  checkTrue["no lambda_R1 <-> lambda_R2 mirror exists", mn > 10^-3]];
Print[""];

(* --- 4. does the omitted phase show up in the transmitted current? - *)
Print[" 4. graphene vs silicene with leads, phase omitted and restored"];
Print["    (closed ring: J_s was blind to lambda_R2 at half filling, N=2 mod 4)"];
Module[{rows, sp},
  sp[v_] := Max[Abs[v - First[v]]];
  rows = Table[
    Module[{p = 2 Pi q/ns},
      {q, w,
       sp[Table[tchg[ns, Sequence @@ si[12/100, x, p], 0, gam, w],
                {x, 0, 3/10, 1/20}]],
       sp[Table[tchg[ns, Sequence @@ si[12/100, x, p], 1, gam, w],
                {x, 0, 3/10, 1/20}]],
       sp[Table[pz[ns, Sequence @@ si[12/100, x, p], 0, gam, w],
                {x, 0, 3/10, 1/20}]],
       sp[Table[pz[ns, Sequence @@ si[12/100, x, p], 1, gam, w],
                {x, 0, 3/10, 1/20}]]}],
    {q, {1/8, 1/4, 2/5}}, {w, {-8/10, -2/10, 3/10, 9/10}}];
  rows = Flatten[rows, 1];
  Print["    phi      w      dT(no ph)   dT(ph)    dPz(no ph)  dPz(ph)"];
  Do[Print["    ", nf[r[[1]]], "  ",
           nf[r[[2]]], "  ",
           sf[r[[3]]], "  ", sf[r[[4]]], "  ",
           sf[r[[5]]], "  ", sf[r[[6]]]],
    {r, rows}];
  Print[""];
  Print["    worst-case blindness to lambda_R2, phase omitted:"];
  Print["      min over points of dT  : ",
        sf[Min[rows[[All, 3]]]]];
  Print["      min over points of dPz : ",
        sf[Min[rows[[All, 5]]]]];
  checkTrue["with leads, lambda_R2 acts on the charge transmission even "
            <> "without its phase", Min[rows[[All, 3]]] > 10^-6]];

Print[""];
Print["  ", $pass, " passed, ", $fail, " failed"];
