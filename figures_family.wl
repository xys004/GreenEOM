(* Data for the two comparative figures.

   (1) the three models side by side, persistent charge and spin current
       against flux, with the range-two phase omitted and restored wherever
       there is a range-two term to omit it from;
   (2) the same family with two reservoirs attached: transmission against
       injection energy, and against lambda_R2 at fixed energy.

   One builder for the whole of Table II, in the canonical bond angle
   theta^(m)_n = phi_n + m Pi/N of Eq. (14). Writes two CSV files.        *)

Get["GreenEOM.wl"]; Get["GreenNEQ.wl"];

ix[n_, s_] := 2 (n - 1) + If[s === 1, 1, 2];
ang[ns_, n_, m_] := 2 Pi (n - 1)/ns + m Pi/ns;

famH[ns_, al_, be_, lEO_, lR1_, lR2_, ph_, phs_, pei_] :=
 Module[{m = ConstantArray[0, {2 ns, 2 ns}], add, nx, n2, p1, p2},
  nx[n_] := Mod[n, ns] + 1; n2[n_] := Mod[n + 1, ns] + 1;
  add[i_, j_, v_] := (m[[i, j]] += v; m[[j, i]] += Conjugate[v]);
  p1[s_] := Exp[I (ph + s phs)];
  p2[s_] := If[pei === 1, Exp[2 I (ph + s phs)], 1];
  Do[add[ix[n, s], ix[nx[n], s], p1[s]];
     add[ix[n, s], ix[nx[n], -s],
         (-I al s Exp[-I s ang[ns, n, 1]] + be Exp[I s ang[ns, n, 1]]) p1[s]];
     add[ix[n, s], ix[nx[n], -s],
         -(lR1/2) s Exp[-I s ang[ns, n, 1]] p1[s]];
     add[ix[n, s], ix[n2[n], s], I lEO (s/2) p2[s]];
     add[ix[n, s], ix[n2[n], -s], -(lR2/2) s Exp[-I s ang[ns, n, 2]] p2[s]],
   {n, ns}, {s, {1, -1}}];
  m];

en[ns_, p__, ph_, phs_, pei_, nf_] :=
  Total[Take[Sort[Re[Eigenvalues[N[famH[ns, p, ph, phs, pei]]]]], nf]];
d = 1/2000;
jc[ns_, p__, ph_, pei_, nf_] :=
  -(en[ns, p, ph + d, 0, pei, nf] - en[ns, p, ph - d, 0, pei, nf])/(2 d);
jsz[ns_, p__, ph_, pei_, nf_] :=
  -(en[ns, p, ph, d, pei, nf] - en[ns, p, ph, -d, pei, nf])/(2 d);

(* the three models of Table II, at couplings large enough to see *)
rd  = {3/10, 1/10, 0, 0, 0};              (* alpha, beta                *)
gra = {0, 0, 2/10, 12/100, 0};            (* lEO, lR1                   *)
sil = {0, 0, 2/10, 12/100, 15/100};       (* lEO, lR1, lR2              *)

(* ---------------- 1. closed ring: currents against flux ------------- *)
Module[{ns = 10, nf, rows, f},
 nf = ns;
 rows = Table[
   Module[{ph = 2 Pi q/ns},
    {N[q],
     jc[ns, Sequence @@ rd, ph, 1, nf],  jsz[ns, Sequence @@ rd, ph, 1, nf],
     jc[ns, Sequence @@ gra, ph, 1, nf], jsz[ns, Sequence @@ gra, ph, 1, nf],
     jc[ns, Sequence @@ gra, ph, 0, nf], jsz[ns, Sequence @@ gra, ph, 0, nf],
     jc[ns, Sequence @@ sil, ph, 1, nf], jsz[ns, Sequence @@ sil, ph, 1, nf],
     jc[ns, Sequence @@ sil, ph, 0, nf], jsz[ns, Sequence @@ sil, ph, 0, nf]}],
   {q, -1, 1, 1/150}];
 f = OpenWrite["data/family_flux.csv"];
 WriteString[f, "phi,rd_jc,rd_js,gra_jc,gra_js,graNP_jc,graNP_js,"
                <> "sil_jc,sil_js,silNP_jc,silNP_js\n"];
 Do[WriteString[f, StringRiffle[ToString[CForm[N[#, 10]]] & /@ r, ","], "\n"],
  {r, rows}];
 Close[f];
 Print["data/family_flux.csv  ", Length[rows], " rows"]];

(* ---------------- 2. open ring: two wide-band reservoirs ------------ *)
$eta = 10^-8;
prj[ns_, n_, s_] := prj[ns, n, s] =
  Normal[SparseArray[{{ix[n, s], ix[n, s]} -> 1}, {2 ns, 2 ns}]];
gam[ns_, g_] := gam[ns, g] = With[{nr = 1 + Quotient[ns, 2]},
  g (prj[ns, 1, 1] + prj[ns, 1, -1] + prj[ns, nr, 1] + prj[ns, nr, -1])];

tm[ns_, p__, ph_, pei_, g_, w_] :=
 Module[{h, sig, gr, ga, nr},
  h = N[famH[ns, p, ph, 0, pei]];
  nr = 1 + Quotient[ns, 2];
  sig = WideBandSelfEnergy[gam[ns, g]];
  gr = RetardedG[h, sig, N[w] + I $eta];
  ga = Dagger[gr];
  Table[Re[Tr[(g prj[ns, nr, so]) . gr . (g prj[ns, 1, si]) . ga]],
   {so, {1, -1}}, {si, {1, -1}}]];
tt[a__] := Total[tm[a], 2];
pz[a__] := With[{m = tm[a]},
  (m[[1, 1]] + m[[1, 2]] - m[[2, 1]] - m[[2, 2]])/Total[m, 2]];

Module[{ns = 6, g = 3/10, ph, rows, f},
 ph = 2 Pi (1/4)/ns;
 (* 2a. transmission against injection energy, the three models *)
 rows = Table[{N[w],
    tt[ns, Sequence @@ rd, ph, 1, g, w],
    tt[ns, Sequence @@ gra, ph, 1, g, w],
    tt[ns, Sequence @@ sil, ph, 1, g, w],
    tt[ns, Sequence @@ sil, ph, 0, g, w]},
   {w, -5/2, 5/2, 1/200}];
 f = OpenWrite["data/family_transmission.csv"];
 WriteString[f, "w,rd,gra,sil,silNP\n"];
 Do[WriteString[f, StringRiffle[ToString[CForm[N[#, 10]]] & /@ r, ","], "\n"],
  {r, rows}];
 Close[f];
 Print["data/family_transmission.csv  ", Length[rows], " rows"];

 (* 2b. the open counterpart of Fig. 1: response to lambda_R2 *)
 rows = Table[{N[x],
    tt[ns, 0, 0, 2/10, 12/100, x, ph, 0, g, 3/10],
    tt[ns, 0, 0, 2/10, 12/100, x, ph, 1, g, 3/10],
    pz[ns, 0, 0, 2/10, 12/100, x, ph, 0, g, 3/10],
    pz[ns, 0, 0, 2/10, 12/100, x, ph, 1, g, 3/10]},
   {x, 0, 3/10, 1/400}];
 f = OpenWrite["data/open_lR2.csv"];
 WriteString[f, "lR2,T_nophase,T_phase,Pz_nophase,Pz_phase\n"];
 Do[WriteString[f, StringRiffle[ToString[CForm[N[#, 10]]] & /@ r, ","], "\n"],
  {r, rows}];
 Close[f];
 Print["data/open_lR2.csv  ", Length[rows], " rows"]];
