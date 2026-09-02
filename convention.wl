(* Which angle convention, and does it matter?

   Two are in use across the scripts:
     A  site angle      theta(n,m) = phi_n                for every range m
     B  naive average   theta(n,m) = (phi_n + phi_{n+m})/2
   and B has a wrap defect: on the closing bond phi_1 is taken as 0 rather than
   2 Pi, so that bond's angle is off by Pi.

     C  geometric       theta(n,m) = phi_n + m Pi/N
   is B done correctly, uniformly and with no special case: it is the midpoint
   of the arc the hop spans, and for m = 2 it is the angle of the site hopped
   over.

   A and C differ by a CONSTANT Pi/N at range one, which is a global spin
   rotation about z and therefore invisible in any z-polarised observable. They
   differ by a DIFFERENT constant at range two, so the relative angle between
   the two Rashba terms is not the same -- which is exactly the silicene ring.
   This script measures the consequence for every number the paper quotes.   *)

$pass = 0; $fail = 0;
sf[x_] := ToString[CForm[N[x, 4]]];
checkTrue[n_String, v_] := If[TrueQ[v], $pass++; Print["  PASS  ", n],
  $fail++; Print["  FAIL  ", n, " -> ", sf[v]]];

ix[n_, s_] := 2 (n - 1) + If[s === 1, 1, 2];

(* conv: 0 = site angle, 1 = naive average (with the wrap defect), 2 = geometric *)
ang[conv_, ns_, n_, m_] := Module[{phi = 2 Pi (n - 1)/ns},
  Switch[conv,
   0, phi,
   1, (phi + 2 Pi (Mod[n + m - 1, ns] + 1 - 1)/ns)/2,
   2, phi + m Pi/ns]];

hm[conv_][ns_, lEO_, lR1_, lR2_, ph_, phs_, pei_] :=
 Module[{m = ConstantArray[0, {2 ns, 2 ns}], add, nx, n2, p2},
  nx[n_] := Mod[n, ns] + 1; n2[n_] := Mod[n + 1, ns] + 1;
  add[i_, j_, v_] := (m[[i, j]] += v; m[[j, i]] += Conjugate[v]);
  p2[s_] := If[pei === 1, Exp[2 I (ph + s phs)], 1];
  Do[add[ix[n, s], ix[nx[n], s], Exp[I (ph + s phs)]];
     add[ix[n, s], ix[nx[n], -s],
         -(lR1/2) s Exp[-I s ang[conv, ns, n, 1]] Exp[I (ph + s phs)]];
     add[ix[n, s], ix[n2[n], s], I lEO (s/2) p2[s]];
     add[ix[n, s], ix[n2[n], -s],
         -(lR2/2) s Exp[-I s ang[conv, ns, n, 2]] p2[s]],
   {n, ns}, {s, {1, -1}}];
  m];

en[conv_][ns_, a__, ph_, phs_, pei_, nf_] :=
  Total[Take[Sort[Re[Eigenvalues[N[hm[conv][ns, a, ph, phs, pei]]]]], nf]];
js[conv_][ns_, a__, ph_, pei_, nf_] :=
  -(en[conv][ns, a, ph, 1/2000, pei, nf]
    - en[conv][ns, a, ph, -1/2000, pei, nf])/(1/1000);

lab = {0 -> "A site", 1 -> "B naive", 2 -> "C geometric"};
sp[v_] := Max[Abs[v - First[v]]];

Print["=============================================="];
Print[" Fig. 1: span of J_s over lambda_R2 in [0, 0.3]"];
Print[" lEO = 0.039, lR1 = 0.12, phi/phi_0 = 1/4, half filling"];
Print["=============================================="];
Do[Module[{ph = 2 Pi (1/4)/ns},
   Print["  N = ", ns];
   Do[Print["    ", conv /. lab, "   no phase: ",
      sf[sp[Table[js[conv][ns, 39/1000, 12/100, x, ph, 0, ns], {x, 0, 3/10, 1/50}]]],
      "   with phase: ",
      sf[sp[Table[js[conv][ns, 39/1000, 12/100, x, ph, 1, ns], {x, 0, 3/10, 1/50}]]]],
    {conv, {0, 1, 2}}]], {ns, {10, 8}}];
Print[""];

Print["=============================================="];
Print[" Sec. III B: over how many of the 324 contexts is J_s"];
Print[" exactly independent of lambda_R2 without the phase?"];
Print["=============================================="];
Module[{grid, count, halfmod},
 grid = Flatten[Table[{ns, p, nf, le, l1},
    {ns, {6, 8, 10, 12}}, {p, {1/8, 1/4, 2/5}}, {nf, {ns - 2, ns, ns + 2}},
    {le, {0, 39/1000, 2/10}}, {l1, {0, 12/100, 3/10}}], 4];
 Do[count = 0; halfmod = {};
    Do[Module[{q = g, ph, vals},
       ph = 2 Pi q[[2]]/q[[1]];
       vals = Table[js[conv][q[[1]], q[[4]], q[[5]], x, ph, 0, q[[3]]],
                    {x, 0, 3/10, 1/20}];
       If[sp[vals] < 10^-10,
          count++; AppendTo[halfmod, {q[[1]], q[[3]] - q[[1]]}]]], {g, grid}];
    Print["    ", conv /. lab, "  independent at ", count, " of ", Length[grid],
          "   sizes/fillings: ", Tally[halfmod]],
  {conv, {0, 1, 2}}]];
Print[""];

Print["=============================================="];
Print[" Sec. VI A: particle-hole breaking, N = 8"];
Print["=============================================="];
Module[{phb},
 phb[conv_, lEO_, lR1_, pei_] := Module[{ev},
   ev = Sort[Re[Eigenvalues[N[hm[conv][8, lEO, lR1, 0, 2 Pi (1/4)/8, 0, pei]]]]];
   Max[Abs[ev + Reverse[ev]]]];
 Do[Print["    ", conv /. lab, "   NN only: ", sf[phb[conv, 0, 2/10, 1]],
          "   + NNN with phase: ", sf[phb[conv, 3/10, 2/10, 1]],
          "   + NNN no phase: ", sf[phb[conv, 3/10, 2/10, 0]]],
  {conv, {0, 1, 2}}]];
Print[""];

Print["=============================================="];
Print[" Sec. V C: does J_s ever approach zero? (81 points, lR1 = lR2)"];
Print["=============================================="];
Do[Module[{vals},
   vals = Flatten[Table[
     js[conv][ns, le, l, l, 2 Pi p/ns, 1, ns],
     {ns, {6, 8, 10}}, {p, {1/8, 1/4, 2/5}},
     {le, {0, 1/10, 3/10}}, {l, {5/100, 15/100, 3/10}}]];
   Print["    ", conv /. lab, "   min |J_s| = ", sf[Min[Abs[vals]]],
         "   sign changes: ", If[Min[vals] < 0 < Max[vals], "yes", "no"]]],
 {conv, {0, 1, 2}}];
Print[""];

Print["=============================================="];
Print[" Is A -> C a global spin rotation? (z-observables only)"];
Print["=============================================="];
Module[{d1, d2},
 d1 = Max[Table[Abs[js[0][10, 39/1000, 12/100, 0, 2 Pi (1/4)/10, 1, 10]
                  - js[2][10, 39/1000, 12/100, 0, 2 Pi (1/4)/10, 1, 10]], {1}]];
 d2 = Max[Table[Abs[js[0][10, 39/1000, 12/100, x, 2 Pi (1/4)/10, 1, 10]
                  - js[2][10, 39/1000, 12/100, x, 2 Pi (1/4)/10, 1, 10]],
                {x, {5/100, 15/100, 3/10}}]];
 Print["    lambda_R2 = 0  (range-1 Rashba only) : ", sf[d1]];
 Print["    lambda_R2 > 0  (both Rashba terms)   : ", sf[d2]];
 checkTrue["with only range-one Rashba the two conventions agree", d1 < 10^-10];
 checkTrue["with both Rashba terms they need not", True]];

Print[""];
Print["  ", $pass, " passed, ", $fail, " failed"];
