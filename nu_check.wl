(* Does the nu convention matter? On a 1D ring there is no honeycomb chirality,
   so nu alternating (inherited from graphene) and nu constant (all second-
   neighbour hops go the same way around the ring) are both defensible. *)
$pass=0; $fail=0;
checkTrue[n_String,v_]:=If[TrueQ[v],$pass++;Print["  PASS  ",n],$fail++;Print["  FAIL  ",n]];
DeclareSpecies[c,"Fermion"];
ns=12; nxt[n_]:=Mod[n,ns]+1; nx2[n_]:=Mod[n+1,ns]+1; th[n_]:=2 Pi (n-1)/ns;
hcOf[e_]:=Expand[e]/.cc_. Mono[{Cre[sp_,i_],Ann[sp2_,j_]}]:>
  (cc/.Complex[re_,im_]:>Complex[re,-im]) Mono[{Cre[sp2,j],Ann[sp,i]}];
herm[e_]:=Expand[NCTimes[e]+hcOf[NCTimes[e]]];
hB[t_,ph_]:=t Sum[Exp[I ph] NCTimes[Cre[c,{n,s}],Ann[c,{nxt[n],s}]]+
   Exp[-I ph] NCTimes[Cre[c,{nxt[n],s}],Ann[c,{n,s}]],{n,ns},{s,{1,-1}}];
hCg[lam_,nuf_]:=herm[I lam Sum[nuf[n] (s/2) NCTimes[Cre[c,{n,s}],Ann[c,{nx2[n],s}]],
   {n,ns},{s,{1,-1}}]];
hD[lr_]:=herm[-(lr/2) Sum[s Exp[-I s th[n]] NCTimes[Cre[c,{n,s}],Ann[c,{nxt[n],-s}]],
   {n,ns},{s,{1,-1}}]];
ops=Flatten[Table[Ann[c,{n,s}],{n,ns},{s,{1,-1}}]];
nuAlt[n_]:=(-1)^(n-1); nuConst[n_]:=1;
subs={tt->1,ph->1/4,lEO->3/10,lR1->1/5};
mA=N[HamiltonianMatrix[hB[tt,ph]+hCg[lEO,nuAlt]+hD[lR1],ops]/.subs];
mC=N[HamiltonianMatrix[hB[tt,ph]+hCg[lEO,nuConst]+hD[lR1],ops]/.subs];
evA=Sort[Re[Eigenvalues[mA]]]; evC=Sort[Re[Eigenvalues[mC]]];
Print["  both Hermitian?  alt: ",Max[Abs[mA-ConjugateTranspose[mA]]]<10^-12,
      "   const: ",Max[Abs[mC-ConjugateTranspose[mC]]]<10^-12];
checkTrue["both nu conventions give a Hermitian H",
  Max[Abs[mA-ConjugateTranspose[mA]]]<10^-12 && Max[Abs[mC-ConjugateTranspose[mC]]]<10^-12];
Print["  spectrum, nu alternating : ",Take[evA,4]," ..."];
Print["  spectrum, nu constant    : ",Take[evC,4]," ..."];
Print["  max |difference|         : ",Max[Abs[evA-evC]]];
checkTrue["the two conventions are NOT equivalent", Max[Abs[evA-evC]]>10^-3];
(* and the missing-term conclusion is independent of the convention *)
Do[Module[{cm,idx},
  cm=Comm[Ann[c,{5,1}],hB[tt,ph]+hCg[lEO,nuf]+hD[lR1]];
  idx=Sort[DeleteDuplicates[Cases[cm,Ann[c,i_]:>i,{0,Infinity}]]];
  Print["  nu ",If[nuf===nuAlt,"alternating","constant  "],
        " -> ",Length[idx]," Green functions, includes site n-1 = 3: ",
        MemberQ[idx,{3,1}]];
  checkTrue["backward hop present regardless of nu convention",MemberQ[idx,{3,1}]]],
 {nuf,{nuAlt,nuConst}}];
Print["  PASSED: ",$pass,"   FAILED: ",$fail];
