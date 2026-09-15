(* ¿Que caracteriza a los contextos donde J_s no depende de lambda_R2? *)
hm[ns_,lEO_,lR1_,lR2_,ph_,phs_,pei_]:=Module[{m=ConstantArray[0,{2 ns,2 ns}],add,nx,n2,thm,p2,ix},
  ix[n_,s_]:=2(n-1)+If[s===1,1,2];
  nx[n_]:=Mod[n,ns]+1; n2[n_]:=Mod[n+1,ns]+1;
  (* canonical bond angle: midpoint of the arc a range-m hop spans,
     written so the closing bond needs no special case *)
  thm[n_,m_]:=2 Pi (n-1)/ns + m Pi/ns;
  add[i_,j_,v_]:=(m[[i,j]]+=v; m[[j,i]]+=Conjugate[v]);
  p2[s_]:=If[pei===1,Exp[2 I (ph+s phs)],1];
  Do[add[ix[n,s],ix[nx[n],s],Exp[I(ph+s phs)]];
     add[ix[n,s],ix[nx[n],-s],-(lR1/2) s Exp[-I s thm[n,1]] Exp[I(ph+s phs)]];
     add[ix[n,s],ix[n2[n],s], I lEO (s/2) p2[s]];
     add[ix[n,s],ix[n2[n],-s],-(lR2/2) s Exp[-I s thm[n,2]] p2[s]],
   {n,ns},{s,{1,-1}}];
  m];
en[ns_,a__,ph_,phs_,pei_,nf_]:=Total[Take[Sort[Re[Eigenvalues[N[hm[ns,a,ph,phs,pei]]]]],nf]];
js[ns_,a__,ph_,pei_,nf_]:=-(en[ns,a,ph,1/2000,pei,nf]-en[ns,a,ph,-1/2000,pei,nf])/(1/1000);

grid=Flatten[Table[{ns,p,nf,le,l1},
   {ns,{6,8,10,12}},{p,{1/8,1/4,2/5}},{nf,{ns-2,ns,ns+2}},
   {le,{0,39/1000,2/10}},{l1,{0,12/100,3/10}}],4];
flat=Select[grid, Module[{q=#,base,vals},
   base=js[q[[1]],q[[4]],q[[5]],0,2 Pi q[[2]]/q[[1]],0,q[[3]]];
   vals=Table[js[q[[1]],q[[4]],q[[5]],x,2 Pi q[[2]]/q[[1]],0,q[[3]]],{x,0,3/10,1/20}];
   Max[Abs[vals-base]]<10^-10]&];
Print["contextos con independencia exacta: ",Length[flat]," de ",Length[grid]];
Print[];
Print["  valores de lR1 presentes : ",Tally[flat[[All,5]]]];
Print["  valores de lEO presentes : ",Tally[flat[[All,4]]]];
Print["  llenados presentes       : ",Tally[MapThread[#1-#2&,{flat[[All,3]],flat[[All,1]]}]]];
Print["  tamanos presentes        : ",Tally[flat[[All,1]]]];
Print["  flujos presentes         : ",Tally[flat[[All,2]]]];
Print[];
Print["  (llenado se muestra como nf - N: 0 = semi-llenado)"];
