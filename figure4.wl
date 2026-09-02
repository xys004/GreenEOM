(* Data for the figure that separates the graphene ring from the silicene one.

   lambda_R2 is the coupling that distinguishes the two: at lambda_R2 = 0 the
   model is Emma's graphene ring, and lambda_R2 > 0 is Dayanna's silicene ring.
   The point of the figure is that without the Peierls phase on the range-two
   hops the spin current cannot tell them apart at half filling of a ring with
   N = 2 (mod 4), which is the case Sec. III B says is "a natural one to plot".

   Same hm/en/js as scans.wl, so this shares its source of truth with the
   numbers already quoted in Sec. III B. Emits CSV on stdout. *)

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

(* the parameter point of Sec. III B *)
lEO = 39/1000; lR1 = 12/100; p = 1/4;
xs = Table[x, {x, 0, 3/10, 1/200}];

row[ns_] := Module[{ph = 2 Pi p/ns, nf = ns},
  Table[{N[x],
         js[ns, lEO, lR1, x, ph, 0, nf],
         js[ns, lEO, lR1, x, ph, 1, nf]}, {x, xs}]];

d10 = row[10];
d8  = row[8];

Print["lambda_R2,js10_nophase,js10_phase,js8_nophase,js8_phase"];
Do[Print[StringRiffle[
    ToString[CForm[N[#, 10]]] & /@ {d10[[i, 1]], d10[[i, 2]], d10[[i, 3]],
                                    d8[[i, 2]], d8[[i, 3]]}, ","]],
  {i, Length[xs]}];

sp[v_] := Max[Abs[v - First[v]]];
Print["# N=10 half filling: span without phase = ",
      ScientificForm[sp[d10[[All, 2]]], 3],
      " , with phase = ", ScientificForm[sp[d10[[All, 3]]], 3]];
Print["# N=8  half filling: span without phase = ",
      ScientificForm[sp[d8[[All, 2]]], 3],
      " , with phase = ", ScientificForm[sp[d8[[All, 3]]], 3]];
