(* Puente con el PRB de Hubbard-Rashba: ellos senalan que la simetria
   particula-hueco se estropea con saltos a segundos vecinos. Ese es
   exactamente el termino que anaden los anillos de grafeno y siliceno. *)
$pass=0; $fail=0;
checkTrue[n_String,v_]:=If[TrueQ[v],$pass++;Print["  PASS  ",n],$fail++;
  Print["  FAIL  ",n," -> ",StringTake[ToString[InputForm[v]],UpTo[120]]]];
ns=8; nxt[n_]:=Mod[n,ns]+1; nx2[n_]:=Mod[n+1,ns]+1; idx[n_,s_]:=2(n-1)+If[s===1,1,2];
(* canonical bond angle: midpoint of the arc a range-m hop spans *)
thm[n_,m_]:=2 Pi (n-1)/ns + m Pi/ns;

hm[lEO_,lR_,ph_,pei_]:=Module[{m=ConstantArray[0,{2 ns,2 ns}],add,p2},
  add[i_,j_,v_]:=(m[[i,j]]+=v; m[[j,i]]+=Conjugate[v]);
  p2=If[pei===1,Exp[2 I ph],1];
  Do[add[idx[n,s],idx[nxt[n],s],Exp[I ph]];
     add[idx[n,s],idx[nxt[n],-s],-(lR/2) s Exp[-I s thm[n,1]] Exp[I ph]];
     add[idx[n,s],idx[nx2[n],s],I lEO (s/2) p2],{n,ns},{s,{1,-1}}];
  m];

(* simetria particula-hueco del espectro: el conjunto {E} debe ser simetrico
   bajo E -> -E *)
phBreak[lEO_,lR_,ph_,pei_]:=Module[{ev},
  ev=Sort[Re[Eigenvalues[N[hm[lEO,lR,ph,pei]]]]];
  Max[Abs[ev+Reverse[ev]]]];

Print["  Ruptura de la simetria E <-> -E  (0 = simetrico):"];
(* Flujo por enlace, 2 Pi (phi/phi_0)/N, el mismo que usa jc mas abajo: el
   script mezclaba antes este convenio con ph=1/4 crudo. Los numeros de la
   Sec. VI A son los de phi/phi_0 = 1/4. *)
ph0 = 2 Pi (1/4)/ns;
Print["    solo NN (su PRB), lR=0.2            : ",phBreak[0,2/10,ph0,1]];
Print["    + NNN con fase de Peierls, lEO=0.3  : ",phBreak[3/10,2/10,ph0,1]];
Print["    + NNN SIN fase, lEO=0.3 : ",phBreak[3/10,2/10,ph0,0]];
checkTrue["NN puro conserva la simetria particula-hueco",phBreak[0,2/10,ph0,1]<10^-10];
checkTrue["el termino a segundos vecinos la rompe (con fase)",phBreak[3/10,2/10,ph0,1]>10^-3];
checkTrue["y la rompe igual SIN fase",phBreak[3/10,2/10,ph0,0]>10^-3];

(* pero la corriente solo ve el termino si lleva la fase *)
en[lEO_,lR_,ph_,pei_]:=Total[Take[Sort[Re[Eigenvalues[N[hm[lEO,lR,ph,pei]]]]],ns]];
jc[lEO_,lR_,ph_,pei_]:=-(en[lEO,lR,ph+1/2000,pei]-en[lEO,lR,ph-1/2000,pei])/(1/1000);
Print[""];
Print["  Contribucion del termino de segundos vecinos a la corriente:"];
Module[{a,b,c},
  a=jc[0,2/10,2 Pi (1/4)/ns,1]; b=jc[3/10,2/10,2 Pi (1/4)/ns,1];
  c=jc[3/10,2/10,2 Pi (1/4)/ns,0];
  Print["    J_c sin NNN                 : ",a];
  Print["    J_c con NNN + fase Peierls  : ",b,"   (cambia en ",Abs[b-a],")"];
  Print["    J_c con NNN sin fase        : ",c,"   (cambia en ",Abs[c-a],")"];
  checkTrue["con fase, el NNN cambia la corriente",Abs[b-a]>10^-4];
  checkTrue["sin fase, tambien la cambia (via el espectro ocupado)",Abs[c-a]>10^-6];
  Print["    -> ojo: sin fase el termino SI altera el espectro y por tanto la"];
  Print["       corriente, pero no aporta ningun canal dependiente de flujo"];
  Print["       propio. Los dos efectos hay que distinguirlos."]];
Print[""];
Print["  PASSED: ",$pass,"   FAILED: ",$fail];
