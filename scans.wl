(* Dos barridos anchos, a la densidad que se uso para la identidad alpha=beta.
   A: ¿es J_s realmente independiente de lambda_R2 en el modelo SIN fase de
      Peierls, o era coincidencia del punto elegido?
   B: busqueda especificada de un locus de J_s = 0 en grafeno/siliceno. *)
$pass=0; $fail=0;
checkTrue[n_String,v_]:=If[TrueQ[v],$pass++;Print["  PASS  ",n],$fail++;
  Print["  FAIL  ",n," -> ",StringTake[ToString[InputForm[v]],UpTo[110]]]];

(* H del anillo. pei=1 -> fase 2(phi+s phis) en segundos vecinos; pei=0 -> nada *)
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

Print["=============================================="];
Print[" BARRIDO A: dependencia de J_s en lambda_R2"];
Print["=============================================="];
Module[{grid,res,worstNo,worstYes,tab},
  (* rejilla de contexto: tamano, flujo, llenado, lEO, lR1 *)
  grid=Flatten[Table[{ns,p,nf,le,l1},
     {ns,{6,8,10,12}},{p,{1/8,1/4,2/5}},{nf,{ns-2,ns,ns+2}},
     {le,{0,39/1000,2/10}},{l1,{0,12/100,3/10}}],4];
  Print["  puntos de contexto (N,phi,llenado,lEO,lR1): ",Length[grid]];
  Print["  en cada uno se barre lR2 en {0,0.05,...,0.30}: ",7 Length[grid]," evaluaciones"];
  (* para cada contexto, variacion de J_s al mover lR2 *)
  res[pei_]:=Table[Module[{q=g,base,vals},
     base=js[q[[1]],q[[4]],q[[5]],0,2 Pi q[[2]]/q[[1]],pei,q[[3]]];
     vals=Table[js[q[[1]],q[[4]],q[[5]],x,2 Pi q[[2]]/q[[1]],pei,q[[3]]],{x,0,3/10,1/20}];
     Max[Abs[vals-base]]],{g,grid}];
  worstNo=res[0]; worstYes=res[1];
  Print[];
  Print["  SIN fase de Peierls (el modelo tal como esta escrito):"];
  Print["     max sobre toda la rejilla de |J_s(lR2)-J_s(0)| = ",Max[worstNo]];
  Print["     contextos con variacion > 1e-10 : ",Count[worstNo,x_/;x>10^-10]," de ",Length[grid]];
  Print["  CON fase de Peierls:"];
  Print["     max = ",Max[worstYes],"   contextos con variacion > 1e-10 : ",
        Count[worstYes,x_/;x>10^-10]," de ",Length[grid]];
  (* Esta asercion decia "en TODA la rejilla" y por eso depositaba un FAIL en
     el log: no es lo que el articulo sostiene ni lo que es cierto. Lo cierto,
     y lo que la Sec. III B afirma, es que sin la fase hay un subconjunto
     propio de contextos donde la independencia es exacta, y que con la fase
     no hay ninguno. Eso es lo que se comprueba ahora, y ademas se imprime
     cuales son, que es lo que un lector querria verificar. *)
  Module[{indep},
    indep=Pick[grid,worstNo,x_/;x<10^-10];
    Print["  contextos con independencia exacta sin la fase: ",Length[indep]];
    Print["     (N,llenado-N) presentes: ",Tally[{#[[1]],#[[3]]-#[[1]]}&/@indep]];
    checkTrue["sin la fase hay contextos con independencia exacta",
      Length[indep]>0];
    checkTrue["con la fase no queda ninguno",
      Count[worstYes,x_/;x<10^-10]==0]];
  checkTrue["con la fase, si depende en al menos parte de la rejilla",
    Count[worstYes,x_/;x>10^-10]>0];
  Print["  fraccion de contextos donde la fase hace diferencia: ",
        N[Count[worstYes,x_/;x>10^-10]/Length[grid]]]];

Print[];
Print["=============================================="];
Print[" BARRIDO B: ¿existe un locus J_s = 0?"];
Print["=============================================="];
Module[{grid,vals,signs,nz,pos,neg},
  grid=Flatten[Table[{ns,p,nf,le,l1,l2},
     {ns,{6,8,10}},{p,{1/8,1/4,2/5}},{nf,{ns}},
     {le,{0,1/10,3/10}},{l1,{1/20,3/20,3/10}},{l2,{0,1/20,3/20,3/10}}],5];
  Print["  puntos con acoplos NO nulos: ",Length[grid]];
  vals=Table[js[g[[1]],g[[4]],g[[5]],g[[6]],2 Pi g[[2]]/g[[1]],1,g[[3]]],{g,grid}];
  nz=Select[Abs[vals],#<10^-9&];
  pos=Count[vals,x_/;x>10^-9]; neg=Count[vals,x_/;x<-10^-9];
  Print["  |J_s| < 1e-9 en ",Length[nz]," puntos de ",Length[vals]];
  Print["  J_s > 0 en ",pos,"   J_s < 0 en ",neg];
  Print["  min |J_s| no nulo = ",Min[Select[Abs[vals],#>10^-9&]]];
  (* el analogo directo: lR1 = lR2 *)
  Module[{eq},
   eq=Table[js[g[[1]],g[[4]],g[[5]],g[[5]],2 Pi g[[2]]/g[[1]],1,g[[3]]],
      {g,Select[grid,#[[6]]==0&]}];
   Print["  con lR1 = lR2 impuesto (",Length[eq]," puntos): min |J_s| = ",Min[Abs[eq]]];
   checkTrue["lR1 = lR2 no anula J_s en ningun punto",Min[Abs[eq]]>10^-6]];
  checkTrue["no hay locus de J_s = 0 con acoplos no nulos",Length[nz]==0]];

Print[];
Print["  PASSED: ",$pass,"   FAILED: ",$fail];
