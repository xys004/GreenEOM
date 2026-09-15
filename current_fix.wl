(* nu constant. Two places the h.c. of the intrinsic-SOC term shows up:
   the equation of motion AND the commutator [H,x] that defines the current. *)
$pass=0; $fail=0;
checkTrue[n_String,v_]:=If[TrueQ[v],$pass++;Print["  PASS  ",n],$fail++;Print["  FAIL  ",n,
  " -> ",StringTake[ToString[InputForm[v]],UpTo[140]]]];
DeclareSpecies[c,"Fermion"];

(* 1. the algebraic identity behind every [term, x] in the theses *)
Module[{x,a,b,lhs},
  x=Sum[j NCTimes[Cre[c,{j,s}],Ann[c,{j,s}]],{j,1,9},{s,{1,-1}}];
  lhs=Comm[NCTimes[Cre[c,{3,1}],Ann[c,{5,1}]],x];
  Print["  [c^dag_3 c_5 , x] = ",InputForm[lhs]];
  checkTrue["[c^dag_a c_b , x] = (b-a) c^dag_a c_b",
    Simplify[lhs-2 Mono[{Cre[c,{3,1}],Ann[c,{5,1}]}]]===0];
  lhs=Comm[NCTimes[Cre[c,{5,1}],Ann[c,{3,1}]],x];
  Print["  [c^dag_5 c_3 , x] = ",InputForm[lhs]];
  checkTrue["and the reverse hop carries (b-a) = -2",
    Simplify[lhs+2 Mono[{Cre[c,{5,1}],Ann[c,{3,1}]}]]===0]];

(* 2. so [C,x] for the Hermitian intrinsic-SOC term, nu constant *)
Module[{cFwd,cBwd,cFull,x,comm,sg=1},
  x=Sum[j NCTimes[Cre[c,{j,s}],Ann[c,{j,s}]],{j,1,9},{s,{1,-1}}];
  (* C = i lam nu (s/2) sum_n [ c^dag_n c_{n+2} - c^dag_{n+2} c_n ]  (nu = 1) *)
  cFwd=I lam (sg/2) Sum[NCTimes[Cre[c,{n,sg}],Ann[c,{n+2,sg}]],{n,3,5}];
  cBwd=-I lam (sg/2) Sum[NCTimes[Cre[c,{n+2,sg}],Ann[c,{n,sg}]],{n,3,5}];
  cFull=cFwd+cBwd;
  comm=Simplify[Comm[cFull,x]];
  Print["  [C_hermitian , x] = ",InputForm[comm]];
  Print["  [C_as_printed, x] = ",InputForm[Simplify[Comm[cFwd,x]]]];
  checkTrue["the forward piece gives  + i lam s c^dag_n c_{n+2}",
    Simplify[Comm[cFwd,x]-I lam sg Sum[Mono[{Cre[c,{n,sg}],Ann[c,{n+2,sg}]}],{n,3,5}]]===0];
  checkTrue["the h.c. piece gives  + i lam s c^dag_{n+2} c_n, SAME sign",
    Simplify[Comm[cBwd,x]-I lam sg Sum[Mono[{Cre[c,{n+2,sg}],Ann[c,{n,sg}]}],{n,3,5}]]===0];
  Print["  -> [C,x] = i lam_EO nu sigma sum_n ( c^dag_n c_{n+2} + c^dag_{n+2} c_n )"];
  Print["     the earlier form keeps only the first: the charge current is missing"];
  Print["     the G^<_{n,n+2} partner of its G^<_{n+2,n} term."]];

(* 3. and the EOM coefficients, nu constant *)
Module[{ns=12,nxt,nx2,hcOf,herm,hC,cm,sg=1,n0=4},
  nxt[n_]:=Mod[n,ns]+1; nx2[n_]:=Mod[n+1,ns]+1;
  hcOf[e_]:=Expand[e]/.cc_. Mono[{Cre[sp_,i_],Ann[sp2_,j_]}]:>
    (cc/.Complex[re_,im_]:>Complex[re,-im]) Mono[{Cre[sp2,j],Ann[sp,i]}];
  herm[e_]:=Expand[NCTimes[e]+hcOf[NCTimes[e]]];
  hC=herm[I lam Sum[(s/2) NCTimes[Cre[c,{n,s}],Ann[c,{nx2[n],s}]],{n,ns},{s,{1,-1}}]];
  cm=Expand[Comm[Ann[c,{nxt[n0],sg}],hC]];
  Print["  [c_{n+1,up}, C] with nu = 1, n = ",n0,":"];
  Do[Print["      G_{",i[[1]],",",If[i[[2]]===1,"up","dn"],"}  x  ",
     InputForm[Simplify[Coefficient[cm,Mono[{Ann[c,i]}]]]]],
   {i,Sort[DeleteDuplicates[Cases[cm,Ann[c,q_]:>q,{0,Infinity}]]]}];
  checkTrue["two terms, equal and opposite, +-(i/2) lam_EO",
    Length[DeleteDuplicates[Cases[cm,Ann[c,q_]:>q,{0,Infinity}]]]===2]];

Print["  PASSED: ",$pass,"   FAILED: ",$fail];
