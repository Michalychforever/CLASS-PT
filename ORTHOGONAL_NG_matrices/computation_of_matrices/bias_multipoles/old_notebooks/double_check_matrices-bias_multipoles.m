(* ::Package:: *)

(* ::Code:: *)
(*ClearAll[Evaluate[Context[]<>"*"]];*)
(*ClearSystemCache[];*)


(* ::Input:: *)
(*\[ScriptCapitalI][n1_,n2_]:=(Gamma[3/2-n1] Gamma[3/2-n2] Gamma[-(3/2)+n1+n2])/(8 \[Pi]^(3/2) Gamma[n1] Gamma[3-n1-n2] Gamma[n2]);*)


(* ::Input:: *)
(*J2[k_,\[Nu]1_,\[Nu]2_]:=1/(2\[Pi])^3 (Gamma[3/2-\[Nu]1]Gamma[3/2-\[Nu]2]Gamma[\[Nu]1+\[Nu]2-3/2])/(Gamma[\[Nu]1]Gamma[\[Nu]2]Gamma[3-\[Nu]1-\[Nu]2]) \[Pi]^(3/2) k^(3-2\[Nu]1-2\[Nu]2);*)
(*A1[\[Nu]1_,\[Nu]2_]:=1/2 (J2[1,\[Nu]1-1,\[Nu]2]-J2[1,\[Nu]1,\[Nu]2-1]+J2[1,\[Nu]1,\[Nu]2]);*)
(*A2[\[Nu]1_,\[Nu]2_]:=1/8 (-J2[1,-2+\[Nu]1,\[Nu]2]+2 J2[1,-1+\[Nu]1,-1+\[Nu]2]+2 J2[1,-1+\[Nu]1,\[Nu]2]-J2[1,\[Nu]1,-2+\[Nu]2]+2 J2[1,\[Nu]1,-1+\[Nu]2]-J2[1,\[Nu]1,\[Nu]2]);*)
(*B2[\[Nu]1_,\[Nu]2_]:=1/8 (3 J2[1,-2+\[Nu]1,\[Nu]2]-6 J2[1,-1+\[Nu]1,-1+\[Nu]2]+2 J2[1,-1+\[Nu]1,\[Nu]2]+3 J2[1,\[Nu]1,-2+\[Nu]2]-6 J2[1,\[Nu]1,-1+\[Nu]2]+3 J2[1,\[Nu]1,\[Nu]2]);*)


(* ::Input:: *)
(*(* NOW THE BASELINE WILL JUST BE \[ScriptCapitalI][n1,n2] \[Rule] I put everything in the python script... *)*)


(* ::DisplayFormula:: *)
(*Subscript[Z, 2]\[Superset] Subscript[b, 2]/2*Subscript[f, NL] -> k^3 as the dimension, \[Nu]=-1.25: notice that the factor of Subscript[b, 2]/2 is outside . I . e . I do not have an overall 1/2...*)
(*\[Mu]^0: Pnewabb2*)


(* ::Input:: *)
(*PNEWabb2[n1_,n2_]:=((Csc[n1 \[Pi]] Csc[n2 \[Pi]] Csc[(n1+n2) \[Pi]] ((-15-8 n1 (-1+n2)-8 (-3+n2) n2) Sin[2 n1 \[Pi]]+(-15+8 n2-8 n1 (-3+n1+n2)) Sin[2 n2 \[Pi]]+(-9+8 n1+8 n2-8 n1 n2) Sin[2 (n1+n2) \[Pi]]))/((-1+2 n1) (-1+2 n2) (-5+2 n1+2 n2)))(*\[ScriptCapitalI][n1,n2]*);*)


(* ::DisplayFormula:: *)
(*Subscript[Z, 2]\[Superset] Subscript[b, Subscript[\[ScriptCapitalG], 2]]*Subscript[f, NL]*Subscript[\[ScriptCapitalG], 2][Subscript[k, 1],Subscript[k, 2]] -> k^3 as the dimension, \[Nu]=-0.8 (notice that it is different from the above)*)
(*\[Mu]^0: Pnewab\[ScriptCapitalG]2*)


(* ::Input:: *)
(*PNEWab\[ScriptCapitalG]2[n1_,n2_]:=((Csc[n1 \[Pi]] Csc[n2 \[Pi]] ((-30+96 n2+16 (-3 (-2+n1) n1+n1 (-43+2 (25-8 n1) n1) n2+(-3+2 n1 (-5+2 n1) (-5+4 n1)) n2^2+16 (-1+n1) n1 n2^3)) Cos[(n1+n2) \[Pi]]+Csc[(n1+n2) \[Pi]] ((15+128 n1^3 (-1+n2) n2-8 (-2+n2) n2 (1+4 (-2+n2) n2)+8 n1 (-3+2 n2) (2+n2 (11+8 (-3+n2) n2))+8 n1^2 (3+2 n2 (23+2 n2 (-21+8 n2)))) Sin[2 n1 \[Pi]]+15 Sin[2 n2 \[Pi]]+8 (-(-2+n1) n1 (1+4 (-2+n1) n1)+(-3+2 n1) (2+n1 (11+8 (-3+n1) n1)) n2+(3+2 n1 (23+2 n1 (-21+8 n1))) n2^2+16 (-1+n1) n1 n2^3) Sin[2 n2 \[Pi]])))/(2 n1 (-1+4 n1^2) n2 (-5+2 n1+2 n2) (-1+4 n2^2)))(*\[ScriptCapitalI][n1,n2]*);*)


(* ::Input:: *)
(*b2mattest=(Csc[nu1 \[Pi]] Csc[nu2 \[Pi]] Csc[(nu1+nu2) \[Pi]] ((-15-8 nu1 (-1+nu2)-8 (-3+nu2) nu2) Sin[2 nu1 \[Pi]]+(-15+8 nu2-8 nu1 (-3+nu1+nu2)) Sin[2 nu2 \[Pi]]+(-9+8 nu1+8 nu2-8 nu1 nu2) Sin[2 (nu1+nu2) \[Pi]]))/((-1+2 nu1) (-1+2 nu2) (-5+2 nu1+2 nu2))*)


(* ::Input:: *)
(*\[ScriptCapitalG]2mattest=(Csc[nu1 \[Pi]] Csc[nu2 \[Pi]] ((-30+96 nu2+16 (-3 (-2+nu1) nu1+nu1 (-43+2 (25-8 nu1) nu1) nu2+(-3+2 nu1 (-5+2 nu1) (-5+4 nu1)) nu2^2+16 (-1+nu1) nu1 nu2^3)) Cos[(nu1+nu2) \[Pi]]+Csc[(nu1+nu2) \[Pi]] ((15+128 nu1^3 (-1+nu2) nu2-8 (-2+nu2) nu2 (1+4 (-2+nu2) nu2)+8 nu1 (-3+2 nu2) (2+nu2 (11+8 (-3+nu2) nu2))+8 nu1^2 (3+2 nu2 (23+2 nu2 (-21+8 nu2)))) Sin[2 nu1 \[Pi]]+15 Sin[2 nu2 \[Pi]]+8 (-((-2+nu1) nu1 (1+4 (-2+nu1) nu1))+(-3+2 nu1) (2+nu1 (11+8 (-3+nu1) nu1)) nu2+(3+2 nu1 (23+2 nu1 (-21+8 nu1))) nu2^2+16 (-1+nu1) nu1 nu2^3) Sin[2 nu2 \[Pi]])))/(2 nu1 (-1+4 nu1^2) nu2 (-5+2 nu1+2 nu2) (-1+4 nu2^2))*)


(* ::Input:: *)
(*Simplify[PNEWabb2[nu1,nu2]/b2mattest]*)
(*Simplify[PNEWab\[ScriptCapitalG]2[nu1,nu2]/\[ScriptCapitalG]2mattest]*)


(* ::Title:: *)
(*Multipoles*)


(* ::Input:: *)
(*b2preZ1=PNEWabb2[nu1,nu2];*)
(*b2=(b1+f*\[Mu]^2)*b2preZ1;*)


(* ::Input:: *)
(*\[ScriptCapitalG]2preZ1=PNEWab\[ScriptCapitalG]2[nu1,nu2];*)
(*\[ScriptCapitalG]2=(b1+f*\[Mu]^2)*\[ScriptCapitalG]2preZ1;*)


(* ::Input:: *)
(*rulesb2=CoefficientRules[b2,\[Mu]];*)
(*rules\[ScriptCapitalG]2=CoefficientRules[\[ScriptCapitalG]2,\[Mu]];*)


(* ::Input:: *)
(*Do[CCb2[j]={j}/.rulesb2;,{j,{0,2,4,6}}];*)
(*Do[CC\[ScriptCapitalG]2[j]={j}/.rules\[ScriptCapitalG]2;,{j,{0,2,4,6}}];*)


(* ::Input:: *)
(*\[Mu]b2coeffs={CCb2[0],CCb2[2],CCb2[4],CCb2[6]};*)
(*\[Mu]\[ScriptCapitalG]2coeffs={CC\[ScriptCapitalG]2[0],CC\[ScriptCapitalG]2[2],CC\[ScriptCapitalG]2[4],CC\[ScriptCapitalG]2[6]};*)


(* ::Input:: *)
(*CCb2[4]*)
(*CCb2[6]*)


(* ::Input:: *)
(*CC\[ScriptCapitalG]2[4]*)
(*CC\[ScriptCapitalG]2[6]*)


(* ::Input:: *)
(*\[Mu]b2coeffs=Simplify[{CCb2[0],CCb2[2]}]*)
(*\[Mu]\[ScriptCapitalG]2coeffs=Simplify[{CC\[ScriptCapitalG]2[0],CC\[ScriptCapitalG]2[2]}]*)


(* ::Input:: *)
(*(* I had a notebook and corresponding .c file for matter only, which is 2Subscript[F, 2]\[Times]shape. Then, the ones WITHOUT the name mu_powers are for the MULTIPOLES before AP, and those with mu_powers are for the \[Mu] powers when I need to get AP... *)*)


(* ::Program:: *)
(*coeffs[0]={1,1/3,1/5};*)
(*coeffs[2]={0,2/3,4/7};*)
(*coeffs[4]={0,0,8/35};*)


(* ::Input:: *)
(*(* need to sum the coefficients of \[Mu]^(2n) times these coefficients to get monopole, quadrupole, etc. Notice that here I copied from a notebook where I forgot Subscript[Z, 1], so I didn't have \[Mu]^6... Since here I go up to \[Mu]^2 only things are easier. *)*)


(* ::Input:: *)
(*(2*0+1)/2*Integrate[LegendreP[0,x]*(A+B*x^2),{x,-1,1}]//Expand*)


(* ::Input:: *)
(*(2*2+1)/2*Integrate[LegendreP[2,x]*(A+B*x^2),{x,-1,1}]//Expand*)


(* ::Input:: *)
(*coeffs[0]={1,1/3};*)
(*coeffs[2]={0,2/3};*)


(* ::Input:: *)
(*total={};*)


(* ::Chapter:: *)
(*Subscript[b, 2]*)


(* ::Input:: *)
(*rulesb2b1\[Mu]0=CoefficientRules[\[Mu]b2coeffs . coeffs[0],b1];*)
(*rulesb2b1\[Mu]2=CoefficientRules[\[Mu]b2coeffs . coeffs[2],b1];*)


(* ::Subtitle:: *)
(*Monopole*)


(* ::Subsection:: *)
(*(* b1^0 *)*)


(* ::Input:: *)
(*vv0=Simplify[{0}/.rulesb2b1\[Mu]0]*)


(* ::Input:: *)
(*\[ScriptF]rules=CoefficientRules[vv0,f]//Simplify;*)


(* ::Input:: *)
(*whichpowers=\[ScriptF]rules[[All,1]]//Flatten//Sort*)


(* ::Input:: *)
(*Do[vv0\[ScriptF][iii][nu1_,nu2_]=Simplify[({iii}/.\[ScriptF]rules)];,{iii,whichpowers}];*)


(* ::Input:: *)
(*Do[AppendTo[total,{"b2, vv0, f = "<>ToString[iii],vv0\[ScriptF][iii][nu1,nu2]}],{iii,whichpowers}];*)


(* ::Program:: *)
(*Do[Export[NotebookDirectory[] <> "b2_vv0_f" <> ToString[iii] <> ".txt", vv0\[ScriptF][iii][nu1, nu2]];, {iii, whichpowers}]*)


(* ::Input:: *)
(*Print[total]*)


(* ::Subsection:: *)
(*(* b1^1 *)*)


(* ::Input:: *)
(*vd0=Simplify[{1}/.rulesb2b1\[Mu]0]*)


(* ::Input:: *)
(*\[ScriptF]rules=CoefficientRules[vd0,f]//Simplify;*)


(* ::Input:: *)
(*whichpowers=\[ScriptF]rules[[All,1]]//Flatten//Sort*)


(* ::Input:: *)
(*Do[vd0\[ScriptF][iii][nu1_,nu2_]=Simplify[({iii}/.\[ScriptF]rules)];,{iii,whichpowers}];*)


(* ::Input:: *)
(*Do[AppendTo[total,{"b2, vd0, f = "<>ToString[iii],vd0\[ScriptF][iii][nu1,nu2]}],{iii,whichpowers}];*)


(* ::Program:: *)
(*Do[Export[NotebookDirectory[] <> "b2_vd0_f" <> ToString[iii] <> ".txt", vd0\[ScriptF][iii][nu1, nu2]];, {iii, whichpowers}]*)


(* ::Input:: *)
(*Print[total]*)


(* ::Subsection:: *)
(*(* b1^2 *)*)


(* ::Input:: *)
(*dd0=Simplify[{2}/.rulesb2b1\[Mu]0]*)


(* ::Subtitle:: *)
(*Quadrupole*)


(* ::Subsection:: *)
(*(* b1^0 *)*)


(* ::Input:: *)
(*vv2=Simplify[{0}/.rulesb2b1\[Mu]2]*)


(* ::Input:: *)
(*\[ScriptF]rules=CoefficientRules[vv2,f]//Simplify;*)


(* ::Input:: *)
(*whichpowers=\[ScriptF]rules[[All,1]]//Flatten//Sort*)


(* ::Input:: *)
(*Do[vv2\[ScriptF][iii][nu1_,nu2_]=Simplify[({iii}/.\[ScriptF]rules)];,{iii,whichpowers}];*)


(* ::Input:: *)
(*Do[AppendTo[total,{"b2, vv2, f = "<>ToString[iii],vv2\[ScriptF][iii][nu1,nu2]}],{iii,whichpowers}];*)


(* ::Program:: *)
(*Do[Export[NotebookDirectory[] <> "b2_vv2_f" <> ToString[iii] <> ".txt", vv2\[ScriptF][iii][nu1, nu2]];, {iii, whichpowers}]*)


(* ::Input:: *)
(*Print[total]*)


(* ::Subsection:: *)
(*(* b1^1 *)*)


(* ::Input:: *)
(*vd2=Simplify[{1}/.rulesb2b1\[Mu]2]*)


(* ::Subsection:: *)
(*(* b1^2 *)*)


(* ::Program:: *)
(*dd2 = Simplify[{2} /. rulesb2b1\[Mu]0]*)


(* ::Program:: *)
(*{2}*)


(* ::Input:: *)
(*(* LUCKY!!! *)*)


(* ::Input:: *)
(*dd2=Simplify[{2}/.rulesb2b1\[Mu]2]*)


(* ::Input:: *)
(*(* NOTICE THAT FOR SOME REASON SOME OF THE NUMERICAL FACTORS IN FRONT WERE NOT DE-INTEGERIZED... DID IT DIRECTLY ON C FOR THE OLD MATRICES... *)*)


(* ::Input:: *)
(*length=Length[total]*)


(* ::Input:: *)
(*Monitor[Do[If[ii<jj&&RationalExpressionQ[Simplify[total[[ii]][[2]]/total[[jj]][[2]]],{nu1,nu2}],Print[{ii,jj}]],{ii,1,length},{jj,1,length}],{ii,jj}]*)


(* ::Input:: *)
(*(* actually it pays to do it... We see that we just need to export -> *)*)


(* ::Input:: *)
(*total[[2]][[2]]/total[[1]][[2]]*)


(* ::Input:: *)
(*total[[2]][[1]]*)
(*total[[1]][[1]]*)


(* ::Input:: *)
(*(*############################*)*)


(* ::Input:: *)
(*total[[3]][[2]]/total[[1]][[2]]*)


(* ::Input:: *)
(*total[[3]][[1]]*)
(*total[[1]][[1]]*)


(* ::Chapter:: *)
(*Subscript[\[ScriptCapitalG], 2]*)


(* ::Input:: *)
(*rules\[ScriptCapitalG]2b1\[Mu]0=CoefficientRules[\[Mu]\[ScriptCapitalG]2coeffs . coeffs[0],b1];*)
(*rules\[ScriptCapitalG]2b1\[Mu]2=CoefficientRules[\[Mu]\[ScriptCapitalG]2coeffs . coeffs[2],b1];*)


(* ::Subtitle:: *)
(*Monopole*)


(* ::Input:: *)
(*total={};*)


(* ::Subsection:: *)
(*(* b1^0 *)*)


(* ::Input:: *)
(*vv0=Simplify[{0}/.rules\[ScriptCapitalG]2b1\[Mu]0]*)


(* ::Input:: *)
(*\[ScriptF]rules=CoefficientRules[vv0,f]//Simplify;*)


(* ::Input:: *)
(*whichpowers=\[ScriptF]rules[[All,1]]//Flatten//Sort*)


(* ::Input:: *)
(*Do[vv0\[ScriptF][iii][nu1_,nu2_]=Simplify[({iii}/.\[ScriptF]rules)];,{iii,whichpowers}];*)


(* ::Input:: *)
(*Do[AppendTo[total,{"bG2, vv0, f = "<>ToString[iii],vv0\[ScriptF][iii][nu1,nu2]}],{iii,whichpowers}];*)


(* ::Program:: *)
(*Do[Export[NotebookDirectory[] <> "bG2_vv0_f" <> ToString[iii] <> ".txt", vv0\[ScriptF][iii][nu1, nu2]];, {iii, whichpowers}]*)


(* ::Input:: *)
(*Print[total]*)


(* ::Subsection:: *)
(*(* b1^1 *)*)


(* ::Input:: *)
(*vd0=Simplify[{1}/.rules\[ScriptCapitalG]2b1\[Mu]0]*)


(* ::Input:: *)
(*\[ScriptF]rules=CoefficientRules[vd0,f]//Simplify;*)


(* ::Input:: *)
(*whichpowers=\[ScriptF]rules[[All,1]]//Flatten//Sort*)


(* ::Input:: *)
(*Do[vd0\[ScriptF][iii][nu1_,nu2_]=Simplify[({iii}/.\[ScriptF]rules)];,{iii,whichpowers}];*)


(* ::Input:: *)
(*Do[AppendTo[total,{"bG2, vd0, f = "<>ToString[iii],vd0\[ScriptF][iii][nu1,nu2]}],{iii,whichpowers}];*)


(* ::Program:: *)
(*Do[Export[NotebookDirectory[] <> "bG2_vd0_f" <> ToString[iii] <> ".txt", vd0\[ScriptF][iii][nu1, nu2]];, {iii, whichpowers}]*)


(* ::Input:: *)
(*Print[total]*)


(* ::Subsection:: *)
(*(* b1^2 *)*)


(* ::Input:: *)
(*dd0=Simplify[{2}/.rules\[ScriptCapitalG]2b1\[Mu]0]*)


(* ::Subtitle:: *)
(*Quadrupole*)


(* ::Subsection:: *)
(*(* b1^0 *)*)


(* ::Input:: *)
(*vv2=Simplify[{0}/.rules\[ScriptCapitalG]2b1\[Mu]2]*)


(* ::Input:: *)
(*\[ScriptF]rules=CoefficientRules[vv2,f]//Simplify;*)


(* ::Input:: *)
(*whichpowers=\[ScriptF]rules[[All,1]]//Flatten//Sort*)


(* ::Input:: *)
(*Do[vv2\[ScriptF][iii][nu1_,nu2_]=Simplify[({iii}/.\[ScriptF]rules)];,{iii,whichpowers}];*)


(* ::Input:: *)
(*Do[AppendTo[total,{"bG2, vv2, f = "<>ToString[iii],vv2\[ScriptF][iii][nu1,nu2]}],{iii,whichpowers}];*)


(* ::Program:: *)
(*Do[Export[NotebookDirectory[] <> "bG2_vv2_f" <> ToString[iii] <> ".txt", vv2\[ScriptF][iii][nu1, nu2]];, {iii, whichpowers}]*)


(* ::Input:: *)
(*Print[total]*)


(* ::Subsection:: *)
(*(* b1^1 *)*)


(* ::Input:: *)
(*vd2=Simplify[{1}/.rules\[ScriptCapitalG]2b1\[Mu]2]*)


(* ::Subsection:: *)
(*(* b1^2 *)*)


(* ::Input:: *)
(*dd2=Simplify[{2}/.rules\[ScriptCapitalG]2b1\[Mu]2]*)


(* ::Input:: *)
(*length=Length[total]*)


(* ::Input:: *)
(*Monitor[Do[If[ii<jj&&RationalExpressionQ[Simplify[total[[ii]][[2]]/total[[jj]][[2]]],{nu1,nu2}],Print[{ii,jj}]],{ii,1,length},{jj,1,length}],{ii,jj}]*)


(* ::Input:: *)
(*(* actually it pays to do it... We see that we just need to export -> *)*)


(* ::Input:: *)
(*total[[2]][[2]]/total[[1]][[2]]*)


(* ::Input:: *)
(*total[[2]][[1]]*)
(*total[[1]][[1]]*)


(* ::Input:: *)
(*(*############################*)*)


(* ::Input:: *)
(*total[[3]][[2]]/total[[1]][[2]]*)


(* ::Input:: *)
(*total[[3]][[1]]*)
(*total[[1]][[1]]*)



