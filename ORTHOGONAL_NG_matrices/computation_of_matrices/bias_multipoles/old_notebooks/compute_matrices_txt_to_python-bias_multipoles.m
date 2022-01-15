(* ::Package:: *)

(* ::Code:: *)
(*ClearAll[Evaluate[Context[]<>"*"]];*)
(*ClearSystemCache[];*)


(* ::Input:: *)
(*(* as detailed in the other notebook, I had screwed up by not actually multiplying by \[ScriptF] when extracting the multipoles? Ah no, actually ... Now \[Rule] BEFORE MULTIPLYING BY Subscript[Z, 1], which will need to be changed if we include FoG, I have (recall that the division by k^2 is done at the level of \[ScriptCapitalM] [THIS IS FUNDAMENTAL, I do not have a simple way to avoid this...], and "NGparam" will be overall, once I decide 2/5 vs. 2/3... Overall there is also \[ScriptCapitalM][k]/k^2 FOR EVERYONE, where I recall all k are in units of h/Mpc...) => *)*)


(* ::Input:: *)
(*(* for the bias operators, I only care about the \[Mu] powers, not multipoles. Also, the \[Mu] powers are trivial since they come from Subscript[Z, 1]: so I will have a power of Subscript[b, 1], and one of f\[Mu]^2. NOTHING MORE... I do not need to keep this into the matrices, I can add it by myself... *)*)


(* ::Input:: *)
(*\[ScriptCapitalI][n1_,n2_]:=(Gamma[3/2-n1] Gamma[3/2-n2] Gamma[-(3/2)+n1+n2])/(8 \[Pi]^(3/2) Gamma[n1] Gamma[3-n1-n2] Gamma[n2]);*)


(* ::Input:: *)
(*J2[k_,\[Nu]1_,\[Nu]2_]:=1/(2\[Pi])^3 (Gamma[3/2-\[Nu]1]Gamma[3/2-\[Nu]2]Gamma[\[Nu]1+\[Nu]2-3/2])/(Gamma[\[Nu]1]Gamma[\[Nu]2]Gamma[3-\[Nu]1-\[Nu]2]) \[Pi]^(3/2) k^(3-2\[Nu]1-2\[Nu]2);*)
(*A1[\[Nu]1_,\[Nu]2_]:=1/2 (J2[1,\[Nu]1-1,\[Nu]2]-J2[1,\[Nu]1,\[Nu]2-1]+J2[1,\[Nu]1,\[Nu]2]);*)
(*A2[\[Nu]1_,\[Nu]2_]:=1/8 (-J2[1,-2+\[Nu]1,\[Nu]2]+2 J2[1,-1+\[Nu]1,-1+\[Nu]2]+2 J2[1,-1+\[Nu]1,\[Nu]2]-J2[1,\[Nu]1,-2+\[Nu]2]+2 J2[1,\[Nu]1,-1+\[Nu]2]-J2[1,\[Nu]1,\[Nu]2]);*)
(*B2[\[Nu]1_,\[Nu]2_]:=1/8 (3 J2[1,-2+\[Nu]1,\[Nu]2]-6 J2[1,-1+\[Nu]1,-1+\[Nu]2]+2 J2[1,-1+\[Nu]1,\[Nu]2]+3 J2[1,\[Nu]1,-2+\[Nu]2]-6 J2[1,\[Nu]1,-1+\[Nu]2]+3 J2[1,\[Nu]1,\[Nu]2]);*)


(* ::Input:: *)
(*PNEWab[n1_,n2_]:=(Csc[n1 \[Pi]] Csc[n2 \[Pi]] (-2 (60+448 n1^4 (-1+n2) n2-3 n2 (43+4 n2 (-15+7 n2))+4 n1^3 (-21+2 n2 (169+4 n2 (-65+28 n2)))+4 n1^2 (45+4 n2 (-79+n2 (169+2 n2 (-65+14 n2))))+n1 (-129-8 n2 (-60+n2 (158+n2 (-169+56 n2))))) Cos[(n1+n2) \[Pi]]-Csc[(n1+n2) \[Pi]] ((-60+3 n1 (43+4 n1 (-15+7 n1))+209 n2-8 n1 (17+(-2+n1) n1 (-43+56 n1)) n2+4 (101+2 n1 (-229+2 n1 (245+4 n1 (-43+7 n1)))) n2^2+4 (-345+2 n1 (463+4 n1 (-137+42 n1))) n2^3+64 (16+3 n1 (-12+7 n1)) n2^4+224 (-1+2 n1) n2^5) Sin[2 n1 \[Pi]]+(-60+n1 (209+4 n1 (101+n1 (-345+8 (32-7 n1) n1)))+129 n2+8 n1 (-17+n1 (-229+n1 (463+8 n1 (-36+7 n1)))) n2+4 (-45+4 n1 (-43+n1 (245-274 n1+84 n1^2))) n2^2+4 (21+2 n1 (155+8 n1 (-43+21 n1))) n2^3+448 (-1+n1) n1 n2^4) Sin[2 n2 \[Pi]])))/(28 n1 (-1+4 n1^2) n2 (-5+2 n1+2 n2) (-1+4 n2^2));*)


(* ::Program:: *)
(*((Csc[n1 \[Pi]] Csc[n2 \[Pi]] (-2 (60+448 n1^4 (-1+n2) n2-3 n2 (43+4 n2 (-15+7 n2))+4 n1^3 (-21+2 n2 (169+4 n2 (-65+28 n2)))+4 n1^2 (45+4 n2 (-79+n2 (169+2 n2 (-65+14 n2))))+n1 (-129-8 n2 (-60+n2 (158+n2 (-169+56 n2))))) Cos[(n1+n2) \[Pi]]-Csc[(n1+n2) \[Pi]] ((-60+3 n1 (43+4 n1 (-15+7 n1))+209 n2-8 n1 (17+(-2+n1) n1 (-43+56 n1)) n2+4 (101+2 n1 (-229+2 n1 (245+4 n1 (-43+7 n1)))) n2^2+4 (-345+2 n1 (463+4 n1 (-137+42 n1))) n2^3+64 (16+3 n1 (-12+7 n1)) n2^4+224 (-1+2 n1) n2^5) Sin[2 n1 \[Pi]]+(-60+n1 (209+4 n1 (101+n1 (-345+8 (32-7 n1) n1)))+129 n2+8 n1 (-17+n1 (-229+n1 (463+8 n1 (-36+7 n1)))) n2+4 (-45+4 n1 (-43+n1 (245-274 n1+84 n1^2))) n2^2+4 (21+2 n1 (155+8 n1 (-43+21 n1))) n2^3+448 (-1+n1) n1 n2^4) Sin[2 n2 \[Pi]])))/(28 n1 (-1+4 n1^2) n2 (-5+2 n1+2 n2) (-1+4 n2^2)))*)
(*(* \[Rule] Subscript[b, 1]Subscript[F, 2] piece from full_multipoles.nb. Exactly as expected it is the same... *)*)


(* ::Input:: *)
(*(* THIS IS FOR 2Subscript[F, 2]^2: CAREFUL (NO, it is that for Misha [see CHECK_M22_PYTHON.nb \[Rule] I checked that indeed M22oneline.dat WITHOUT BASIC is exactly the term from 2Subscript[F, 2]^2... And maybe this is what I had in mind when I wrote this!]... For me it is the matter-only piece, the 2Subscript[F, 2]\[Times]bispectrum, with the 2 of the permutations. FUCK ALL THESE TYPOS IN COMMENTS \[Rule] I think the issue was just that I wanted to say 2Subscript[F, 2]\[Times]shape as in the C code)!!! This is the dd contribution at \[Mu]=0 for matter, and indeed b1\[Mu]powers[1,0] is 1 now. This will be the baseline... *)*)


(* ::Input:: *)
(*(*matrix=PNEWab[nu1,nu2]*)*)


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
(*(* I also need to look at the multipoles, before the AP... *)*)


(* ::Input:: *)
(*(* if I look only at matter and multiply by Subscript[Z, 1] \[Rule] b1\[Mu]powers => everything is multiplied by Subscript[f, NL], k^3, and the NGparam, and \[ScriptCapitalM][k]/k^2... *)*)


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


(* ::Input:: *)
(*Do[Export[NotebookDirectory[]<>"b2_vv0_f"<>ToString[iii]<>".txt",vv0\[ScriptF][iii][nu1,nu2]];,{iii,whichpowers}]*)


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


(* ::Input:: *)
(*Do[Export[NotebookDirectory[]<>"b2_vd0_f"<>ToString[iii]<>".txt",vd0\[ScriptF][iii][nu1,nu2]];,{iii,whichpowers}]*)


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


(* ::Input:: *)
(*Do[Export[NotebookDirectory[]<>"b2_vv2_f"<>ToString[iii]<>".txt",vv2\[ScriptF][iii][nu1,nu2]];,{iii,whichpowers}]*)


(* ::Subsection:: *)
(*(* b1^1 *)*)


(* ::Input:: *)
(*vd2=Simplify[{1}/.rulesb2b1\[Mu]2]*)


(* ::Subsection:: *)
(*(* b1^2 *)*)


(* ::Input:: *)
(*dd2=Simplify[{2}/.rulesb2b1\[Mu]0]*)


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


(* ::Input:: *)
(*Do[Export[NotebookDirectory[]<>"bG2_vv0_f"<>ToString[iii]<>".txt",vv0\[ScriptF][iii][nu1,nu2]];,{iii,whichpowers}]*)


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


(* ::Input:: *)
(*Do[Export[NotebookDirectory[]<>"bG2_vd0_f"<>ToString[iii]<>".txt",vd0\[ScriptF][iii][nu1,nu2]];,{iii,whichpowers}]*)


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


(* ::Input:: *)
(*Do[Export[NotebookDirectory[]<>"bG2_vv2_f"<>ToString[iii]<>".txt",vv2\[ScriptF][iii][nu1,nu2]];,{iii,whichpowers}]*)


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



