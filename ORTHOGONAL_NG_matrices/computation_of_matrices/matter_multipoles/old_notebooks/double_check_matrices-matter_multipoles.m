(* ::Package:: *)

(* ::Code:: *)
(*ClearAll[Evaluate[Context[]<>"*"]];*)
(*ClearSystemCache[];*)


(* ::Input:: *)
(*(* Subscript[k, min] and Subscript[k, max] enter logarithmically? No, Subscript[k, min] also overall... But NOT in the matrices... Because they just eat Subscript[\[Eta], m] or whatever, that does NOT depend in any way but logarithmic... *)*)


(* ::Program:: *)
(*Nmax = 512;*)
(*b = -8/10;*)


(* ::Program:: *)
(*kmax=100;*)
(*k0=5*10^-5;*)
(*Delta=Log[kmax/k0]/(Nmax-1);*)
(*jsNm=Table[el,{el,-Nmax/2,Nmax/2}];*)


(* ::Program:: *)
(*etam = b + 2*I*\[Pi]*(jsNm)/Nmax/Delta;*)


(* ::Input:: *)
(*(* the issue is that I cannot combine in a single term because the growth rate gets spread around. I need the monopole matrix as an expansion in powers of the growth rate, if I want to read it directly from a table. So there will be many terms needed... *)*)


(* ::Section:: *)
(*Matrices -- 1*)


(* ::Input:: *)
(*(* as detailed in the other notebook, I had screwed up by not actually multiplying by \[ScriptF] when extracting the multipoles? Ah no, actually ... Now \[Rule] BEFORE MULTIPLYING BY Subscript[Z, 1], which will need to be changed if we include FoG, I have (recall that the division by k^2 is done at the level of \[ScriptCapitalM] [THIS IS FUNDAMENTAL, I do not have a simple way to avoid this...], and "NGparam" will be overall, once I decide 2/5 vs. 2/3... Overall there is also \[ScriptCapitalM][k]/k^2 FOR EVERYONE, where I recall all k are in units of h/Mpc...) => *)*)


(* ::Input:: *)
(*\[ScriptCapitalI][n1_,n2_]:=(Gamma[3/2-n1] Gamma[3/2-n2] Gamma[-(3/2)+n1+n2])/(8 \[Pi]^(3/2) Gamma[n1] Gamma[3-n1-n2] Gamma[n2]);*)


(* ::Input:: *)
(*J2[k_,\[Nu]1_,\[Nu]2_]:=1/(2\[Pi])^3 (Gamma[3/2-\[Nu]1]Gamma[3/2-\[Nu]2]Gamma[\[Nu]1+\[Nu]2-3/2])/(Gamma[\[Nu]1]Gamma[\[Nu]2]Gamma[3-\[Nu]1-\[Nu]2]) \[Pi]^(3/2) k^(3-2\[Nu]1-2\[Nu]2);*)
(*A1[\[Nu]1_,\[Nu]2_]:=1/2 (J2[1,\[Nu]1-1,\[Nu]2]-J2[1,\[Nu]1,\[Nu]2-1]+J2[1,\[Nu]1,\[Nu]2]);*)
(*A2[\[Nu]1_,\[Nu]2_]:=1/8 (-J2[1,-2+\[Nu]1,\[Nu]2]+2 J2[1,-1+\[Nu]1,-1+\[Nu]2]+2 J2[1,-1+\[Nu]1,\[Nu]2]-J2[1,\[Nu]1,-2+\[Nu]2]+2 J2[1,\[Nu]1,-1+\[Nu]2]-J2[1,\[Nu]1,\[Nu]2]);*)
(*B2[\[Nu]1_,\[Nu]2_]:=1/8 (3 J2[1,-2+\[Nu]1,\[Nu]2]-6 J2[1,-1+\[Nu]1,-1+\[Nu]2]+2 J2[1,-1+\[Nu]1,\[Nu]2]+3 J2[1,\[Nu]1,-2+\[Nu]2]-6 J2[1,\[Nu]1,-1+\[Nu]2]+3 J2[1,\[Nu]1,\[Nu]2]);*)


(* ::Input:: *)
(*PNEWab[n1_,n2_]:=(Csc[n1 \[Pi]] Csc[n2 \[Pi]] (-2 (60+448 n1^4 (-1+n2) n2-3 n2 (43+4 n2 (-15+7 n2))+4 n1^3 (-21+2 n2 (169+4 n2 (-65+28 n2)))+4 n1^2 (45+4 n2 (-79+n2 (169+2 n2 (-65+14 n2))))+n1 (-129-8 n2 (-60+n2 (158+n2 (-169+56 n2))))) Cos[(n1+n2) \[Pi]]-Csc[(n1+n2) \[Pi]] ((-60+3 n1 (43+4 n1 (-15+7 n1))+209 n2-8 n1 (17+(-2+n1) n1 (-43+56 n1)) n2+4 (101+2 n1 (-229+2 n1 (245+4 n1 (-43+7 n1)))) n2^2+4 (-345+2 n1 (463+4 n1 (-137+42 n1))) n2^3+64 (16+3 n1 (-12+7 n1)) n2^4+224 (-1+2 n1) n2^5) Sin[2 n1 \[Pi]]+(-60+n1 (209+4 n1 (101+n1 (-345+8 (32-7 n1) n1)))+129 n2+8 n1 (-17+n1 (-229+n1 (463+8 n1 (-36+7 n1)))) n2+4 (-45+4 n1 (-43+n1 (245-274 n1+84 n1^2))) n2^2+4 (21+2 n1 (155+8 n1 (-43+21 n1))) n2^3+448 (-1+n1) n1 n2^4) Sin[2 n2 \[Pi]])))/(28 n1 (-1+4 n1^2) n2 (-5+2 n1+2 n2) (-1+4 n2^2));*)


(* ::Input:: *)
(*(* THIS IS FOR 2Subscript[F, 2]^2: CAREFUL!!! This is the dd contribution at \[Mu]=0 for matter, and indeed b1\[Mu]powers[1,0] is 1 now. This will be the baseline... *)*)


(* ::Input:: *)
(*(*matrix=PNEWab[nu1,nu2]*)*)


(* ::DisplayFormula:: *)
(*Subscript[Z, 2]\[Superset] Subscript[b, 2]/2*Subscript[f, NL] -> k^3 as the dimension, \[Nu]=-1.25*)
(*\[Mu]^0: Pnewabb2*)


(* ::Input:: *)
(*PNEWabb2[n1_,n2_]:=((Csc[n1 \[Pi]] Csc[n2 \[Pi]] Csc[(n1+n2) \[Pi]] ((-15-8 n1 (-1+n2)-8 (-3+n2) n2) Sin[2 n1 \[Pi]]+(-15+8 n2-8 n1 (-3+n1+n2)) Sin[2 n2 \[Pi]]+(-9+8 n1+8 n2-8 n1 n2) Sin[2 (n1+n2) \[Pi]]))/((-1+2 n1) (-1+2 n2) (-5+2 n1+2 n2)))(*\[ScriptCapitalI][n1,n2]*);*)


(* ::DisplayFormula:: *)
(*Subscript[Z, 2]\[Superset] Subscript[b, Subscript[\[ScriptCapitalG], 2]]*Subscript[f, NL]*Subscript[\[ScriptCapitalG], 2][Subscript[k, 1],Subscript[k, 2]] -> k^3 as the dimension, \[Nu]=-0.8 (notice that it is different from the above)*)
(*\[Mu]^0: Pnewab\[ScriptCapitalG]2*)


(* ::Input:: *)
(*PNEWab\[ScriptCapitalG]2[n1_,n2_]:=((Csc[n1 \[Pi]] Csc[n2 \[Pi]] ((-30+96 n2+16 (-3 (-2+n1) n1+n1 (-43+2 (25-8 n1) n1) n2+(-3+2 n1 (-5+2 n1) (-5+4 n1)) n2^2+16 (-1+n1) n1 n2^3)) Cos[(n1+n2) \[Pi]]+Csc[(n1+n2) \[Pi]] ((15+128 n1^3 (-1+n2) n2-8 (-2+n2) n2 (1+4 (-2+n2) n2)+8 n1 (-3+2 n2) (2+n2 (11+8 (-3+n2) n2))+8 n1^2 (3+2 n2 (23+2 n2 (-21+8 n2)))) Sin[2 n1 \[Pi]]+15 Sin[2 n2 \[Pi]]+8 (-(-2+n1) n1 (1+4 (-2+n1) n1)+(-3+2 n1) (2+n1 (11+8 (-3+n1) n1)) n2+(3+2 n1 (23+2 n1 (-21+8 n1))) n2^2+16 (-1+n1) n1 n2^3) Sin[2 n2 \[Pi]])))/(2 n1 (-1+4 n1^2) n2 (-5+2 n1+2 n2) (-1+4 n2^2)))(*\[ScriptCapitalI][n1,n2]*);*)


(* ::DisplayFormula:: *)
(*Subscript[Z, 2]\[Superset] Subscript[b, 1]*Subscript[f, NL]*Subscript[F, 2][Subscript[k, 1],Subscript[k, 2]] -> k^3 as the dimension, \[Nu]=-0.8*)
(*\[Mu]^0: PNEWabb1F2*)


(* ::Input:: *)
(*PNEWabb1F2[n1_,n2_]:=((Csc[n1 \[Pi]] Csc[n2 \[Pi]] (-2 (60+448 n1^4 (-1+n2) n2-3 n2 (43+4 n2 (-15+7 n2))+4 n1^3 (-21+2 n2 (169+4 n2 (-65+28 n2)))+4 n1^2 (45+4 n2 (-79+n2 (169+2 n2 (-65+14 n2))))+n1 (-129-8 n2 (-60+n2 (158+n2 (-169+56 n2))))) Cos[(n1+n2) \[Pi]]-Csc[(n1+n2) \[Pi]] ((-60+3 n1 (43+4 n1 (-15+7 n1))+209 n2-8 n1 (17+(-2+n1) n1 (-43+56 n1)) n2+4 (101+2 n1 (-229+2 n1 (245+4 n1 (-43+7 n1)))) n2^2+4 (-345+2 n1 (463+4 n1 (-137+42 n1))) n2^3+64 (16+3 n1 (-12+7 n1)) n2^4+224 (-1+2 n1) n2^5) Sin[2 n1 \[Pi]]+(-60+n1 (209+4 n1 (101+n1 (-345+8 (32-7 n1) n1)))+129 n2+8 n1 (-17+n1 (-229+n1 (463+8 n1 (-36+7 n1)))) n2+4 (-45+4 n1 (-43+n1 (245-274 n1+84 n1^2))) n2^2+4 (21+2 n1 (155+8 n1 (-43+21 n1))) n2^3+448 (-1+n1) n1 n2^4) Sin[2 n2 \[Pi]])))/(28 n1 (-1+4 n1^2) n2 (-5+2 n1+2 n2) (-1+4 n2^2)))(*\[ScriptCapitalI][n1,n2]*);*)


(* ::DisplayFormula:: *)
(*Subscript[Z, 2]\[Superset] Subscript[f, NL]*\[ScriptF]*Subscript[G, 2][Subscript[k, 1],Subscript[k, 2]] -> k^3 as the dimension, \[Nu]=-0.8*)
(*\[Mu]^2: PNEWabG2*)


(* ::Input:: *)
(*PNEWabG2[n1_,n2_]:=((Csc[n1 \[Pi]] Csc[n2 \[Pi]] (2 (-120+3 n1 (107+4 n1 (-23+7 n1))+321 n2+8 (-2+n1) n1 (116+n1 (-121+56 n1)) n2-4 (69+4 n1 (-179+n1 (289+2 n1 (-81+14 n1)))) n2^2+4 (21-2 n1 (233+4 n1 (-81+28 n1))) n2^3-448 (-1+n1) n1 n2^4) Cos[(n1+n2) \[Pi]]-Csc[(n1+n2) \[Pi]] ((-120+448 n1^4 (-1+n2) n2+12 n1^3 (7+2 n2 (73+8 n2 (-17+7 n2)))+n2 (145+4 n2 (237+n2 (-473+8 (36-7 n2) n2)))+4 n1^2 (-69+4 n2 (-135+n2 (413-338 n2+84 n2^2)))+n1 (321+8 n2 (99+n2 (-605+n2 (751+8 n2 (-44+7 n2)))))) Sin[2 n1 \[Pi]]+(-120+n1 (145+4 n1 (237+n1 (-473+8 (36-7 n1) n1)))+321 n2+8 n1 (99+n1 (-605+n1 (751+8 n1 (-44+7 n1)))) n2+4 (-69+4 n1 (-135+n1 (413-338 n1+84 n1^2))) n2^2+12 (7+2 n1 (73+8 n1 (-17+7 n1))) n2^3+448 (-1+n1) n1 n2^4) Sin[2 n2 \[Pi]])))/(28 n1 (-1+4 n1^2) n2 (-5+2 n1+2 n2) (-1+4 n2^2)))(*\[ScriptCapitalI][n1,n2]*);*)


(* ::DisplayFormula:: *)
(*Subscript[Z, 2]\[Superset] Subscript[b, 1]*Subscript[f, NL]*((\[ScriptF]*\[Mu]*k)/2)*(Subscript[\[Mu], 1]/Subscript[k, 1]+Subscript[\[Mu], 2]/Subscript[k, 2]) -> \[Nu]=-0.8, overall becomes Subscript[b, 1]*Subscript[f, NL]*\[ScriptF]/2 with again k^3 as dimension!*)
(*\[Mu]^2: PNEWabb1\[ScriptF] -> PNEWabb1\[ScriptF]SIMPLE*)


(* ::Input:: *)
(*PNEWabb1\[ScriptF]A[n1_,n2_]:=(((-3+n1+n2) Csc[n1 \[Pi]] Csc[n2 \[Pi]] Csc[(n1+n2) \[Pi]] ((64 n1^5 n2+n2 (-3+2 n2) (1+2 n2) (13+8 (-3+n2) n2)+8 n1^4 (-1+16 (-2+n2) n2)+4 n1^3 (5+78 n2-72 n2^2)-2 n1^2 (5+2 n2 (-5+2 n2)) (1+4 n2 (3+4 n2))+n1 (-11+2 n2 (5+2 n2 (33-2 n2 (47+4 n2 (-9+2 n2)))))) Sin[2 n1 \[Pi]]+(n1 (1+2 n1) (101+4 (-3+n1) n1 (29+8 (-3+n1) n1))+(9+2 n1 (-67+2 n1 (-23+2 n1 (81+4 n1 (-17+4 n1))))) n2+4 (-1+n1) (-3+2 n1 (-9+4 n1)) n2^2-4 (3+2 n1 (17+4 n1 (-9+4 n1))) n2^3-64 (-1+n1) n1 n2^4) Sin[2 n2 \[Pi]]+(n1 (1+2 n1) (11+4 (-3+n1) n1)+(-9+2 n1 (-25+2 n1 (-23+2 n1 (41+8 (-4+n1) n1)))) n2+4 (-1+n1) (3+2 n1 (-11+8 (-1+n1) n1)) n2^2-4 (-3+2 n1) (1+4 n1 (-3+2 n1)) n2^3-64 (-1+n1) n1 n2^4) Sin[2 (n1+n2) \[Pi]]))/(n1 (-3+2 n1) (-1+2 n1) (1+2 n1) n2 (-5+2 n1+2 n2) (-1+4 n2^2)))A1[n1,n2];*)


(* ::Input:: *)
(*PNEWabb1\[ScriptF]B[n1_,n2_]:=(-(1/((-1+2 n1) n2 (-1+4 n2^2)))(-2+n1+n2) Csc[n1 \[Pi]] Csc[n2 \[Pi]] Csc[(n1+n2) \[Pi]] ((-1+8 n2 (-1+n1+n2)) Sin[2 n1 \[Pi]]+(7-8 n2+8 n1 (-2+n1+n2)) Sin[2 n2 \[Pi]]+(1+8 (-1+n1) n2) Sin[2 (n1+n2) \[Pi]])) (Gamma[3/2-n1] Gamma[3/2-n2] Gamma[-(3/2)+n1+n2])/(8 \[Pi]^(3/2) Gamma[n1] Gamma[3-n1-n2] Gamma[n2]);*)


(* ::Input:: *)
(*PNEWabb1\[ScriptF][n1_,n2_]:=PNEWabb1\[ScriptF]A[n1,n2]+PNEWabb1\[ScriptF]B[n1,n2];*)


(* ::Program:: *)
(*PNEWabb1\[ScriptF]SIMPLE[n1_,n2_]=FullSimplify[FunctionExpand[PNEWabb1\[ScriptF][n1,n2]/\[ScriptCapitalI][n1,n2]]]*)


(* ::Program:: *)
(*(Csc[n1 \[Pi]] Csc[n2 \[Pi]] (-2 (64 n1^4 (-1+n2) n2+3 n2 (3-4 (-1+n2) n2)+4 n1^3 (-3+2 n2 (15+4 n2 (-7+4 n2)))+4 n1^2 (3+4 n2 (3+n2 (7+2 n2 (-7+2 n2))))+n1 (9-8 n2 (16+n2 (-6+n2 (-15+8 n2))))) Cos[(n1+n2) \[Pi]]-Csc[(n1+n2) \[Pi]] ((39 n2+64 n1^4 (-1+n2) n2-4 n2^2 (5+n2 (31+8 (-4+n2) n2))+4 n1^2 (-3+4 (-1+n2) n2 (-7+6 n2 (-3+2 n2)))+4 n1^3 (3+2 n2 (13+8 n2 (-5+3 n2)))+n1 (-9+8 n2 (-19+n2 (21+n2 (25+8 (-4+n2) n2))))) Sin[2 n1 \[Pi]]+(-n1 (-3+2 n1) (1+2 n1) (13+8 (-3+n1) n1)+(-9+8 n1 (-19+n1 (21+n1 (25+8 (-4+n1) n1)))) n2+4 (-3+4 (-1+n1) n1 (-7+6 n1 (-3+2 n1))) n2^2+4 (3+2 n1 (13+8 n1 (-5+3 n1))) n2^3+64 (-1+n1) n1 n2^4) Sin[2 n2 \[Pi]])))/(2 n1 (-1+4 n1^2) n2 (-5+2 n1+2 n2) (-1+4 n2^2))*)


(* ::Input:: *)
(*PNEWabb1\[ScriptF]SIMPLE[n1_,n2_]:=(Csc[n1 \[Pi]] Csc[n2 \[Pi]] (-2 (64 n1^4 (-1+n2) n2+3 n2 (3-4 (-1+n2) n2)+4 n1^3 (-3+2 n2 (15+4 n2 (-7+4 n2)))+4 n1^2 (3+4 n2 (3+n2 (7+2 n2 (-7+2 n2))))+n1 (9-8 n2 (16+n2 (-6+n2 (-15+8 n2))))) Cos[(n1+n2) \[Pi]]-Csc[(n1+n2) \[Pi]] ((39 n2+64 n1^4 (-1+n2) n2-4 n2^2 (5+n2 (31+8 (-4+n2) n2))+4 n1^2 (-3+4 (-1+n2) n2 (-7+6 n2 (-3+2 n2)))+4 n1^3 (3+2 n2 (13+8 n2 (-5+3 n2)))+n1 (-9+8 n2 (-19+n2 (21+n2 (25+8 (-4+n2) n2))))) Sin[2 n1 \[Pi]]+(-n1 (-3+2 n1) (1+2 n1) (13+8 (-3+n1) n1)+(-9+8 n1 (-19+n1 (21+n1 (25+8 (-4+n1) n1)))) n2+4 (-3+4 (-1+n1) n1 (-7+6 n1 (-3+2 n1))) n2^2+4 (3+2 n1 (13+8 n1 (-5+3 n1))) n2^3+64 (-1+n1) n1 n2^4) Sin[2 n2 \[Pi]])))/(2 n1 (-1+4 n1^2) n2 (-5+2 n1+2 n2) (-1+4 n2^2));*)


(* ::DisplayFormula:: *)
(*Subscript[Z, 2]\[Superset] Subscript[f, NL]*((\[ScriptF]*\[Mu]*k)/2)*\[ScriptF]*((Subscript[\[Mu], 1]*\!\(\*SubsuperscriptBox[\(\[Mu]\), \(2\), \(2\)]\))/Subscript[k, 1]+(Subscript[\[Mu], 2]*\!\(\*SubsuperscriptBox[\(\[Mu]\), \(1\), \(2\)]\))/Subscript[k, 2]) -> \[Nu]=-0.8, overall becomes Subscript[f, NL]*\[ScriptF]^2/2 with again k^3 as dimension!*)
(*\[Mu]^2: PNEWab\[Mu]1 -> PNEWab\[Mu]2SIMPLE*)
(*\[Mu]^4: PNEWab\[Mu]3 -> PNEWab\[Mu]4SIMPLE*)


(* ::Input:: *)
(*PNEWab\[ScriptF]2A1[n1_,n2_]:=(-(((-3+n1+n2) (-2+n1+n2) Csc[n1 \[Pi]] Csc[n2 \[Pi]] Csc[(n1+n2) \[Pi]] ((-3+4 n1-4 n1^2+2 (-1+2 n1) (1+8 (-1+n1) n1) n2+8 (3+8 (-1+n1) n1) n2^2+16 (-1+2 n1) n2^3) Sin[2 n1 \[Pi]]+(-3+6 n2+4 n1 (-1+n1+n2) (-3-8 n2+8 n1 (-1+n1+n2))) Sin[2 n2 \[Pi]]+(3-6 n2+4 n1 (-1+n1+(7+4 n1 (-3+2 n1)) n2+8 (-1+n1) n2^2)) Sin[2 (n1+n2) \[Pi]]))/(n1 (-3+2 n1) (-1+2 n1) (1+2 n1) n2 (-1+4 n2^2))))A1[n1,n2];*)


(* ::Input:: *)
(*PNEWab\[ScriptF]2A2[n1_,n2_]:=(-(((-4+n1+n2) (-3+n1+n2) Csc[n1 \[Pi]] Csc[n2 \[Pi]] Csc[(n1+n2) \[Pi]] ((15+128 n1^3 (-1+n2) n2-8 (-2+n2) n2 (1+4 (-2+n2) n2)+8 n1 (-3+2 n2) (2+n2 (11+8 (-3+n2) n2))+8 n1^2 (3+2 n2 (23+2 n2 (-21+8 n2)))) Sin[2 n1 \[Pi]]+(15-48 n2+8 (-(-2+n1) n1 (1+4 (-2+n1) n1)+n1 (-29+2 n1 (47+4 n1 (-9+2 n1))) n2+(3+2 n1 (23+2 n1 (-21+8 n1))) n2^2+16 (-1+n1) n1 n2^3)) Sin[2 n2 \[Pi]]+(-3 (5+8 (-2+n1) n1)-8 (-3+2 n1) (2+n1 (-13+8 n1)) n2+8 (-3+2 n1 (-5+2 n1) (-5+4 n1)) n2^2+128 (-1+n1) n1 n2^3) Sin[2 (n1+n2) \[Pi]]))/(n1 (-3+2 n1) (-1+2 n1) (1+2 n1) n2 (-3+2 n2) (-1+2 n2) (1+2 n2))))A2[n1,n2];*)


(* ::Input:: *)
(*PNEWab\[ScriptF]2B2[n1_,n2_]:=(-(((-4+n1+n2) (-3+n1+n2) Csc[n1 \[Pi]] Csc[n2 \[Pi]] Csc[(n1+n2) \[Pi]] ((-45-48 n2+4 (-(-2+n1) n1 (13+4 (-2+n1) n1)+n1 (103+2 n1 (-135+8 n1 (19+2 (-5+n1) n1))) n2+2 (51+n1 (-207+2 n1 (157+4 n1 (-25+6 n1)))) n2^2+8 (-12+n1 (43+4 n1 (-10+3 n1))) n2^3+8 (3+2 n1 (-5+2 n1)) n2^4)) Sin[2 n1 \[Pi]]+(4 (-2+n1) n1 (-29+4 (-2+n1) n1 (3+8 (-2+n1) n1))+4 n1 (31+2 n1 (-263+4 n1 (107+2 n1 (-31+6 n1)))) n2+8 (-9+n1 (-95+2 n1 (135+8 n1 (-13+3 n1)))) n2^2+64 (-1+n1) n1 (-5+2 n1) n2^3+9 (-5+16 n2)) Sin[2 n2 \[Pi]]+(45+4 (-2+n1) n1 (13+4 (-2+n1) n1)-144 n2+4 n1 (169+2 n1 (-161+16 n1 (10+(-5+n1) n1))) n2+8 (9+n1 (-113+2 n1 (101+4 n1 (-17+4 n1)))) n2^2+64 (-1+n1) n1 (-5+2 n1) n2^3) Sin[2 (n1+n2) \[Pi]]))/(n1 (-5+2 n1) (-3+2 n1) (-1+2 n1) (1+2 n1) n2 (-1+2 n2) (1+2 n2) (-5+2 n1+2 n2))))B2[n1,n2];*)


(* ::Input:: *)
(*PNEWab\[Mu]1[n1_,n2_]:=-PNEWab\[ScriptF]2A2[n1,n2];*)
(*PNEWab\[Mu]3[n1_,n2_]:=PNEWab\[ScriptF]2A1[n1,n2]-PNEWab\[ScriptF]2B2[n1,n2];*)


(* ::Program:: *)
(*PNEWab\[Mu]2SIMPLE[n1_,n2_]=FullSimplify[FunctionExpand[PNEWab\[Mu]1[n1,n2]/\[ScriptCapitalI][n1,n2]]]*)


(* ::Program:: *)
(*(Csc[n1 \[Pi]] Csc[n2 \[Pi]] ((-30+96 n2+16 (-3 (-2+n1) n1+n1 (-43+2 (25-8 n1) n1) n2+(-3+2 n1 (-5+2 n1) (-5+4 n1)) n2^2+16 (-1+n1) n1 n2^3)) Cos[(n1+n2) \[Pi]]+Csc[(n1+n2) \[Pi]] ((15+128 n1^3 (-1+n2) n2-8 (-2+n2) n2 (1+4 (-2+n2) n2)+8 n1 (-3+2 n2) (2+n2 (11+8 (-3+n2) n2))+8 n1^2 (3+2 n2 (23+2 n2 (-21+8 n2)))) Sin[2 n1 \[Pi]]+15 Sin[2 n2 \[Pi]]+8 (-(-2+n1) n1 (1+4 (-2+n1) n1)+(-3+2 n1) (2+n1 (11+8 (-3+n1) n1)) n2+(3+2 n1 (23+2 n1 (-21+8 n1))) n2^2+16 (-1+n1) n1 n2^3) Sin[2 n2 \[Pi]])))/(4 n1 (-1+4 n1^2) n2 (-5+2 n1+2 n2) (-1+4 n2^2))*)


(* ::Input:: *)
(*PNEWab\[Mu]2SIMPLE[n1_,n2_]:=(Csc[n1 \[Pi]] Csc[n2 \[Pi]] ((-30+96 n2+16 (-3 (-2+n1) n1+n1 (-43+2 (25-8 n1) n1) n2+(-3+2 n1 (-5+2 n1) (-5+4 n1)) n2^2+16 (-1+n1) n1 n2^3)) Cos[(n1+n2) \[Pi]]+Csc[(n1+n2) \[Pi]] ((15+128 n1^3 (-1+n2) n2-8 (-2+n2) n2 (1+4 (-2+n2) n2)+8 n1 (-3+2 n2) (2+n2 (11+8 (-3+n2) n2))+8 n1^2 (3+2 n2 (23+2 n2 (-21+8 n2)))) Sin[2 n1 \[Pi]]+15 Sin[2 n2 \[Pi]]+8 (-(-2+n1) n1 (1+4 (-2+n1) n1)+(-3+2 n1) (2+n1 (11+8 (-3+n1) n1)) n2+(3+2 n1 (23+2 n1 (-21+8 n1))) n2^2+16 (-1+n1) n1 n2^3) Sin[2 n2 \[Pi]])))/(4 n1 (-1+4 n1^2) n2 (-5+2 n1+2 n2) (-1+4 n2^2));*)


(* ::Program:: *)
(*PNEWab\[Mu]4SIMPLE[n1_,n2_]=FullSimplify[FunctionExpand[PNEWab\[Mu]3[n1,n2]/\[ScriptCapitalI][n1,n2]]]*)


(* ::Program:: *)
(*(Csc[n1 \[Pi]] Csc[n2 \[Pi]] (-2 (15-6 n1 (5+4 (-2+n1) n1)-30 n2+8 n1 (11-2 n1 (19+n1 (-23+8 n1))) n2+16 (-1+n1) (-3+4 n1 (4+n1 (-7+2 n1))) n2^2+8 (-3+2 n1 (23+4 n1 (-9+4 n1))) n2^3+128 (-1+n1) n1 n2^4) Cos[(n1+n2) \[Pi]]-Csc[(n1+n2) \[Pi]] ((-15+6 n1 (5+4 (-2+n1) n1)+62 n2-8 n1 (9+2 n1 (9+n1 (-21+8 n1))) n2+32 (3+n1 (-13+4 (-4+n1) (-2+n1) n1)) n2^2+8 (-47+2 n1 (61+4 n1 (-19+6 n1))) n2^3+32 (9+4 n1 (-5+3 n1)) n2^4+64 (-1+2 n1) n2^5) Sin[2 n1 \[Pi]]+(-15+30 n2+2 (n1 (31+4 n1 (12+n1 (-47+4 (9-2 n1) n1)))+4 n1 (-9+2 (-2+n1) n1 (13+8 (-3+n1) n1)) n2+8 (-1+n1) (3+4 n1 (3+n1 (-13+6 n1))) n2^2+12 (1+2 n1 (7+8 (-2+n1) n1)) n2^3+64 (-1+n1) n1 n2^4)) Sin[2 n2 \[Pi]])))/(4 n1 (-1+4 n1^2) n2 (-5+2 n1+2 n2) (-1+4 n2^2))*)


(* ::Input:: *)
(*PNEWab\[Mu]4SIMPLE[n1_,n2_]:=(Csc[n1 \[Pi]] Csc[n2 \[Pi]] (-2 (15-6 n1 (5+4 (-2+n1) n1)-30 n2+8 n1 (11-2 n1 (19+n1 (-23+8 n1))) n2+16 (-1+n1) (-3+4 n1 (4+n1 (-7+2 n1))) n2^2+8 (-3+2 n1 (23+4 n1 (-9+4 n1))) n2^3+128 (-1+n1) n1 n2^4) Cos[(n1+n2) \[Pi]]-Csc[(n1+n2) \[Pi]] ((-15+6 n1 (5+4 (-2+n1) n1)+62 n2-8 n1 (9+2 n1 (9+n1 (-21+8 n1))) n2+32 (3+n1 (-13+4 (-4+n1) (-2+n1) n1)) n2^2+8 (-47+2 n1 (61+4 n1 (-19+6 n1))) n2^3+32 (9+4 n1 (-5+3 n1)) n2^4+64 (-1+2 n1) n2^5) Sin[2 n1 \[Pi]]+(-15+30 n2+2 (n1 (31+4 n1 (12+n1 (-47+4 (9-2 n1) n1)))+4 n1 (-9+2 (-2+n1) n1 (13+8 (-3+n1) n1)) n2+8 (-1+n1) (3+4 n1 (3+n1 (-13+6 n1))) n2^2+12 (1+2 n1 (7+8 (-2+n1) n1)) n2^3+64 (-1+n1) n1 n2^4)) Sin[2 n2 \[Pi]])))/(4 n1 (-1+4 n1^2) n2 (-5+2 n1+2 n2) (-1+4 n2^2));*)


(* ::Input:: *)
(*(* if I look only at matter and multiply by Subscript[Z, 1] \[Rule] b1\[Mu]powers => everything is multiplied by Subscript[f, NL], k^3, and the NGparam, and \[ScriptCapitalM][k]/k^2... *)*)


(* ::Input:: *)
(*(* this is the same as before... Again I compute w.r.t. the first of these below... *)*)


(* ::DisplayFormula:: *)
(*Subscript[Z, 2]\[Superset] Subscript[b, 1]*Subscript[f, NL]*Subscript[F, 2][Subscript[k, 1],Subscript[k, 2]] -> k^3 as the dimension, \[Nu]=-0.8*)
(*\[Mu]^0: PNEWabb1F2*)


(* ::DisplayFormula:: *)
(*Subscript[Z, 2]\[Superset] Subscript[f, NL]*\[ScriptF]*Subscript[G, 2][Subscript[k, 1],Subscript[k, 2]] -> k^3 as the dimension, \[Nu]=-0.8*)
(*\[Mu]^2: PNEWabG2*)


(* ::DisplayFormula:: *)
(*Subscript[Z, 2]\[Superset] Subscript[b, 1]*Subscript[f, NL]*((\[ScriptF]*\[Mu]*k)/2)*(Subscript[\[Mu], 1]/Subscript[k, 1]+Subscript[\[Mu], 2]/Subscript[k, 2]) -> \[Nu]=-0.8, overall becomes Subscript[b, 1]*Subscript[f, NL]*\[ScriptF]/2 with again k^3 as dimension!*)
(*\[Mu]^2: PNEWabb1\[ScriptF] -> PNEWabb1\[ScriptF]SIMPLE*)


(* ::DisplayFormula:: *)
(*Subscript[Z, 2]\[Superset] Subscript[f, NL]*((\[ScriptF]*\[Mu]*k)/2)*\[ScriptF]*((Subscript[\[Mu], 1]*\!\(\*SubsuperscriptBox[\(\[Mu]\), \(2\), \(2\)]\))/Subscript[k, 1]+(Subscript[\[Mu], 2]*\!\(\*SubsuperscriptBox[\(\[Mu]\), \(1\), \(2\)]\))/Subscript[k, 2]) -> \[Nu]=-0.8, overall becomes Subscript[f, NL]*\[ScriptF]^2/2 with again k^3 as dimension!*)
(*\[Mu]^2: PNEWab\[Mu]1 -> PNEWab\[Mu]2SIMPLE*)
(*\[Mu]^4: PNEWab\[Mu]3 -> PNEWab\[Mu]4SIMPLE*)


(* ::Input:: *)
(*b1\[Mu]powers[0,0]=0;*)
(*b1\[Mu]powers[0,2]=Simplify[(\[ScriptF]*PNEWabG2[n1,n2]+\[ScriptF]^2/2*PNEWab\[Mu]2SIMPLE[n1,n2])/PNEWab[n1,n2]];*)
(*(*Print["Done!"];*)*)
(*b1\[Mu]powers[0,4]=Simplify[(\[ScriptF]^2/2*PNEWab\[Mu]4SIMPLE[n1,n2])/PNEWab[n1,n2]];*)
(*(*Print["Done!"];*)*)
(*b1\[Mu]powers[1,0]=Simplify[PNEWabb1F2[n1,n2]/PNEWab[n1,n2]];*)
(*(*Print["Done!"];*)*)
(*b1\[Mu]powers[1,2]=Simplify[(\[ScriptF]/2*PNEWabb1\[ScriptF]SIMPLE[n1,n2])/PNEWab[n1,n2]];*)
(*(*Print["Done!"];*)*)
(*b1\[Mu]powers[1,4]=0;*)


(* ::Input:: *)
(*b1\[Mu]powers[1,0]*)


(* ::Input:: *)
(*(* a check! *)*)


(* ::Input:: *)
(*totpreZ1=Sum[b1\[Mu]powers[i,j]*b1^i*\[Mu]^j,{i,0,1},{j,{0,2,4}}];*)
(*tot=(b1+\[ScriptF]*\[Mu]^2)*totpreZ1;*)


(* ::Input:: *)
(*rules=CoefficientRules[tot,\[Mu]]*)


(* ::Input:: *)
(*Do[CC[j]={j}/.rules;,{j,{0,2,4,6}}];*)


(* ::Input:: *)
(*\[Mu]coeffs={CC[0],CC[2],CC[4],CC[6]};*)


(* ::Input:: *)
(*F=c[0]+c[2]\[Mu]^2+c[4]\[Mu]^4+c[6]\[Mu]^6;*)


(* ::Input:: *)
(*basis[n_,\[Mu]_]=(2n+1)/2 LegendreP[n,\[Mu]];*)


(* ::Input:: *)
(*Integrate[F basis[0,\[Mu]],{\[Mu],-1,1}];*)
(*Integrate[F basis[2,\[Mu]],{\[Mu],-1,1}];*)
(*Integrate[F basis[4,\[Mu]],{\[Mu],-1,1}];*)
(*Integrate[F basis[6,\[Mu]],{\[Mu],-1,1}];*)


(* ::Input:: *)
(*coeffs[0]=D[Integrate[F basis[0,\[Mu]],{\[Mu],-1,1}],{{c[0],c[2],c[4],c[6]}}]*)
(*coeffs[2]=D[Integrate[F basis[2,\[Mu]],{\[Mu],-1,1}],{{c[0],c[2],c[4],c[6]}}]*)
(*coeffs[4]=D[Integrate[F basis[4,\[Mu]],{\[Mu],-1,1}],{{c[0],c[2],c[4],c[6]}}]*)


(* ::Input:: *)
(*(* here I can do like this. For AP, I need numerical integration since I have all powers of \[Mu]... *)*)


(* ::Input:: *)
(*\[Mu]0=Sum[coeffs[0][[i]]\[Mu]coeffs[[i]],{i,1,4}];*)
(*\[Mu]2=Sum[coeffs[2][[i]]\[Mu]coeffs[[i]],{i,1,4}];*)
(*\[Mu]4=Sum[coeffs[4][[i]]\[Mu]coeffs[[i]],{i,1,4}];*)


(* ::Input:: *)
(*rulesb1\[Mu]0=CoefficientRules[\[Mu]0,b1];*)
(*rulesb1\[Mu]2=CoefficientRules[\[Mu]2,b1];*)
(*rulesb1\[Mu]4=CoefficientRules[\[Mu]4,b1];*)


(* ::Input:: *)
(*(* monopole *)*)


(* ::Input:: *)
(*(* b1^0 *)*)


(* ::Input:: *)
(*vv0=Simplify[{0}/.rulesb1\[Mu]0];*)


(* ::Input:: *)
(*(* b1^1 *)*)


(* ::Input:: *)
(*vd0=Simplify[{1}/.rulesb1\[Mu]0];*)


(* ::Input:: *)
(*(* b1^2 *)*)


(* ::Input:: *)
(*dd0=Simplify[{2}/.rulesb1\[Mu]0];*)


(* ::Input:: *)
(*(* quadrupole *)*)


(* ::Input:: *)
(*(* b1^0 *)*)


(* ::Input:: *)
(*vv2=Simplify[{0}/.rulesb1\[Mu]2];*)


(* ::Input:: *)
(*(* b1^1 *)*)


(* ::Input:: *)
(*vd2=Simplify[{1}/.rulesb1\[Mu]2];*)


(* ::Input:: *)
(*(* b1^2 *)*)


(* ::Input:: *)
(*dd2=Simplify[{2}/.rulesb1\[Mu]2];*)


(* ::Input:: *)
(*(* hexadecapole *)*)


(* ::Input:: *)
(*(* b1^0 *)*)


(* ::Input:: *)
(*vv4=Simplify[{0}/.rulesb1\[Mu]4];*)


(* ::Input:: *)
(*(* b1^1 *)*)


(* ::Input:: *)
(*vd4=Simplify[{1}/.rulesb1\[Mu]4];*)


(* ::Input:: *)
(*(* b1^2 *)*)


(* ::Input:: *)
(*Simplify[{2}/.rulesb1\[Mu]4]*)


(* ::Input:: *)
(*dd4=0;*)


(* ::Section:: *)
(*Summary*)


(* ::Text:: *)
(*We have vv0, vd0, dd0; vv2, vd2, dd2; vv4, vd4 BUT NOT dd4... They are called via these names, and \[ScriptF] is the growth rate, and are functions of n1, n2 as names...*)


(* ::Section:: *)
(*Matrices -- 2*)


(* ::Input:: *)
(*total={};*)


(* ::Subsection:: *)
(*vv0*)


(* ::Input:: *)
(*\[ScriptF]rules=CoefficientRules[vv0,\[ScriptF]]//Simplify;*)


(* ::Input:: *)
(*whichpowers=\[ScriptF]rules[[All,1]]//Flatten//Sort*)


(* ::Input:: *)
(*Do[vv0\[ScriptF][iii][n1_,n2_]=Simplify[({iii}/.\[ScriptF]rules)*PNEWab[n1,n2]];,{iii,whichpowers}];*)


(* ::Input:: *)
(*Do[AppendTo[total,{"vv0, f = "<>ToString[iii],vv0\[ScriptF][iii][n1,n2]}],{iii,whichpowers}];*)


(* ::Program:: *)
(*Do[Export[NotebookDirectory[] <> "vv0_f" <> ToString[iii](*<>"-"<>ToString[Nmax]*)<> ".txt", vv0\[ScriptF][iii][n1, n2]];, {iii, whichpowers}]*)


(* ::Input:: *)
(*total*)


(* ::Subsection:: *)
(*vv2*)


(* ::Input:: *)
(*\[ScriptF]rules=CoefficientRules[vv2,\[ScriptF]]//Simplify;*)


(* ::Input:: *)
(*whichpowers=\[ScriptF]rules[[All,1]]//Flatten//Sort*)


(* ::Input:: *)
(*Do[vv2\[ScriptF][iii][n1_,n2_]=Simplify[({iii}/.\[ScriptF]rules)*PNEWab[n1,n2]];,{iii,whichpowers}];*)


(* ::Input:: *)
(*Do[AppendTo[total,{"vv2, f = "<>ToString[iii],vv2\[ScriptF][iii][n1,n2]}],{iii,whichpowers}];*)


(* ::Program:: *)
(*Do[Export[NotebookDirectory[] <> "vv2_f" <> ToString[iii](*<>"-"<>ToString[Nmax]*)<> ".txt", vv2\[ScriptF][iii][n1, n2]];, {iii, whichpowers}]*)


(* ::Input:: *)
(*total*)


(* ::Subsection:: *)
(*vv4*)


(* ::Input:: *)
(*\[ScriptF]rules=CoefficientRules[vv4,\[ScriptF]]//Simplify;*)


(* ::Input:: *)
(*whichpowers=\[ScriptF]rules[[All,1]]//Flatten//Sort*)


(* ::Input:: *)
(*Do[vv4\[ScriptF][iii][n1_,n2_]=Simplify[({iii}/.\[ScriptF]rules)*PNEWab[n1,n2]];,{iii,whichpowers}];*)


(* ::Input:: *)
(*Do[AppendTo[total,{"vv4, f = "<>ToString[iii],vv4\[ScriptF][iii][n1,n2]}],{iii,whichpowers}];*)


(* ::Program:: *)
(*Do[Export[NotebookDirectory[] <> "vv4_f" <> ToString[iii](*<>"-"<>ToString[Nmax]*)<> ".txt", vv4\[ScriptF][iii][n1, n2]];, {iii, whichpowers}]*)


(* ::Input:: *)
(*total*)


(* ::Subsection:: *)
(*vd0*)


(* ::Input:: *)
(*\[ScriptF]rules=CoefficientRules[vd0,\[ScriptF]]//Simplify;*)


(* ::Input:: *)
(*whichpowers=\[ScriptF]rules[[All,1]]//Flatten//Sort*)


(* ::Input:: *)
(*Do[vd0\[ScriptF][iii][n1_,n2_]=Simplify[({iii}/.\[ScriptF]rules)*PNEWab[n1,n2]];,{iii,whichpowers}];*)


(* ::Input:: *)
(*Do[AppendTo[total,{"vd0, f = "<>ToString[iii],vd0\[ScriptF][iii][n1,n2]}],{iii,whichpowers}];*)


(* ::Program:: *)
(*Do[Export[NotebookDirectory[] <> "vd0_f" <> ToString[iii](*<>"-"<>ToString[Nmax]*)<> ".txt", vd0\[ScriptF][iii][n1, n2]];, {iii, whichpowers}]*)


(* ::Input:: *)
(*total*)


(* ::Subsection:: *)
(*vd2*)


(* ::Input:: *)
(*\[ScriptF]rules=CoefficientRules[vd2,\[ScriptF]]//Simplify;*)


(* ::Input:: *)
(*whichpowers=\[ScriptF]rules[[All,1]]//Flatten//Sort*)


(* ::Input:: *)
(*Do[vd2\[ScriptF][iii][n1_,n2_]=Simplify[({iii}/.\[ScriptF]rules)*PNEWab[n1,n2]];,{iii,whichpowers}];*)


(* ::Input:: *)
(*Do[AppendTo[total,{"vd2, f = "<>ToString[iii],vd2\[ScriptF][iii][n1,n2]}],{iii,whichpowers}];*)


(* ::Program:: *)
(*Do[Export[NotebookDirectory[] <> "vd2_f" <> ToString[iii](*<>"-"<>ToString[Nmax]*)<> ".txt", vd2\[ScriptF][iii][n1, n2]];, {iii, whichpowers}]*)


(* ::Input:: *)
(*total*)


(* ::Subsection:: *)
(*vd4*)


(* ::Input:: *)
(*\[ScriptF]rules=CoefficientRules[vd4,\[ScriptF]]//Simplify;*)


(* ::Input:: *)
(*whichpowers=\[ScriptF]rules[[All,1]]//Flatten//Sort*)


(* ::Input:: *)
(*Do[vd4\[ScriptF][iii][n1_,n2_]=Simplify[({iii}/.\[ScriptF]rules)*PNEWab[n1,n2]];,{iii,whichpowers}];*)


(* ::Input:: *)
(*Do[AppendTo[total,{"vd4, f = "<>ToString[iii],vd4\[ScriptF][iii][n1,n2]}],{iii,whichpowers}];*)


(* ::Program:: *)
(*Do[Export[NotebookDirectory[] <> "vd4_f" <> ToString[iii](*<>"-"<>ToString[Nmax]*)<> ".txt", vd4\[ScriptF][iii][n1, n2]];, {iii, whichpowers}]*)


(* ::Input:: *)
(*total*)


(* ::Subsection:: *)
(*dd0*)


(* ::Input:: *)
(*\[ScriptF]rules=CoefficientRules[dd0,\[ScriptF]]//Simplify;*)


(* ::Input:: *)
(*whichpowers=\[ScriptF]rules[[All,1]]//Flatten//Sort*)


(* ::Input:: *)
(*Do[dd0\[ScriptF][iii][n1_,n2_]=Simplify[({iii}/.\[ScriptF]rules)*PNEWab[n1,n2]];,{iii,whichpowers}];*)


(* ::Input:: *)
(*Do[AppendTo[total,{"dd0, f = "<>ToString[iii],dd0\[ScriptF][iii][n1,n2]}],{iii,whichpowers}];*)


(* ::Program:: *)
(*Do[Export[NotebookDirectory[] <> "dd0_f" <> ToString[iii](*<>"-"<>ToString[Nmax]*)<> ".txt", dd0\[ScriptF][iii][n1, n2]];, {iii, whichpowers}]*)


(* ::Input:: *)
(*total*)


(* ::Subsection:: *)
(*dd2*)


(* ::Input:: *)
(*\[ScriptF]rules=CoefficientRules[dd2,\[ScriptF]]//Simplify;*)


(* ::Input:: *)
(*whichpowers=\[ScriptF]rules[[All,1]]//Flatten//Sort*)


(* ::Input:: *)
(*Do[dd2\[ScriptF][iii][n1_,n2_]=Simplify[({iii}/.\[ScriptF]rules)*PNEWab[n1,n2]];,{iii,whichpowers}];*)


(* ::Input:: *)
(*Do[AppendTo[total,{"dd2, f = "<>ToString[iii],dd2\[ScriptF][iii][n1,n2]}],{iii,whichpowers}];*)


(* ::Program:: *)
(*Do[Export[NotebookDirectory[] <> "dd2_f" <> ToString[iii](*<>"-"<>ToString[Nmax]*)<> ".txt", dd2\[ScriptF][iii][n1, n2]];, {iii, whichpowers}]*)


(* ::Input:: *)
(*total*)


(* ::Subsection:: *)
(*dd4*)


(* ::Input:: *)
(*\[ScriptF]rules=CoefficientRules[dd4,\[ScriptF]]//Simplify;*)


(* ::Input:: *)
(*whichpowers=\[ScriptF]rules[[All,1]]//Flatten//Sort*)


(* ::Input:: *)
(*(* they are 14... Then technically I need the powers of \[Mu] as well, for how the code is written... There are less there... I move to MontePython for fun now... *)*)


(* ::Input:: *)
(*total*)


(* ::Input:: *)
(*total//Length*)


(* ::Section:: *)
(*Matrices -- 3*)


(* ::Input:: *)
(*length=Length[total]*)


(* ::Input:: *)
(*Monitor[Do[If[ii<jj&&RationalExpressionQ[Simplify[total[[ii]][[2]]/total[[jj]][[2]]],{n1,n2}],Print[{ii,jj}]],{ii,1,length},{jj,1,length}],{ii,jj}]*)


(* ::Input:: *)
(*(* actually it pays to do it... We see that we just need to export -> *)*)


(* ::Input:: *)
(*total[[3]][[2]]/total[[1]][[2]]*)


(* ::Input:: *)
(*total[[3]][[1]]*)
(*total[[1]][[1]]*)


(* ::Input:: *)
(*(*/////////////////////////////*)*)


(* ::Input:: *)
(*total[[5]][[2]]/total[[1]][[2]]*)


(* ::Input:: *)
(*total[[5]][[1]]*)
(*total[[1]][[1]]*)


(* ::Input:: *)
(*(*/////////////////////////////*)*)


(* ::Input:: *)
(*total[[9]][[2]]/total[[7]][[2]]*)


(* ::Input:: *)
(*total[[9]][[1]]*)
(*total[[7]][[1]]*)


(* ::Input:: *)
(*(*/////////////////////////////*)*)


(* ::Input:: *)
(*total[[14]][[2]]/total[[13]][[2]]*)


(* ::Input:: *)
(*total[[14]][[1]]*)
(*total[[13]][[1]]*)


(* ::Input:: *)
(*(* again, I will recheck these once I do the final code checks... *)*)


(* ::CodeText::Regular:: *)
(*(openblas_test) gcabass@mbp156 ~/Documents/work/CLASS-PT-master-fNL/pt_matrices/compute_matrices_python/matter_multipoles % more simplifications.txt*)
(*-) vv2, f = 2 is [ 20/7 ] of vv0, f = 2*)
(*-) vv4, f = 2 is [ 8/7  ] of vv0, f = 2*)
(*-) vd2, f = 1 is [ 2    ] of vd0, f = 1*)
(*-) dd2, f = 1 is [ 2    ] of dd0, f = 1*)
(*(openblas_test) gcabass@mbp156 ~/Documents/work/CLASS-PT-master-fNL/pt_matrices/compute_matrices_python/matter_multipoles % more whichdict.py *)
(*#--> vv0_f2.txt, vv0_f3.txt, vv2_f2.txt, vv2_f3.txt, vv4_f2.txt, vv4_f3.txt, vd0_f1.txt, vd0_f2.txt, vd2_f1.txt, vd2_f2.txt, vd4_f2.txt, dd0_f0.txt, dd0_f1.txt, dd2_f1.txt*)
(**)
(*#-) vv2, f = 2 is [ 20/7 ] of vv0, f = 2*)
(*#-) vv4, f = 2 is [ 8/7  ] of vv0, f = 2*)
(*#-) vd2, f = 1 is [ 2    ] of vd0, f = 1*)
(*#-) dd2, f = 1 is [ 2    ] of dd0, f = 1*)
(**)
(*importdict={"1":"vv0_f2.txt","2":"vv0_f3.txt","3":"vv2_f3.txt","4":"vv4_f3.txt","5":"vd0_f1.txt","6":"vd0_f2.txt","7":"vd2_f2.txt","8":"vd4_f2.txt","9":"dd0_f0.txt","10":"dd0_f1.txt"}*)


(* ::Input:: *)
(*(* recall I had all apart dd4... And -> *)*)


(* ::Program:: *)
(*vv0*)
(*whichpowers=\[ScriptF]rules[[All,1]]//Flatten//Sort*)
(*{2,3}*)
(*vv2*)
(*whichpowers=\[ScriptF]rules[[All,1]]//Flatten//Sort*)
(*{2,3}*)
(*vv4*)
(*whichpowers=\[ScriptF]rules[[All,1]]//Flatten//Sort*)
(*{2,3}*)
(*vd0*)
(*whichpowers=\[ScriptF]rules[[All,1]]//Flatten//Sort*)
(*{1,2}*)
(*vd2*)
(*whichpowers=\[ScriptF]rules[[All,1]]//Flatten//Sort*)
(*{1,2}*)
(*vd4*)
(*whichpowers=\[ScriptF]rules[[All,1]]//Flatten//Sort*)
(*{2}*)
(*dd0*)
(*whichpowers=\[ScriptF]rules[[All,1]]//Flatten//Sort*)
(*{0,1}*)
(*dd2*)
(*whichpowers=\[ScriptF]rules[[All,1]]//Flatten//Sort*)
(*{1}*)
