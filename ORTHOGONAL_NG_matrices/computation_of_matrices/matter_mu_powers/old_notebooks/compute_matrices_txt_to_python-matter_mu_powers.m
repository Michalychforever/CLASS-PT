(* ::Package:: *)

(* ::Code:: *)
(*ClearAll[Evaluate[Context[]<>"*"]];*)
(*ClearSystemCache[];*)


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
(*(* THIS IS FOR 2Subscript[F, 2]^2: CAREFUL!!! This is the dd contribution at \[Mu]=0 for matter, and indeed b1\[Mu]powers[1,0] is 1 now. This will be the baseline... Should it be b1\[Mu]powers[2,0]? Whatever... Below rules... *)*)


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
(*b1\[Mu]powers[0,0]=0;*)
(*b1\[Mu]powers[0,2]=(f*PNEWabG2[nu1,nu2]+f^2/2*PNEWab\[Mu]2SIMPLE[nu1,nu2])/PNEWab[nu1,nu2];*)
(*b1\[Mu]powers[0,4]=(f^2/2*PNEWab\[Mu]4SIMPLE[nu1,nu2])/PNEWab[nu1,nu2];*)
(*b1\[Mu]powers[1,0]=PNEWabb1F2[nu1,nu2]/PNEWab[nu1,nu2];*)
(*b1\[Mu]powers[1,2]=(f/2*PNEWabb1\[ScriptF]SIMPLE[nu1,nu2])/PNEWab[nu1,nu2];*)
(*b1\[Mu]powers[1,4]=0;*)


(* ::Input:: *)
(*totpreZ1=Sum[b1\[Mu]powers[i,j]*b1^i*\[Mu]^j,{i,0,1},{j,{0,2,4}}];*)
(*tot=(b1+f*\[Mu]^2)*totpreZ1;*)


(* ::Input:: *)
(*rules=CoefficientRules[tot,\[Mu]]*)


(* ::Input:: *)
(*Do[CC[j]={j}/.rules;,{j,{0,2,4,6}}];*)


(* ::Input:: *)
(*\[Mu]coeffs={CC[0],CC[2],CC[4],CC[6]};*)


(* ::Input:: *)
(*rulesb1\[Mu]0=CoefficientRules[CC[0],b1];*)
(*rulesb1\[Mu]2=CoefficientRules[CC[2],b1];*)
(*rulesb1\[Mu]4=CoefficientRules[CC[4],b1];*)
(*rulesb1\[Mu]6=CoefficientRules[CC[6],b1];*)


(* ::Input:: *)
(*(* \[Mu]^0 *)*)


(* ::Input:: *)
(*(* b1^0 *)*)


(* ::Input:: *)
(*vv0=Simplify[{0}/.rulesb1\[Mu]0]*)
(*vv0=0;*)


(* ::Input:: *)
(*(* b1^1 *)*)


(* ::Input:: *)
(*vd0=Simplify[{1}/.rulesb1\[Mu]0]*)
(*vd0=0;*)


(* ::Input:: *)
(*(* b1^2 *)*)


(* ::Input:: *)
(*dd0=Simplify[{2}/.rulesb1\[Mu]0]*)


(* ::Input:: *)
(*(* so only need to start from \[Mu] squared... *)*)


(* ::Input:: *)
(*(* \[Mu]^2 *)*)


(* ::Input:: *)
(*(* b1^0 *)*)


(* ::Input:: *)
(*vv2=Simplify[{0}/.rulesb1\[Mu]2]*)
(*vv2=0;*)


(* ::Input:: *)
(*(* b1^1 *)*)


(* ::Input:: *)
(*vd2=Simplify[{1}/.rulesb1\[Mu]2]*)


(* ::Input:: *)
(*(* b1^2 *)*)


(* ::Input:: *)
(*dd2=Simplify[{2}/.rulesb1\[Mu]2]*)


(* ::Input:: *)
(*(* \[Mu]^4 *)*)


(* ::Input:: *)
(*(* b1^0 *)*)


(* ::Input:: *)
(*vv4=Simplify[{0}/.rulesb1\[Mu]4]*)


(* ::Input:: *)
(*(* b1^1 *)*)


(* ::Input:: *)
(*vd4=Simplify[{1}/.rulesb1\[Mu]4]*)


(* ::Input:: *)
(*(* b1^2 *)*)


(* ::Input:: *)
(*dd4=Simplify[{2}/.rulesb1\[Mu]4]*)
(*dd4=0;*)


(* ::Input:: *)
(*(* \[Mu]^6 *)*)


(* ::Input:: *)
(*(* b1^0 *)*)


(* ::Input:: *)
(*vv6=Simplify[{0}/.rulesb1\[Mu]6]*)


(* ::Input:: *)
(*(* b1^1 *)*)


(* ::Input:: *)
(*vd6=Simplify[{1}/.rulesb1\[Mu]6]*)
(*vd6=0;*)


(* ::Input:: *)
(*(* b1^2 *)*)


(* ::Input:: *)
(*dd6=Simplify[{2}/.rulesb1\[Mu]6]*)
(*dd6=0;*)


(* ::Input:: *)
(*total={};*)


(* ::Section:: *)
(*\[Mu]^0*)


(* ::Subsection:: *)
(*b1^0*)


(* ::Input:: *)
(*vv0*)


(* ::Subsection:: *)
(*b1^1*)


(* ::Input:: *)
(*vd0*)


(* ::Subsection:: *)
(*b1^2*)


(* ::Input:: *)
(*dd0*)


(* ::Input:: *)
(*(* so nothing needed... The .c code already accounts for this... *)*)


(* ::Section:: *)
(*\[Mu]^2*)


(* ::Subsection:: *)
(*b1^0*)


(* ::Input:: *)
(*vv2*)


(* ::Subsection:: *)
(*b1^1*)


(* ::Input:: *)
(*vd2*)


(* ::Input:: *)
(*\[ScriptF]rules=CoefficientRules[vd2,f]//Simplify;*)


(* ::Input:: *)
(*whichpowers=\[ScriptF]rules[[All,1]]//Flatten//Sort*)


(* ::Input:: *)
(*Do[vd2\[ScriptF][iii][nu1_,nu2_]=Simplify[({iii}/.\[ScriptF]rules)*PNEWab[nu1,nu2]];,{iii,whichpowers}];*)


(* ::Input:: *)
(*Do[AppendTo[total,{"vd2, f = "<>ToString[iii],vd2\[ScriptF][iii][nu1,nu2]}],{iii,whichpowers}];*)


(* ::Input:: *)
(*Do[Export[NotebookDirectory[]<>"vd2_f"<>ToString[iii]<>".txt",vd2\[ScriptF][iii][nu1,nu2]];,{iii,whichpowers}]*)


(* ::Subsection:: *)
(*b1^2*)


(* ::Input:: *)
(*dd2*)


(* ::Input:: *)
(*\[ScriptF]rules=CoefficientRules[dd2,f]//Simplify;*)


(* ::Input:: *)
(*whichpowers=\[ScriptF]rules[[All,1]]//Flatten//Sort*)


(* ::Input:: *)
(*Do[dd2\[ScriptF][iii][nu1_,nu2_]=Simplify[({iii}/.\[ScriptF]rules)*PNEWab[nu1,nu2]];,{iii,whichpowers}];*)


(* ::Input:: *)
(*Do[AppendTo[total,{"dd2, f = "<>ToString[iii],dd2\[ScriptF][iii][nu1,nu2]}],{iii,whichpowers}];*)


(* ::Input:: *)
(*Do[Export[NotebookDirectory[]<>"dd2_f"<>ToString[iii]<>".txt",dd2\[ScriptF][iii][nu1,nu2]];,{iii,whichpowers}]*)


(* ::Section:: *)
(*\[Mu]^4*)


(* ::Subsection:: *)
(*b1^0*)


(* ::Input:: *)
(*vv4*)


(* ::Input:: *)
(*\[ScriptF]rules=CoefficientRules[vv4,f]//Simplify;*)


(* ::Input:: *)
(*whichpowers=\[ScriptF]rules[[All,1]]//Flatten//Sort*)


(* ::Input:: *)
(*Do[vv4\[ScriptF][iii][nu1_,nu2_]=Simplify[({iii}/.\[ScriptF]rules)*PNEWab[nu1,nu2]];,{iii,whichpowers}];*)


(* ::Input:: *)
(*Do[AppendTo[total,{"vv4, f = "<>ToString[iii],vv4\[ScriptF][iii][nu1,nu2]}],{iii,whichpowers}];*)


(* ::Input:: *)
(*Do[Export[NotebookDirectory[]<>"vv4_f"<>ToString[iii]<>".txt",vv4\[ScriptF][iii][nu1,nu2]];,{iii,whichpowers}]*)


(* ::Subsection:: *)
(*b1^1*)


(* ::Input:: *)
(*vd4*)


(* ::Input:: *)
(*\[ScriptF]rules=CoefficientRules[vd4,f]//Simplify;*)


(* ::Input:: *)
(*whichpowers=\[ScriptF]rules[[All,1]]//Flatten//Sort*)


(* ::Input:: *)
(*Do[vd4\[ScriptF][iii][nu1_,nu2_]=Simplify[({iii}/.\[ScriptF]rules)*PNEWab[nu1,nu2]];,{iii,whichpowers}];*)


(* ::Input:: *)
(*Do[AppendTo[total,{"vd4, f = "<>ToString[iii],vd4\[ScriptF][iii][nu1,nu2]}],{iii,whichpowers}];*)


(* ::Input:: *)
(*Do[Export[NotebookDirectory[]<>"vd4_f"<>ToString[iii]<>".txt",vd4\[ScriptF][iii][nu1,nu2]];,{iii,whichpowers}]*)


(* ::Subsection:: *)
(*b1^2*)


(* ::Input:: *)
(*dd4*)


(* ::Section:: *)
(*\[Mu]^6*)


(* ::Subsection:: *)
(*b1^0*)


(* ::Input:: *)
(*vv6*)


(* ::Input:: *)
(*\[ScriptF]rules=CoefficientRules[vv6,f]//Simplify;*)


(* ::Input:: *)
(*whichpowers=\[ScriptF]rules[[All,1]]//Flatten//Sort*)


(* ::Input:: *)
(*Do[vv6\[ScriptF][iii][nu1_,nu2_]=Simplify[({iii}/.\[ScriptF]rules)*PNEWab[nu1,nu2]];,{iii,whichpowers}];*)


(* ::Input:: *)
(*Do[AppendTo[total,{"vv6, f = "<>ToString[iii],vv6\[ScriptF][iii][nu1,nu2]}],{iii,whichpowers}];*)


(* ::Input:: *)
(*Do[Export[NotebookDirectory[]<>"vv6_f"<>ToString[iii]<>".txt",vv6\[ScriptF][iii][nu1,nu2]];,{iii,whichpowers}]*)


(* ::Subsection:: *)
(*b1^1*)


(* ::Input:: *)
(*vd6*)


(* ::Subsection:: *)
(*b1^2*)


(* ::Input:: *)
(*dd6*)


(* ::Section:: *)
(*Check*)


(* ::Input:: *)
(*length=Length[total]*)


(* ::Input:: *)
(*Monitor[Do[If[ii<jj&&RationalExpressionQ[Simplify[total[[ii]][[2]]/total[[jj]][[2]]],{nu1,nu2}],Print[{ii,jj}]],{ii,1,length},{jj,1,length}],{ii,jj}]*)


(* ::Input:: *)
(*(* actually it pays to do it... We see that we just need to export -> *)*)


(* ::Input:: *)
(*total[[5]][[2]]/total[[2]][[2]]*)


(* ::Input:: *)
(*total[[5]][[1]]*)
(*total[[2]][[1]]*)


(* ::Input:: *)
(*(*/////////////////////////////*)*)


(* ::Input:: *)
(*(* notice that I had screwed up because I was using the matrix M22basic_oneline_complex of Misha, that uses a different bias! The factors of 2 from the permutation were included in the matrix "reduction" formula, but for me I need two different biases for Subscript[\[ScriptCapitalG], 2] and Subscript[b, 2]... So I will need two matrices in real space... For the case of multipoles I will check what to do. If they were written in terms of \[Mu] powers it would have been better. I doubt that a great reduction in numbers will occur... We will see... *)*)
