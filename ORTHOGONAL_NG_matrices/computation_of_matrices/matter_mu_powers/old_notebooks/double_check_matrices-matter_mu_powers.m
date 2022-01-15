(* ::Package:: *)

(* ::Code:: *)
(*ClearAll[Evaluate[Context[]<>"*"]];*)
(*ClearSystemCache[];*)


(* ::Title:: *)
(*OLD*)


(* ::Input:: *)
(*\[ScriptCapitalI][n1_,n2_]:=(Gamma[3/2-n1] Gamma[3/2-n2] Gamma[-(3/2)+n1+n2])/(8 \[Pi]^(3/2) Gamma[n1] Gamma[3-n1-n2] Gamma[n2]);*)


(* ::Input:: *)
(*J2[k_,\[Nu]1_,\[Nu]2_]:=1/(2\[Pi])^3 (Gamma[3/2-\[Nu]1]Gamma[3/2-\[Nu]2]Gamma[\[Nu]1+\[Nu]2-3/2])/(Gamma[\[Nu]1]Gamma[\[Nu]2]Gamma[3-\[Nu]1-\[Nu]2]) \[Pi]^(3/2) k^(3-2\[Nu]1-2\[Nu]2);*)
(*A1[\[Nu]1_,\[Nu]2_]:=1/2 (J2[1,\[Nu]1-1,\[Nu]2]-J2[1,\[Nu]1,\[Nu]2-1]+J2[1,\[Nu]1,\[Nu]2]);*)
(*A2[\[Nu]1_,\[Nu]2_]:=1/8 (-J2[1,-2+\[Nu]1,\[Nu]2]+2 J2[1,-1+\[Nu]1,-1+\[Nu]2]+2 J2[1,-1+\[Nu]1,\[Nu]2]-J2[1,\[Nu]1,-2+\[Nu]2]+2 J2[1,\[Nu]1,-1+\[Nu]2]-J2[1,\[Nu]1,\[Nu]2]);*)
(*B2[\[Nu]1_,\[Nu]2_]:=1/8 (3 J2[1,-2+\[Nu]1,\[Nu]2]-6 J2[1,-1+\[Nu]1,-1+\[Nu]2]+2 J2[1,-1+\[Nu]1,\[Nu]2]+3 J2[1,\[Nu]1,-2+\[Nu]2]-6 J2[1,\[Nu]1,-1+\[Nu]2]+3 J2[1,\[Nu]1,\[Nu]2]);*)


(* ::DisplayFormula:: *)
(*Subscript[Z, 2]\[Superset] Subscript[b, 1]*Subscript[f, NL]*Subscript[F, 2][Subscript[k, 1],Subscript[k, 2]] -> k^3 as the dimension, \[Nu]=-0.8*)
(*\[Mu]^0: PNEWabb1F2 -> ALREADY CHECKED...*)


(* ::Input:: *)
(*PNEWabb1F2[n1_,n2_]:=((Csc[n1 \[Pi]] Csc[n2 \[Pi]] (-2 (60+448 n1^4 (-1+n2) n2-3 n2 (43+4 n2 (-15+7 n2))+4 n1^3 (-21+2 n2 (169+4 n2 (-65+28 n2)))+4 n1^2 (45+4 n2 (-79+n2 (169+2 n2 (-65+14 n2))))+n1 (-129-8 n2 (-60+n2 (158+n2 (-169+56 n2))))) Cos[(n1+n2) \[Pi]]-Csc[(n1+n2) \[Pi]] ((-60+3 n1 (43+4 n1 (-15+7 n1))+209 n2-8 n1 (17+(-2+n1) n1 (-43+56 n1)) n2+4 (101+2 n1 (-229+2 n1 (245+4 n1 (-43+7 n1)))) n2^2+4 (-345+2 n1 (463+4 n1 (-137+42 n1))) n2^3+64 (16+3 n1 (-12+7 n1)) n2^4+224 (-1+2 n1) n2^5) Sin[2 n1 \[Pi]]+(-60+n1 (209+4 n1 (101+n1 (-345+8 (32-7 n1) n1)))+129 n2+8 n1 (-17+n1 (-229+n1 (463+8 n1 (-36+7 n1)))) n2+4 (-45+4 n1 (-43+n1 (245-274 n1+84 n1^2))) n2^2+4 (21+2 n1 (155+8 n1 (-43+21 n1))) n2^3+448 (-1+n1) n1 n2^4) Sin[2 n2 \[Pi]])))/(28 n1 (-1+4 n1^2) n2 (-5+2 n1+2 n2) (-1+4 n2^2)))(*\[ScriptCapitalI][n1,n2]*);*)


(* ::Input:: *)
(*PNEWabb1F2CHECK[n1_,n2_]:=(Csc[n1 \[Pi]] Csc[n2 \[Pi]] (-2 (60+448 n1^4 (-1+n2) n2-3 n2 (43+4 n2 (-15+7 n2))+4 n1^3 (-21+2 n2 (169+4 n2 (-65+28 n2)))+4 n1^2 (45+4 n2 (-79+n2 (169+2 n2 (-65+14 n2))))+n1 (-129-8 n2 (-60+n2 (158+n2 (-169+56 n2))))) Cos[(n1+n2) \[Pi]]-Csc[(n1+n2) \[Pi]] ((-60+3 n1 (43+4 n1 (-15+7 n1))+(209-8 n1 (17+(-2+n1) n1 (-43+56 n1))) n2+4 (101+2 n1 (-229+2 n1 (245+4 n1 (-43+7 n1)))) n2^2+4 (-345+2 n1 (463+4 n1 (-137+42 n1))) n2^3+64 (16+3 n1 (-12+7 n1)) n2^4+224 (-1+2 n1) n2^5) Sin[2 n1 \[Pi]]+(-60+n1 (209+4 n1 (101+n1 (-345+8 (32-7 n1) n1)))+129 n2+8 n1 (-17+n1 (-229+n1 (463+8 n1 (-36+7 n1)))) n2+4 (-45+4 n1 (-43+n1 (245-274 n1+84 n1^2))) n2^2+4 (21+2 n1 (155+8 n1 (-43+21 n1))) n2^3+448 (-1+n1) n1 n2^4) Sin[2 n2 \[Pi]])))/(28 n1 (-1+4 n1^2) n2 (-5+2 n1+2 n2) (-1+4 n2^2));*)


(* ::Input:: *)
(*(* checked... *)*)


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


(* ::Title:: *)
(*CHECK*)


(* ::Input:: *)
(*perm12=Permutations[{q1,q2}];*)
(*perm123=Permutations[{q1,q2,q3}];*)


(* ::Input:: *)
(*Clear[Fn,Gn,n,m];*)
(*Fn[n_,v_]=If[n==1,1,Sum[Gn[m,v[[1;;m]]]/((2n+3)(n-1)) ((2n+1)al[Total[v[[1;;m]]],Total[v[[m+1;;n]]]]Fn[n-m,v[[m+1;;n]]]+2be[Total[v[[1;;m]]],Total[v[[m+1;;n]]]]Gn[n-m,v[[m+1;;n]]]),{m,1,n-1}]];*)
(*Gn[n_,v_]=If[n==1,1,Sum[Gn[m,v[[1;;m]]]/((2n+3)(n-1)) (3al[Total[v[[1;;m]]],Total[v[[m+1;;n]]]]Fn[n-m,v[[m+1;;n]]]+2n be[Total[v[[1;;m]]],Total[v[[m+1;;n]]]]Gn[n-m,v[[m+1;;n]]]),{m,1,n-1}]];*)


(* ::Input:: *)
(*F2s[q1_,q2_]=Simplify[Sum[1/Length[perm12]Fn[2,{perm12[[i,1]],perm12[[i,2]]}],{i,1,Length[perm12]}]];*)
(*F3s[q1_,q2_,q3_]=Simplify[Sum[1/Length[perm123]Fn[3,{perm123[[i,1]],perm123[[i,2]],perm123[[i,3]]}],{i,1,Length[perm123]}]];*)


(* ::Input:: *)
(*alf[k1_,k2_]=(k1+k2) . k1/(k1 . k1);*)
(*bef[k1_,k2_]=(k1+k2) . (k1+k2)(k1 . k2)/2/(k1 . k1)/(k2 . k2) ;*)
(*alfm[k1_,k2_]=1+1/2 (mag[k1+k2]^2-mag[k1]^2-mag[k2]^2)/mag[k1]^2;*)
(*befm[k1_,k2_]=(mag[k1+k2]^2 (mag[k1+k2]^2-mag[k1]^2-mag[k2]^2))/(4mag[k1]^2 mag[k2]^2);*)
(*sigma2[k1_,k2_]=1/4 (mag[k1+k2]^2-mag[k1]^2-mag[k2]^2)^2/(mag[k1]^2 mag[k2]^2)-1;*)


(* ::Input:: *)
(*G2s[q1_,q2_]=Simplify[Sum[1/Length[perm12]Gn[2,{perm12[[i,1]],perm12[[i,2]]}],{i,1,Length[perm12]}]];*)


(* ::Input:: *)
(*(* 1/4 is for cosine squared in terms of norms... *)*)


(* ::Input:: *)
(*(* NG shape *)*)


(* ::Input:: *)
(*\[ScriptCapitalS][k1_,k2_,k3_]=(mag[k1]/mag[k2]+mag[k1]/mag[k3]+mag[k2]/mag[k1]+mag[k2]/mag[k3]+mag[k3]/mag[k1]+mag[k3]/mag[k2]-mag[k1]^2/(mag[k2]mag[k3])-mag[k2]^2/(mag[k1]mag[k3])-mag[k3]^2/(mag[k1]mag[k2])-2);*)


(* ::Section:: *)
(*(********)*)


(* ::Input:: *)
(*PNEWintegrand=(2G2s[q,k-q]\[ScriptCapitalS][q,k-q,k]/.{al->alfm,be->befm})/.{mag[0]->0,mag[-q]->q,mag[q]->q,mag[k]->k,mag[k-q]->kmq}//Expand*)


(* ::Input:: *)
(*(* check no k+q... *)*)


(* ::Input:: *)
(*PNEWintegrand=(2G2s[q,k-q]\[ScriptCapitalS][q,k-q,k]/.{al->alfm,be->befm})/.{mag[0]->0,mag[-q]->q,mag[q]->q,mag[k]->k,mag[k-q]->kmq,mag[k+q]->kmq}//Expand*)


(* ::Input:: *)
(*TabNEW=Table[{-(1/2)q D[PNEWintegrand[[i]] ,q]/PNEWintegrand[[i]],-(1/2)kmq D[PNEWintegrand[[i]] ,kmq]/PNEWintegrand[[i]],PNEWintegrand[[i]]/.{q->1,kmq->1,k->1}},{i,1,Length[PNEWintegrand]}]*)


(* ::Input:: *)
(*sum=Sum[TabNEW[[i, 3]] J2[1, TabNEW[[i, 1]] + n1, TabNEW[[i, 2]] + n2]/J2[1, n1, n2], {i, 1, Length[TabNEW]}] // FunctionExpand // FullSimplify*)


(* ::Input:: *)
(*PNEWabG2CHECK[n1_,n2_]=sum;*)


(* ::Input:: *)
(*Simplify[PNEWabG2CHECK[n1,n2]/PNEWabG2[n1,n2]]*)


(* ::Input:: *)
(*Simplify[PNEWabG2CHECK[n1,n2]-PNEWabG2CHECK[n2,n1]]*)


(* ::Input:: *)
(*(* of course... Also Subscript[F, 2] and Subscript[b, 2] and Subscript[\[ScriptCapitalG], 2] are obviously symmetric... *)*)


(* ::Section:: *)
(*(********)*)


(* ::Subsection:: *)
(*WRONG ATTEMPTS*)


(* ::Program:: *)
(*PNEWintegrand=(2*((1/2 (mag[k+q]^2-mag[k]^2-mag[q]^2)/((mag[k]mag[q])^0)(*careful!*))/mag[q]^2+(mag[k]^2-1/2 (mag[k+q]^2-mag[k]^2-mag[q]^2)/((mag[k]mag[q])^0)(*careful!*))/mag[k-q]^2)*\[ScriptCapitalS][q,k-q,k]/.{al->alfm,be->befm})/.{mag[0]->0,mag[-q]->q,mag[q]->q,mag[k]->k,mag[k-q]->kmq,mag[k+q]->kmq}//Expand*)


(* ::Program:: *)
(*4+(3 k^3)/kmq^3-(6 k^2)/kmq^2+(2 k)/kmq-(3 kmq)/k-k^3/q^3+k^4/(kmq q^3)-(2 k^2 kmq)/q^3+(2 k kmq^2)/q^3+kmq^3/q^3-kmq^4/(k q^3)+(2 k^2)/q^2-k^3/(kmq q^2)-(2 kmq^2)/q^2+kmq^3/(k q^2)-(6 k)/q-(3 k^4)/(kmq^3 q)+(3 k^3)/(kmq^2 q)+(4 k^2)/(kmq q)-kmq/q+(3 kmq^2)/(k q)-(3 q)/k+(2 k^2 q)/kmq^3+(4 k q)/kmq^2-q/kmq-(2 k q^2)/kmq^3-(2 q^2)/kmq^2+(3 q^2)/(k kmq)+q^3/kmq^3+q^3/(k kmq^2)-q^4/(k kmq^3)*)


(* ::Program:: *)
(*PNEWintegrand=(2*((1/2 (mag[k+q]^2-mag[k]^2-mag[q]^2)/((mag[k]mag[q])^0)(*careful!*))/mag[q]^2+(mag[k]^2-1/2 (mag[k+q]^2-mag[k]^2-mag[q]^2)/((mag[k]mag[q])^0)(*careful!*))/mag[k-q]^2)*\[ScriptCapitalS][q,k-q,k]/.{al->alfm,be->befm})/.{mag[0]->0,mag[-q]->q,mag[q]->q,mag[k]->k,mag[k-q]->kmq}//Expand*)


(* ::Program:: *)
(*2+(3 k^3)/kmq^3-(6 k^2)/kmq^2+(3 k)/kmq-kmq/k-k^3/q^3+k^4/(kmq q^3)-(k^2 kmq)/q^3+(k kmq^2)/q^3+(2 k^2)/q^2-k^3/(kmq q^2)-(k kmq)/q^2-(5 k)/q-(3 k^4)/(kmq^3 q)+(3 k^3)/(kmq^2 q)+(3 k^2)/(kmq q)-kmq/q+kmq^2/(k q)-(2 q)/k+(2 k^2 q)/kmq^3+(4 k q)/kmq^2-(2 k q^2)/kmq^3-(2 q^2)/kmq^2+(2 q^2)/(k kmq)+q^3/kmq^3+q^3/(k kmq^2)-q^4/(k kmq^3)-(k mag[k+q]^2)/kmq^3+(2 mag[k+q]^2)/kmq^2-(2 mag[k+q]^2)/(k kmq)+(k mag[k+q]^2)/q^3-(k^2 mag[k+q]^2)/(kmq q^3)+(kmq mag[k+q]^2)/q^3-(kmq^2 mag[k+q]^2)/(k q^3)-(2 mag[k+q]^2)/q^2+(k mag[k+q]^2)/(kmq q^2)+(kmq mag[k+q]^2)/(k q^2)+(2 mag[k+q]^2)/(k q)+(k^2 mag[k+q]^2)/(kmq^3 q)-(k mag[k+q]^2)/(kmq^2 q)-(q mag[k+q]^2)/kmq^3-(q mag[k+q]^2)/(k kmq^2)+(q^2 mag[k+q]^2)/(k kmq^3)*)


(* ::Program:: *)
(*(* NO, this does not work since I have BOTH k+q and k-q! I ran into this issue also when I was deriving the matrices on the old Mac... If I use this integrand, I do not match what I am using in CLASS-PT and I find an order of magnitude larger thing... Notice that in the other way I never run into this issue. Not with Subscript[G, 2], or Subscript[\[ScriptCapitalG], 2], or no-one else... Notice that I have also k-q in the NG shape!!! *)*)


(* ::Program:: *)
(*PNEWintegrand=(2*(-(-mag[k]^2+1/2 (mag[-q]^2-mag[k]^2-mag[k-q]^2))/mag[q]^2+(1/2 (mag[-q]^2-mag[k]^2-mag[k-q]^2))/mag[k-q]^2)*\[ScriptCapitalS][q,k-q,k]/.{al->alfm,be->befm})/.{mag[0]->0,mag[-q]->q,mag[q]->q,mag[k]->k,mag[k-q]->kmq,mag[k+q]->kmq}//Expand*)


(* ::Program:: *)
(*4-k^3/kmq^3+(2 k^2)/kmq^2-(6 k)/kmq-(3 kmq)/k+(3 k^3)/q^3-(3 k^4)/(kmq q^3)+(2 k^2 kmq)/q^3-(2 k kmq^2)/q^3+kmq^3/q^3-kmq^4/(k q^3)-(6 k^2)/q^2+(3 k^3)/(kmq q^2)+(4 k kmq)/q^2-(2 kmq^2)/q^2+kmq^3/(k q^2)+(2 k)/q+k^4/(kmq^3 q)-k^3/(kmq^2 q)+(4 k^2)/(kmq q)-kmq/q+(3 kmq^2)/(k q)-(3 q)/k-(2 k^2 q)/kmq^3-q/kmq+(2 k q^2)/kmq^3-(2 q^2)/kmq^2+(3 q^2)/(k kmq)+q^3/kmq^3+q^3/(k kmq^2)-q^4/(k kmq^3)*)


(* ::Program:: *)
(*PNEWintegrand=(2*(-(-mag[k]^2+1/2 (mag[-q]^2-mag[k]^2-mag[k-q]^2))/mag[q]^2+(1/2 (mag[-q]^2-mag[k]^2-mag[k-q]^2))/mag[k-q]^2)*\[ScriptCapitalS][q,k-q,k]/.{al->alfm,be->befm})/.{mag[0]->0,mag[-q]->q,mag[q]->q,mag[k]->k,mag[k-q]->kmq}//Expand*)


(* ::Program:: *)
(*4-k^3/kmq^3+(2 k^2)/kmq^2-(6 k)/kmq-(3 kmq)/k+(3 k^3)/q^3-(3 k^4)/(kmq q^3)+(2 k^2 kmq)/q^3-(2 k kmq^2)/q^3+kmq^3/q^3-kmq^4/(k q^3)-(6 k^2)/q^2+(3 k^3)/(kmq q^2)+(4 k kmq)/q^2-(2 kmq^2)/q^2+kmq^3/(k q^2)+(2 k)/q+k^4/(kmq^3 q)-k^3/(kmq^2 q)+(4 k^2)/(kmq q)-kmq/q+(3 kmq^2)/(k q)-(3 q)/k-(2 k^2 q)/kmq^3-q/kmq+(2 k q^2)/kmq^3-(2 q^2)/kmq^2+(3 q^2)/(k kmq)+q^3/kmq^3+q^3/(k kmq^2)-q^4/(k kmq^3)*)


(* ::Program:: *)
(*(* the second piece is k.(k-q) = ((|k+k-q(|^2)-|k(|^2)-|k-q(|^2))/2). In the first piece I have k.q which I rewrite as -(-k.q)=-(-k^2+k^2-k.q)=-(-k^2+k.(k-q))... *)*)


(* ::Program:: *)
(*(* NO these formulas are wrong... (|k+k-q(|^2)-|k(|^2)-|k-q(|^2))/2 is NOT k.(k-q)... I have mag[k-q]^2=k^2+q^2-2\[Times]k.q... FFS... *)*)


(* ::Subsection:: *)
(*CORRECT ONE...*)


(* ::Input:: *)
(*PNEWintegrand=(2*(((1/2 (mag[k]^2+mag[q]^2-mag[k-q]^2))/mag[q]^2+(mag[k]^2-(1/2 (mag[k]^2+mag[q]^2-mag[k-q]^2)))/mag[k-q]^2)*\[ScriptCapitalS][q,k-q,k])/.{al->alfm,be->befm})/.{mag[0]->0,mag[-q]->q,mag[q]->q,mag[k]->k,mag[k-q]->kmq}//Expand*)


(* ::Input:: *)
(*TabNEW=Table[{-(1/2)q D[PNEWintegrand[[i]] ,q]/PNEWintegrand[[i]],-(1/2)kmq D[PNEWintegrand[[i]] ,kmq]/PNEWintegrand[[i]],PNEWintegrand[[i]]/.{q->1,kmq->1,k->1}},{i,1,Length[PNEWintegrand]}]*)


(* ::Input:: *)
(*sum=Sum[TabNEW[[i, 3]] J2[1, TabNEW[[i, 1]] + n1, TabNEW[[i, 2]] + n2]/J2[1, n1, n2], {i, 1, Length[TabNEW]}] // FunctionExpand // FullSimplify*)


(* ::Input:: *)
(*(* put here to check, since I had to go through many iterations before actually fixing this... *)*)


(* ::Input:: *)
(*PNEWabb1\[ScriptF]SIMPLE[n1,n2]*)


(* ::Input:: *)
(*PNEWabb1\[ScriptF]SIMPLECHECK[n1_,n2_]=sum;*)


(* ::Input:: *)
(*Simplify[PNEWabb1\[ScriptF]SIMPLECHECK[n1,n2]/PNEWabb1\[ScriptF]SIMPLE[n1,n2]]*)


(* ::Input:: *)
(*Simplify[PNEWabb1\[ScriptF]SIMPLECHECK[n1,n2]-PNEWabb1\[ScriptF]SIMPLECHECK[n2,n1]]*)


(* ::Input:: *)
(*(* THE FULL MATRIX IS SYMMETRIC!!! *)*)


(* ::Section:: *)
(*(********)*)


(* ::Input:: *)
(*kdotq=1/2 (mag[k]^2+mag[q]^2-mag[k-q]^2);*)
(*qdotkmq=1/2 (mag[k]^2+mag[q]^2-mag[k-q]^2)-mag[q]^2;*)
(*kdotkmq=mag[k]^2-kdotq;*)


(* ::Input:: *)
(*\[ScriptCapitalA]prime=2/3*((qdotkmq*kdotkmq)/(mag[q]^2 mag[k-q]^2)+(qdotkmq*kdotq)/(mag[k-q]^2 mag[q]^2))+1/3*(kdotq/mag[q]^2+kdotkmq/mag[k-q]^2);*)


(* ::Input:: *)
(*\[ScriptCapitalB]prime=kdotq/mag[q]^2*kdotkmq^2/(mag[k-q]^2 mag[k]^2)+kdotkmq/mag[k-q]^2*kdotq^2/(mag[k]^2 mag[q]^2);*)


(* ::Subsection:: *)
(*\[Mu]^2*)


(* ::Input:: *)
(*PNEWintegrand=(2*(3/2*(\[ScriptCapitalA]prime-\[ScriptCapitalB]prime)*\[ScriptCapitalS][q,k-q,k])/.{al->alfm,be->befm})/.{mag[0]->0,mag[-q]->q,mag[q]->q,mag[k]->k,mag[k-q]->kmq}//Expand*)


(* ::Input:: *)
(*TabNEW=Table[{-(1/2)q D[PNEWintegrand[[i]] ,q]/PNEWintegrand[[i]],-(1/2)kmq D[PNEWintegrand[[i]] ,kmq]/PNEWintegrand[[i]],PNEWintegrand[[i]]/.{q->1,kmq->1,k->1}},{i,1,Length[PNEWintegrand]}]*)


(* ::Input:: *)
(*sum=Sum[TabNEW[[i, 3]] J2[1, TabNEW[[i, 1]] + n1, TabNEW[[i, 2]] + n2]/J2[1, n1, n2], {i, 1, Length[TabNEW]}] // FunctionExpand // FullSimplify*)


(* ::Input:: *)
(*PNEWab\[Mu]2SIMPLECHECK[n1_,n2_]=sum*)


(* ::Input:: *)
(*Simplify[PNEWab\[Mu]2SIMPLECHECK[n1,n2]/PNEWab\[Mu]2SIMPLE[n1,n2]]*)


(* ::Input:: *)
(*Simplify[PNEWab\[Mu]2SIMPLECHECK[n1,n2]-PNEWab\[Mu]2SIMPLECHECK[n2,n1]]*)


(* ::Input:: *)
(*(* THE FULL MATRIX IS SYMMETRIC!!! *)*)


(* ::Subsection:: *)
(*\[Mu]^4*)


(* ::Input:: *)
(*PNEWintegrand=(2*(-(3/2)*(\[ScriptCapitalA]prime-5/3*\[ScriptCapitalB]prime)*\[ScriptCapitalS][q,k-q,k])/.{al->alfm,be->befm})/.{mag[0]->0,mag[-q]->q,mag[q]->q,mag[k]->k,mag[k-q]->kmq}//Expand*)


(* ::Input:: *)
(*TabNEW=Table[{-(1/2)q D[PNEWintegrand[[i]] ,q]/PNEWintegrand[[i]],-(1/2)kmq D[PNEWintegrand[[i]] ,kmq]/PNEWintegrand[[i]],PNEWintegrand[[i]]/.{q->1,kmq->1,k->1}},{i,1,Length[PNEWintegrand]}]*)


(* ::Input:: *)
(*sum=Sum[TabNEW[[i, 3]] J2[1, TabNEW[[i, 1]] + n1, TabNEW[[i, 2]] + n2]/J2[1, n1, n2], {i, 1, Length[TabNEW]}] // FunctionExpand // FullSimplify*)


(* ::Input:: *)
(*PNEWab\[Mu]4SIMPLECHECK[n1_,n2_]=sum*)


(* ::Input:: *)
(*Simplify[PNEWab\[Mu]4SIMPLECHECK[n1,n2]/PNEWab\[Mu]4SIMPLE[n1,n2]]*)


(* ::Input:: *)
(*Simplify[PNEWab\[Mu]4SIMPLECHECK[n1,n2]-PNEWab\[Mu]4SIMPLECHECK[n2,n1]]*)


(* ::Input:: *)
(*(* THE FULL MATRIX IS SYMMETRIC!!! *)*)


(* ::Title:: *)
(*POWERS of \[Mu]*)


(* ::Input:: *)
(*(* if I look only at matter and multiply by Subscript[Z, 1] \[Rule] b1\[Mu]powers => everything is multiplied by Subscript[f, NL], k^3, and the NGparam, and \[ScriptCapitalM][k]/k^2... *)*)


(* ::DisplayFormula:: *)
(*Subscript[Z, 2]\[Superset] Subscript[b, 1]*Subscript[f, NL]*Subscript[F, 2][Subscript[k, 1],Subscript[k, 2]] -> k^3 as the dimension, \[Nu]=-0.8*)
(*\[Mu]^0: PNEWabb1F2 -> ALREADY CHECKED...*)


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
(*(* order here is [b1,\[Mu]] powers... *)*)


(* ::Input:: *)
(*PNEWab[n1_,n2_]:=(Csc[n1 \[Pi]] Csc[n2 \[Pi]] (-2 (60+448 n1^4 (-1+n2) n2-3 n2 (43+4 n2 (-15+7 n2))+4 n1^3 (-21+2 n2 (169+4 n2 (-65+28 n2)))+4 n1^2 (45+4 n2 (-79+n2 (169+2 n2 (-65+14 n2))))+n1 (-129-8 n2 (-60+n2 (158+n2 (-169+56 n2))))) Cos[(n1+n2) \[Pi]]-Csc[(n1+n2) \[Pi]] ((-60+3 n1 (43+4 n1 (-15+7 n1))+209 n2-8 n1 (17+(-2+n1) n1 (-43+56 n1)) n2+4 (101+2 n1 (-229+2 n1 (245+4 n1 (-43+7 n1)))) n2^2+4 (-345+2 n1 (463+4 n1 (-137+42 n1))) n2^3+64 (16+3 n1 (-12+7 n1)) n2^4+224 (-1+2 n1) n2^5) Sin[2 n1 \[Pi]]+(-60+n1 (209+4 n1 (101+n1 (-345+8 (32-7 n1) n1)))+129 n2+8 n1 (-17+n1 (-229+n1 (463+8 n1 (-36+7 n1)))) n2+4 (-45+4 n1 (-43+n1 (245-274 n1+84 n1^2))) n2^2+4 (21+2 n1 (155+8 n1 (-43+21 n1))) n2^3+448 (-1+n1) n1 n2^4) Sin[2 n2 \[Pi]])))/(28 n1 (-1+4 n1^2) n2 (-5+2 n1+2 n2) (-1+4 n2^2));*)


(* ::Input:: *)
(*(* this is the same as PNEWabb1F2. I define them w.r.t. this, since it allows me to save a matrix import for \[Mu] powers... *)*)


(* ::Input:: *)
(*Simplify[PNEWabb1F2[nu1,nu2]/PNEWab[nu1,nu2]]*)


(* ::Input:: *)
(*(*############################################*)*)
(*(*############################################*)*)
(*(*############################################*)*)


(* ::Input:: *)
(*b1\[Mu]powers[0,0]=0;*)
(*b1\[Mu]powers[0,2]=(f*PNEWabG2[nu1,nu2]+f^2/2*PNEWab\[Mu]2SIMPLE[nu1,nu2])/PNEWab[nu1,nu2];*)
(*b1\[Mu]powers[0,4]=(f^2/2*PNEWab\[Mu]4SIMPLE[nu1,nu2])/PNEWab[nu1,nu2];*)
(*b1\[Mu]powers[1,0]=PNEWabb1F2[nu1,nu2]/PNEWab[nu1,nu2];*)
(*b1\[Mu]powers[1,2]=(f/2*PNEWabb1\[ScriptF]SIMPLE[nu1,nu2])/PNEWab[nu1,nu2];*)
(*b1\[Mu]powers[1,4]=0;*)


(* ::Input:: *)
(*(* factors of 1/2 for the "RSD shift pieces" are correct... *)*)


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
(*(* see that I re-multiply by PNEWab... *)*)


(* ::Input:: *)
(*Do[AppendTo[total,{"vd2, f = "<>ToString[iii],vd2\[ScriptF][iii][nu1,nu2]}],{iii,whichpowers}];*)


(* ::Program:: *)
(*Do[Export[NotebookDirectory[] <> "vd2_f" <> ToString[iii] <> ".txt", vd2\[ScriptF][iii][nu1, nu2]];, {iii, whichpowers}]*)


(* ::Input:: *)
(*total*)


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


(* ::Program:: *)
(*Do[Export[NotebookDirectory[] <> "dd2_f" <> ToString[iii] <> ".txt", dd2\[ScriptF][iii][nu1, nu2]];, {iii, whichpowers}]*)


(* ::Input:: *)
(*total*)


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


(* ::Program:: *)
(*Do[Export[NotebookDirectory[] <> "vv4_f" <> ToString[iii] <> ".txt", vv4\[ScriptF][iii][nu1, nu2]];, {iii, whichpowers}]*)


(* ::Input:: *)
(*total*)


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


(* ::Program:: *)
(*Do[Export[NotebookDirectory[] <> "vd4_f" <> ToString[iii] <> ".txt", vd4\[ScriptF][iii][nu1, nu2]];, {iii, whichpowers}]*)


(* ::Input:: *)
(*total*)


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


(* ::Program:: *)
(*Do[Export[NotebookDirectory[] <> "vv6_f" <> ToString[iii] <> ".txt", vv6\[ScriptF][iii][nu1, nu2]];, {iii, whichpowers}]*)


(* ::Input:: *)
(*total*)


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
(*(* there is no need to check more after here... I have checked against the Mathematica implementation, which does NOT use all these extractions of powers of f and \[Mu], and it does everything itself... Likely I will do a check again for the multipoles, after AP and IR resummation... We'll see... For now I need to check against the brute-force integration... *)*)


(* ::Input:: *)
(*(*/////////////////////////////*)*)


(* ::Input:: *)
(*(* notice that I had screwed up because I was using the matrix M22basic_oneline_complex of Misha, that uses a different bias! The factors of 2 from the permutation were included in the matrix "reduction" formula, but for me I need two different biases for Subscript[\[ScriptCapitalG], 2] and Subscript[b, 2]... So I will need two matrices in real space... For the case of multipoles I will check what to do. If they were written in terms of \[Mu] powers it would have been better. I doubt that a great reduction in numbers will occur... We will see... *)*)


(* ::Input:: *)
(*(* WHATEVER THE COMMENT ABOVE MEANS, it is irrelevant... Now I import all matrices, computed by myself... For Subscript[b, 2], I have the correct bias of -1.25 in the bias script... Also, this comment is for some reason here, instead that in the notebooks for bias terms... Whatever... *)*)


(* ::Text:: *)
(*(openblas_test) gcabass@mbp156 ~/Documents/work/CLASS-PT-master-fNL/pt_matrices/compute_matrices_python/bias_real_space % more matrix_gen-bias_real_space.py*)
(*import numpy as np*)
(**)
(*#from whichdict import importdict*)
(**)
(*from sympy.parsing.mathematica import mathematica*)
(*from sympy import **)
(**)
(*from mpmath import **)
(**)
(*mp.dps = 32*)
(*mp.pretty = True*)
(**)
(*nu1 = var('nu1')*)
(*nu2 = var('nu2')*)
(**)
(*def J(nu1,nu2):*)
(*        return (gamma(1.5 - nu1) * gamma(1.5 - nu2) * gamma(nu1 + nu2 - 1.5) / (gamma(nu1) * gamma(nu2) * gamma(3. - nu2 - nu1))) / (8. * pi**(1.5))*)
(**)
(*#b=-0.8*)
(*#btab={-1.25,-0.8}*)
(*kmax = 1.e2*)
(*k0 = 0.00005*)
(**)
(*Nmaxtab=[128,256,512]*)
(**)
(*#namelist={"b2.txt","G2.txt"}*)
(**)
(*dictio={"b2.txt":-1.25,"G2.txt":-0.8}*)
(**)
(*for el in dictio:*)
(**)
(*    name = el*)
(*    b = dictio[el]*)


(* ::Text:: *)
(*(openblas_test) gcabass@mbp156 ~/Documents/work/CLASS-PT-master-fNL/pt_matrices/compute_matrices_python/bias_multipoles % more matrix_gen-bias_multipoles.py*)
(*import numpy as np*)
(**)
(*from whichdict import importdict*)
(**)
(*from sympy.parsing.mathematica import mathematica*)
(*from sympy import **)
(**)
(*from mpmath import **)
(**)
(*mp.dps = 32*)
(*mp.pretty = True*)
(**)
(*nu1 = var('nu1')*)
(*nu2 = var('nu2')*)
(**)
(*def J(nu1,nu2):*)
(*        return (gamma(1.5 - nu1) * gamma(1.5 - nu2) * gamma(nu1 + nu2 - 1.5) / (gamma(nu1) * gamma(nu2) * gamma(3. - nu2 - nu1))) / (8. * pi**(1.5))*)
(**)
(*kmax = 1.e2*)
(*k0 = 0.00005*)
(**)
(*Nmaxtab=[128,256,512]*)
(**)
(*for el in importdict:*)
(**)
(*    name = importdict[el][0]*)
(*    b = importdict[el][1]*)
(**)
(*    print("[][][][][][][][][][][][][][][][][][][][][][]")*)
(*    print(name)*)
(*    print("bias -> "+str(b))*)
(*    print("[][][][][][][][][][][][][][][][][][][][][][]\n")*)
(**)
(*    with open(name,"r") as file:*)
(*        expr = file.read()*)
(**)
(*    mathexpr=mathematica(expr)*)
(**)
(*    M12temp=lambdify([nu1,nu2],mathexpr,"mpmath")*)
(**)
(*(openblas_test) gcabass@mbp156 ~/Documents/work/CLASS-PT-master-fNL/pt_matrices/compute_matrices_python/bias_multipoles % more whichdict.py *)
(*#--> b2_vd0_f0.txt, b2_vv0_f1.txt, b2_vv2_f1.txt*)
(*#--> bG2_vd0_f0.txt, bG2_vv0_f1.txt, bG2_vv2_f1.txt*)
(**)
(*#-) vd0, f = 0 is [ 3 ] of vv0, f = 1*)
(*#-) vv2, f = 1 is [ 2 ] of vv0, f = 1*)
(**)
(*importdict={"1":["b2_vv0_f1.txt",-1.25],"2":["bG2_vv0_f1.txt",-0.8]}*)


(* ::CodeText:: *)
(*(openblas_test) gcabass@mbp156 ~/Documents/work/CLASS-PT-master-fNL/pt_matrices/compute_matrices_python/matter_mu_powers % more simplifications.txt*)
(*-) vv4, f = 3 is [ 1 ] of vd2, f = 2*)
(*(openblas_test) gcabass@mbp156 ~/Documents/work/CLASS-PT-master-fNL/pt_matrices/compute_matrices_python/matter_mu_powers % more whichdict.py *)
(*#--> vd2_f1.txt, vd2_f2.txt, dd2_f1.txt, vv4_f2.txt, vv4_f3.txt, vd4_f2.txt, vv6_f3.txt*)
(**)
(*#-) vv4, f = 3 is [ 1 ] of vd2, f = 2*)
(**)
(*importdict={"1":"vd2_f1.txt","2":"vd2_f2.txt","3":"dd2_f1.txt","4":"vv4_f2.txt","5":"vd4_f2.txt","6":"vv6_f3.txt"}*)


(* ::CodeText:: *)
(*\[Mu]^0*)
(*b1^0*)
(*vv0*)
(*0*)
(*b1^1*)
(*vd0*)
(*0*)
(*b1^2*)
(*dd0*)
(*1*)
(*\[Mu]^2*)
(*b1^0*)
(*vv2*)
(*0*)
(*b1^1*)
(*\[ScriptF]rules=CoefficientRules[vd2,f]//Simplify;*)
(*whichpowers=\[ScriptF]rules[[All,1]]//Flatten//Sort*)
(*{1,2}*)
(*b1^2*)
(*\[ScriptF]rules=CoefficientRules[dd2,f]//Simplify;*)
(*whichpowers=\[ScriptF]rules[[All,1]]//Flatten//Sort*)
(*{1}*)
(*\[Mu]^4*)
(*b1^0*)
(*\[ScriptF]rules=CoefficientRules[vv4,f]//Simplify;*)
(*whichpowers=\[ScriptF]rules[[All,1]]//Flatten//Sort*)
(*{2,3}*)
(*b1^1*)
(*\[ScriptF]rules=CoefficientRules[vd4,f]//Simplify;*)
(*whichpowers=\[ScriptF]rules[[All,1]]//Flatten//Sort*)
(*{2}*)
(*b1^2*)
(*dd4*)
(*0*)
(*\[Mu]^6*)
(*b1^0*)
(*\[ScriptF]rules=CoefficientRules[vv6,f]//Simplify;*)
(*whichpowers=\[ScriptF]rules[[All,1]]//Flatten//Sort*)
(*{3}*)
(*b1^1*)
(*vd6*)
(*0*)
(*b1^2*)
(*dd6*)
(*0*)
