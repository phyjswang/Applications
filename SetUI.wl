(* ::Package:: *)

BeginPackage["SetUI`"]


SetUI::usage="
value `0`: minimal UI
value `1`: default UI
all others: returns error info
";


Begin["`Private`"]


SetUI[1]:=SetOptions[$FrontEnd,WindowElements->(WindowElements/. Options[$DefaultFrontEnd,WindowElements])];
SetUI[0]:=SetOptions[$FrontEnd, WindowElements -> {"VerticalScrollBar"}];
SetUI[x_]:=Print["only 1 or 0 are allowed: default or clean UI"]


End[]


EndPackage[]
