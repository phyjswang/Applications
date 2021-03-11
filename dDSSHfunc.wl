(* ::Package:: *)

(* Get["/Users/jw-mbp/Documents/Wolfram Mathematica/jwPac.wl"] *)
(* <<PauliAlgebra` *)
equalQ[m1_,m2_]:=Which[
Norm[m1]>10^-8,Norm[m1-m2]/Norm[m1]<10^-3,
Norm[m2]>10^-8,Norm[m1-m2]/Norm[m2]<10^-3,
True,True];
SpinQ[S_]:=IntegerQ[2S]&&S>=0;
splus[0]={{0}}//SparseArray;
splus[S_?SpinQ]:=splus[S]= SparseArray[Band[{1,2}]->Table[Sqrt[S(S+1)-M(M+1)],{M,S-1,-S,-1}],{2S+1,2S+1}];
sminus[S_?SpinQ]:=Transpose[splus[S]];
sx[S_?SpinQ]:=sx[S]=(splus[S]+sminus[S])/2;
sy[S_?SpinQ]:=sy[S]=(splus[S]-sminus[S])/(2I);
sz[S_?SpinQ]:=sz[S]=SparseArray[Band[{1,1}]->Range[S,-S,-1],{2S+1,2S+1}];
id[S_?SpinQ]:=id[S]=IdentityMatrix[2S+1,SparseArray];
s[1]=sx[1];
s[2]=sy[1];
s[3]=sz[1];
sp=splus[1];
sm=sminus[1];
ket={E^(I \[Eta]/2) Sin[\[Theta]/2]E^(-I \[Phi]/2) Cos[\[Chi]/2],Cos[\[Theta]/2],E^(I \[Eta]/2) Sin[\[Theta]/2]E^(I \[Phi]/2) Sin[\[Chi]/2]};
expVal[o_,k_]:=k\[Conjugate] . o . k;
(*
f[\[Alpha]_,\[Beta]_,\[Gamma]_,\[Sigma]_]:=Sum[tm\[ConjugateTranspose]\[LeftDoubleBracket]\[Alpha],;;\[RightDoubleBracket].s[\[Mu]].tm\[LeftDoubleBracket];;,\[Beta]\[RightDoubleBracket]*tm\[ConjugateTranspose]\[LeftDoubleBracket]\[Gamma],;;\[RightDoubleBracket].s[\[Mu]].tm\[LeftDoubleBracket];;,\[Sigma]\[RightDoubleBracket],{\[Mu],2}];
(tm={{Sin[\[Theta]/2]/Sqrt[2],Cos[\[Theta]/2],Sin[\[Theta]/2]/Sqrt[2]},{-(Cos[\[Theta]/2]/Sqrt[2]),Sin[\[Theta]/2],-(Cos[\[Theta]/2]/Sqrt[2])},{1/Sqrt[2],0,-1/Sqrt[2]}}\[Transpose]);
*)
Clear[mt];
mt[1,1]=SparseArray[{1,2}->-t2,{2,2}];
mt[1,0]=SparseArray[{2,1}->-t1,{2,2}]+SparseArray[{2,1}->-t1,{2,2}]\[ConjugateTranspose];
mt[2,1]=SparseArray[{{1,2}->-t2,{4,3}->-t2},{4,4}];
mt[2,2]=SparseArray[{{2,3}->-t2,{1,4}->-t2},{4,4}];
mt[2,0]=SparseArray[{{2,1}->-t1,{3,2}->-t1,{4,3}->-t1,{1,4}->-t1},{4,4}]+
SparseArray[{{2,1}->-t1,{3,2}->-t1,{4,3}->-t1,{1,4}->-t1},{4,4}]\[ConjugateTranspose];
