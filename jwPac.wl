(* ::Package:: *)

(* Wolfram Language Package *)

BeginPackage["jwPackage`"]
(* Exported symbols added here with SymbolName::usage *)  

FullComplexExpand::usage = "FullComplexExpand[tensor] gives ComplexExpand to the infinite level
by default. One can also use the option level -> n to specify the level";

Eigen::usage = "Eigen[m] gives eigenvalues and eigenvectors like Eigensystem,
but it also checks reality of eigenvalues and normalizes eigenvectors.";

EigenB::usage = "EigenB[m] gives eigenvalues and eigenvectors like Eigen, but for a bosonic BdG matrix.
Eigenvalues are sorted in the order like 1,2,3,...,-1,-2,-3,...
Eigenvectors are normalized even in the presence of degeneracy";

VOISimplify::usage = "VOISimplify[vars_,expr_,assum_] simplifies expr using assum.";

Begin["`Private`"] (* Begin Private Context *) 

FullComplexExpand[m_,level_:Infinity]:=Map[ComplexExpand,m,level];

Eigen[m_,OptionsPattern[{check->True,small->10^-3}]] :=
    Module[ {eval,evec},
        {eval,evec} = Transpose@Sort[Transpose@Eigensystem@N@m,Re[#[[1]]]<Re[#2[[1]]]&];
        If[ OptionValue[check],
            If[ Max[Abs@Im@eval]>OptionValue[small],
                Print["Instable! eval=",
                Chop@eval];
                Abort[]
            ]
        ];
        evec = Normalize/@evec;
        (* (try to) fix the phase *)
        Do[evec[[l]] = evec[[l]]*E^(-I Arg[SelectFirst[evec[[l]],Abs[#]>=Max[Abs[evec[[l]]]]&]]),{l,Length@evec}];
        {eval,evec}
    ];
 
EigenB[m_,OptionsPattern[{check->True,small->10^-3}]] :=
    Module[ {group, dim, list, g, mG, eval1, evec1, mF, eval2, 
      evec2, eval3, evec3, eval, evec},
        group[l_List] :=
            Map[{#, 1} &, 
              l] //. {head___, {x_, n1_}, {x_, n2_}, 
               tail___} :> {head, {x, n1 + n2}, tail};
        dim = Length@m;
        {eval, evec} = Eigen[m, check->OptionValue[check],small->OptionValue[small]];
        (* there is a very tricky point that I ignored before, 
        group eigenvectors based on the REAL part of the eigenvalue !!! *)
        list = 
        TakeList[evec, group[Round[Re[eval], 10^-12]]\[Transpose][[2]]];
        Do[g = 
          If[ If[ i == 1,
                  0,
                  Length[Flatten[list[[1 ;; i - 1]], 1]]
              ] < dim/2,
              1,
              -1
          ];
           mG = Inverse[
             list[[i]]\[Conjugate] . 
              KroneckerProduct[PauliMatrix[3], IdentityMatrix[dim/2]] . 
              list[[i]]\[Transpose]];
           {eval1, evec1} = Eigen[mG,  check->OptionValue[check],small->OptionValue[small]];
           mF = evec1\[Transpose] . DiagonalMatrix[Sqrt[g*eval1]] . 
             Inverse[evec1\[Transpose]];
           list[[i]] = (list[[i]]\[Transpose] . mF)\[Transpose],
               {i, Length@list}];
        evec = Flatten[list, 1];
        {eval2, evec2} = 
         Transpose@
          Sort[Transpose@{eval, evec}, Re[#1[[1]]] < Re[#2[[1]]] &];
        {eval3, evec3} = 
         Transpose@
          Sort[Transpose@{eval, evec}, Re[#1[[1]]] > Re[#2[[1]]] &];
        eval = Join[eval2[[dim/2 + 1 ;; dim]], eval3[[dim/2 + 1 ;; dim]]];
        evec = Join[evec2[[dim/2 + 1 ;; dim]], evec3[[dim/2 + 1 ;; dim]]];
        (* try to make `evec` real, possibly fail *)
        Do[evec[[l]] = 
          evec[[l]]*
           E^(-I Arg[
             SelectFirst[evec[[l]], Abs[#] >= Max[Abs[evec[[l]]]] &]]),
                 {l, Length@evec}];
        {eval, evec}
    ];

VOISimplify[vars_,expr_,assum_:True]:=Module[{perm,ee,best},perm=Permutations[vars];
ee=(FullSimplify@@({expr,assum}/.Thread[vars->#]))&/@perm;
best=Sort[Transpose[{LeafCount/@ee,ee,perm}]][[1]];
best[[2]]/.Thread[best[[3]]->vars]];

End[] (* End Private Context *)
EndPackage[]



