(* Wolfram Language Package *)

BeginPackage["UsefulFunctions`"]

(* Exported symbols added here with SymbolName::usage *)

FullComplexExpand::usage = "FullComplexExpand[tensor] gives ComplexExpand of tensor to the infinite level
by default, level can be changed by `level -> n`"

EigenSystem::usage = "EigenSystem[m] returns {eigenvalues,eigenvectors} like System`Eigensystem.
It additionally does the following things:
- check reality of eigenvalues, disabled by `RealQ->False`
- the threshold of imaginary part is 10^-3, can be changed, e.g., by `ImValue->10^-1`
- eigenvalues are sorted in ascending order based on the real part, can be changed to based on the absolute value by `SortByAbs->True`
- eigenvectors are normalized wrt the metric 1";

EigenB::usage = "EigenB[m] returns {eval,evec} like System`Eigensystem.

";

ExpVal::usage = "ExpVal[vec,op] returns the expectation value of an operator (op) wrt the state (vec)";

Trunc::usage = "Trunc[m,d] returns the matrix by deleting d levels in top, bottom, left and right side of m";

AnniOpB::usage = "AnniOpB[n0,dMax] returns the bosonic annihilation operator in the Fock basis, centered at n0, starting (ending) at n0 + (-) dMax ";

CreaOpB::usage = "CreaOpB[n0,dMax] returns the bosonic creation operator in the Fock basis, centered at n0, starting (ending) at n0 + (-) dMax ";

NumbOpB::usage = "NumbOpB[n0,dMax] returns the bosonic number operator in the Fock basis, centered at n0, starting (ending) at n0 + (-) dMax ";



Begin["`Private`"] (* Begin Private Context *) 

FullComplexExpand[m_,OptionsPattern[Level->Infinity]] :=
    Map[ComplexExpand,m,OptionValue[Level]];

EigenSystem[m_?MatrixQ,OptionsPattern[{RealQ->True,ImValue->10^-3,SortByAbs->False}]] :=
    Module[ {eval,evec},
        If[ OptionValue[SortByAbs],
            {eval,evec} = Transpose@Sort[Transpose@Eigensystem@N@m,Abs[#[[1]]]<Abs[#2[[1]]]&],
            {eval,evec} = Transpose@Sort[Transpose@Eigensystem@N@m,Re[#[[1]]]<Re[#2[[1]]]&]
        ];
        If[ OptionValue[RealQ],
            If[ Max[Abs@Im@eval]>OptionValue[ImValue],
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
    
EigenB[m_,OptionsPattern[{RealQ->True,ImValue->10^-3}]] :=
    Module[ {group, dim, list, g, mG, eval1, evec1, mF, eval2, 
      evec2, eval3, evec3, eval, evec},
        group[l_List] :=
            Map[{#, 1} &, 
              l] //. {head___, {x_, n1_}, {x_, n2_}, 
               tail___} :> {head, {x, n1 + n2}, tail};
        dim = Length@m;
        If[ OddQ[m],
            Print["Wrong dimensions of the matrix !"];
            Abort[]
        ];
        {eval, evec} = EigenSystem[m, RealQ->OptionValue[RealQ],ImValue->OptionValue[ImValue]];
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
           {eval1, evec1} = EigenSystem[mG, RealQ->OptionValue[RealQ],ImValue->OptionValue[ImValue]];
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
        (* try to make `evec` real *)
        Do[evec[[l]] = 
          evec[[l]]*
           E^(-I Arg[
             SelectFirst[evec[[l]], Abs[#] >= Max[Abs[evec[[l]]]] &]]),
                 {l, Length@evec}];
        {eval, evec}
    ];

ExpVal[vec_?VectorQ,op_?MatrixQ] :=
    Conjugate[vec].op.vec;
    
Trunc[m_?MatrixQ,d_?(Positive[#]&&IntegerQ[#]&)] :=
    Delete[Delete[
    Delete[Delete[m, d]\[Transpose], 
    d]\[Transpose], -d]\[Transpose], -d]\[Transpose];

AnniOpB[n0_?NonNegative, dMax_?Positive] :=
    SparseArray[{i_, j_} /; i == j - 1 :> Sqrt[
      n0 - dMax + j - 1], {2 dMax + 1, 2 dMax + 1}];
    
CreaOpB[n0_?NonNegative, dMax_?Positive] :=
    SparseArray[{i_, j_} /; i == j + 1 :> Sqrt[
      n0 - dMax + j], {2 dMax + 1, 2 dMax + 1}];
    
NumbOpB[n0_?NonNegative, dMax_?Positive] :=
    SparseArray[{i_, j_} /; i == j :> n0 - dMax + j - 1, {2 dMax + 1, 
      2 dMax + 1}];

End[] (* End Private Context *)

EndPackage[]