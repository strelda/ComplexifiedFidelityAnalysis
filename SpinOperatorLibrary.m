(* ::Package:: *)

(* ::Package:: *)
(**)


BeginPackage["SpinOperatorLibrary`"];

(* Exported symbols *)
Je::usage = "Je[sn, j, m] calculates the matrix element for spin operators.";
Jx::usage = "Jx[n] generates the Jx matrix for dimension n.";
Jx2::usage = "Jx2[n] generates the \!\(\*SuperscriptBox[\(Jx\), \(2\)]\) matrix for dimension n.";
Jy::usage = "Jy[n] generates the Jy matrix for dimension n.";
Jz::usage = "Jz[n] generates the Jz matrix for dimension n.";

Begin["`Private`"]

Je[sn_, j_, m_] := Sqrt[j*(j + 1) - m*(m + sn*1)]
 
Jx[(n_)?IntegerQ] := Module[{diag, b1, b2, m}, b1 = Je[1, n/2, m]/2; 
      SparseArray[{Band[{1, 2}, {n, n + 1}] -> Table[b1, {m, n/2 - 1, -n/2, 
           -1}], Band[{2, 1}, {n + 1, n}] -> Table[Conjugate[b1], 
          {m, n/2 - 1, -n/2, -1}]}, {n + 1, n + 1}]] 
 
Jy[(n_)?IntegerQ] := Module[{diag, b1, b2, m}, b1 = Je[1, n/2, m]/(2*I); 
      SparseArray[{Band[{1, 2}, {n, n + 1}] -> Table[b1, {m, n/2 - 1, -n/2, 
           -1}], Band[{2, 1}, {n + 1, n}] -> Table[Conjugate[b1], 
          {m, n/2 - 1, -n/2, -1}]}, {n + 1, n + 1}]]
 
Jz[(n_)?IntegerQ] := Module[{diag, b1, b2, m}, 
     diag = m; SparseArray[{Band[{1, 1}, {n + 1, n + 1}] -> 
         Table[diag, {m, n/2, -n/2, -1}]}, {n + 1, n + 1}]]
         
Jx2[(n_)?IntegerQ] := Module[{diag,b1,b2,m},
diag=1/4 (Je[-1,n/2,m]Je[1,n/2,m-1]+Je[1,n/2,m]Je[-1,n/2,m+1]);(*diagonal*)
b2=1/4 ( Je[1,n/2,m+1]Je[1,n/2,m]);(*2. upper offdiagonal*)

SparseArray[{
Band[{1,1},{n+1,n+1}]->Table[diag,{m,n/2,-n/2,-1}](*diagonal*),
Band[{1,3},{n-1,n+1}]->Table[b2,{m,n/2-2,-n/2,-1}] (*2. upper offdiagonal*),Band[{3,1},{n+1,n-1}]->Table[b2,{m,n/2-2,-n/2,-1}] (*2. lower offdiagonal*)
},{n+1,n+1}]
]

End[];
EndPackage[];
