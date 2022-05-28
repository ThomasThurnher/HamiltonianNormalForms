function [PmF] = polynomial_addition(P,F)
% This function subtracts two homogenous polynomials P - F
tol = 1e-10;

Fcoeffs =   [F.coeffs];
Find    =   [F.ind];
Pcoeffs =   [P.coeffs];
Pind    =   [P.ind];

PmF.coeffs = [Pcoeffs,Fcoeffs];
PmF.ind    = [Pind   ,Find   ];

PmF = polynomial_collapse(PmF,tol);
end