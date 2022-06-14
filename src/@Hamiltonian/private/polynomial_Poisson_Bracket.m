function [PB] = polynomial_Poisson_Bracket(P,F)
% Computes the Poisson Bracket of two homogenous polynomials P and F which
% both take inputs x_1 , ... , x_n, y_1, .... y_n.
tol = 1e-10;
n = size(P.ind,1)/2;
PB_j = polynomial_initialisation();


for j = 1:n
    DPDxi = polynomial_derivative(P,j);
    DPDyi = polynomial_derivative(P,n+j);
    DFDxi = polynomial_derivative(F,j);
    DFDyi = polynomial_derivative(F,n+j);
       
    tmp1  = polynomial_multiplication(DPDxi,DFDyi);
    tmp2  = polynomial_multiplication(DFDxi,DPDyi);

    PB_j(j) = polynomial_subtraction(tmp1,tmp2);
end

PB.coeffs = [PB_j.coeffs];
PB.ind    = [PB_j.ind];

PB = polynomial_collapse(PB,tol);
end