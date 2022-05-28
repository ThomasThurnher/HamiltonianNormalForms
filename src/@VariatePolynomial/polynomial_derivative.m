function [DjP] = polynomial_derivative(obj,P,j)
% This function computes the derivative dP/Dx_j

if isempty(P.coeffs)
    DjP = polynomial_initialisation();
    return
end

n = size(P.ind,1)/2;
if j > 2*n
    error('Derivative direction exceeds dimension')
end


P_ind    = [P.ind];
P_coeffs = [P.coeffs];


ej = ((1:2*n).' == j);
[DjPind, idx] = multi_subtraction(P_ind,ej,'Derivative');

% Get coefficients
Dcoeff = P_coeffs(idx) .* P_ind(j,idx);

DjP.coeffs = Dcoeff;
DjP.ind    = DjPind;

DjP = polynomial_collapse(DjP);
end