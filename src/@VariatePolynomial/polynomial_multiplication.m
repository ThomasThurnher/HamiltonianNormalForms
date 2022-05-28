function[PF] = polynomial_multiplication(obj,P,F)
% Multiplication of two homogenous polynomials
tol = 1e-10;

if isempty(P.coeffs) || isempty(F.coeffs)
    PF = polynomial_initialisation();
    return 
end

N = size(P.ind,1);

ind = combvec([P.ind],[F.ind]);
ind = ind(1:N,:) + ind(N+1:end,:); % added multi_indices of all possible combinations


coeffs = F.coeffs.'*P.coeffs;
coeffs = reshape(coeffs.',1,[]);

PF.ind = ind;
PF.coeffs = coeffs;

PF = polynomial_collapse(PF, tol);
end