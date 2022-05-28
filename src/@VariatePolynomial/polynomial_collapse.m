function [Pcoll] = polynomial_collapse(obj,P, varargin)
% Collapses polynomial to only contain terms of unique monomials
% where any coefficient smaller than tol is discarded

if isempty(varargin)
    tol = 1e-10;
else
    tol = varargin{1};
end

if isempty([P.coeffs]) 
    Pcoll = polynomial_initialisation();
    return
end

[unInd, ~, idx] = unique([P.ind].' , 'rows');

if ~isempty(idx)
    coeffs = accumarray(idx, [P.coeffs]).'; % add coefficients
    
    keeper = abs(coeffs) > tol;
    
    Pcoll.coeffs = coeffs(keeper);
    Pcoll.ind = unInd(keeper,:).'; % unique monomials
    
    % Output in reverse lexicographic ordering
    Pcoll.coeffs = flip(Pcoll.coeffs,2);
    Pcoll.ind    = flip(Pcoll.ind ,2 );
else
    % The polynomial is empty
    Pcoll = polynomial_initialisation();
end
end