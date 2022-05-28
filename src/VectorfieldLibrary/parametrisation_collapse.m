function [Pcoll] = parametrisation_collapse(P, varargin)
% Collapsesparametrisation to only contain terms of unique monomials
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

[unInd, ~, idx] = unique(P.ind.' , 'rows');

N = size(P.coeffs,1);
coeffs = zeros(N, size(unInd,1)); % unInd is transposed

if ~isempty(idx)
    for i = 1:N
    coeffs(i,:) = accumarray(idx, P.coeffs(i,:)).'; % add coefficients
    end
    [~,col] = find(abs(coeffs) > tol);

    keeper = unique(col);
    Pcoll.coeffs = coeffs(:,keeper);
    Pcoll.ind = unInd(keeper,:).'; % unique monomials
else
    % The polynomial is empty
    Pcoll = polynomial_initialisation();
end
end