function [Param] = vectorfield2parametrisation(VF,n)
% This function converts a vectorfield VF with homogenous polynomials of
% the same order as coefficients into the representation used to store SSM
% parametrisation coefficients.

if isempty([VF.coeffs])
    Param = polynomial_initialisation(); 
    return
end


[un_ind,~ ,col_idx] = unique([VF.ind].','rows');
un_ind = un_ind.';

row_idx = [];
for i = 1:2*n
    if ~isempty(VF(i).coeffs)
    num = numel(VF(i).coeffs);
    row_idx = [row_idx, i * ones(1,num)];
    end
    
end

coeffs = [VF.coeffs];
coeffs = full(sparse(row_idx,col_idx,coeffs,2*n,size(un_ind,2)));

ind   = un_ind;

Param.coeffs = coeffs;
Param.ind = ind;
end