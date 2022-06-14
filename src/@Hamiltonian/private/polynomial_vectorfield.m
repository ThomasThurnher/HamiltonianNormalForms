function [JDP] = polynomial_vectorfield(P)
% This function computes the hamiltonian vector field produced by the homogenous
% polynomial function P
% It stores the i-th component of the vectorfield in its ith entry

if isempty(P) || isempty(P.coeffs)
    JDP = polynomial_initialisation();
    return
end
n   = size(P.ind,1)/2;
DP  = polynomial_gradient(P);

JDP = J_DP(DP,n);

end


function[JDP] = J_DP(DP,n)
% Transforms DP with complex structure
DPc = cell(n,1);
for i = 1:n
    DPc{i} = -DP(i).coeffs;
end

[DP(1:n).coeffs] = deal(DPc{:});

JDP = repmat(polynomial_initialisation(),2*n,1);

JDP(1:n) = DP(n+1:2*n);
JDP(n+1:2*n) = DP(1:n);
end