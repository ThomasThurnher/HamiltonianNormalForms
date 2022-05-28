function [p_revlex] = polynomial_2_revlex(p)
% Reorders coefficients and indices of homogenous polynomial to be in revlex
% ordering. The linear part of phi is the identity which already is in
% revlex ordering.

num = numel(p);
p_revlex = polynomial_initialisation();
% order of polynomial
i = sum(p.ind(:,1));
p_revlex = reorder_polynomial(p,i);

end


function [out] = reorder_polynomial(in,i)

N = size(in.ind,1);
z_i = nchoosek(i+N-1,N-1);
% Positions of indices in revlex ordered set, -> col indices
[I] = multi_index_2_ordering(in.ind,'revlex',[]);

outind      = flip(sortrows(nsumk(N,i,'nonnegative')).',2);

outcoeffs   = zeros(1,z_i);
outcoeffs(I) = in.coeffs;

out.ind = outind;
out.coeffs = outcoeffs;
end