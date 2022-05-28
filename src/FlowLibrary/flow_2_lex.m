function [phi_revlex] = flow_2_lex(phi,varargin)
% Reorders coefficients and indices in every entry of phi to be in revlex
% ordering. 

num = numel(phi);
phi_revlex = repmat(polynomial_initialisation(),num,1);

%{
if isdiag(phi(1).ind)
    phi_revlex(1).ind = flip(phi(1).ind,2);
    phi_revlex(1).coeffs = flip(phi(1).coeffs,2);
else
    phi_revlex(1).ind = phi(1).ind;
    phi_revlex(1).coeffs = phi(1).coeffs;
end
%}
for i = 1:num
    if ~isempty(phi(i)) && ~isempty(phi(i).coeffs)
            phi_revlex(i) = reorder_flow(phi(i));
    end
end
end


function [out] = reorder_flow(in)
N = size(in.coeffs,1);
l = size(in.ind,1);

i = sum(in.ind(:,1));
z_i = nchoosek(i+l-1,l-1);
% Positions of indices in revlex ordered set, -> col indices
[I] = multi_index_2_ordering(in.ind,'lex',[]);

outind      = sortrows(nsumk(l,i,'nonnegative')).';

outcoeffs   = zeros(N,z_i);
outcoeffs(:,I) = in.coeffs;

out.ind = outind;
out.coeffs = outcoeffs;

end

