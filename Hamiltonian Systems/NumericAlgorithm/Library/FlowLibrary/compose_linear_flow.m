function [phi] = compose_linear_flow(A,phi)

for k = 1:numel(phi)
if ~isempty(phi(k)) && ~isempty(phi(k).coeffs)
 phi(k).coeffs = A*phi(k).coeffs;
end
end


end
