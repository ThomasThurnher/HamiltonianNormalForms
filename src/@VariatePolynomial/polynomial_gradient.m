function [DP] = polynomial_gradient(obj,P)
% This function computes the gradient of a homogenous polynomial function

if isempty(P.coeffs)
    DP = polynomial_initialisation();
    return
end

n = size(P.ind,1)/2;

DP = repmat(polynomial_initialisation(),2*n,1); 

for i = 1:2*n
    DP(i) = polynomial_derivative(P,i);
end

end
