function [phi] = hamiltonian_flow(X_H,m,order,n,varargin)
% This function computes the time 1 flow induced by the hamiltonian vector field
% X_H given as phi = exp(X_H)
% X_H   - Hamiltonian Vector field JDH
% m     - Order of the homogenous polynomial H
% order - order of polynomials up to which the flow should be computed

if isempty(varargin)
    reduction = false;
else
    reduction = varargin{1};
    'Reducing over first modepair'
end
i = 1;
VF_order = m-1; % Polynomial order of current vector field coefficients


% Flow Phi consists of nonlinear functions
phi = repmat(polynomial_initialisation(),order,1);
X_H_i = X_H;

% Near identity character of transformation
tmp.coeffs = eye(2*n);
tmp.ind    = eye(2*n);
phi(1) = tmp;

while VF_order <= order

    tmp = vectorfield2parametrisation(X_H_i,n);
    W.coeffs = 1/factorial(i) * tmp.coeffs;
    W.ind    = tmp.ind;
    
    if reduction
    W = reduce_W(W);
    end
    
    phi(VF_order) = W;
    VF_order = (m-1) + (VF_order-1);
    X_H_i = vectorfield_multiplication(X_H,X_H_i,2*n);
    
    i = i+1;
end

phi = flow_2_revlex(phi);
end

function [Wout] = reduce_W(Win)
if ~isempty(Win.ind)
ind = Win.ind([2,4],:);

%find all zero entries
idx = find( sum(ind,1) == 0);

Wout.coeffs = Win.coeffs(:,idx);
Wout.ind = Win.ind(:,idx);
else
    Wout = polynomial_initialisation();
end
end