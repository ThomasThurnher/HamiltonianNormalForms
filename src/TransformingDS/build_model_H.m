function [A,B,F] = build_model(H,N)
% Given a Hamiltonian, this function creates the corresponding nonlinear
% autonomous dynamical system

B = eye(N);

VF_lin = polynomial_vectorfield(H(2));

A = zeros(N);
for i = 1:N
    [row, ~] = find(VF_lin(i).ind);
    A(i,row) = VF_lin(i).coeffs;
end


for i = 3:numel(H)
    VFi = polynomial_vectorfield(H(i));
    F(i-1) = vectorfield2parametrisation(VFi,N/2);
    F(i-1).ind = F(i-1).ind.';
end

end