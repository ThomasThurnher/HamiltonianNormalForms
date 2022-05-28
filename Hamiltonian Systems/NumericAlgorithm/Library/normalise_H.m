function [H_n] = normalise_H(H,m,P)
% This function performs normalisation of a Hamiltonian when given a
% generating function of order m

exp_order = numel(H);

H_n = repmat(polynomial_initialisation(),exp_order,1);
for i = 2:exp_order
    if i < m
        H_n(i) = H(i);
    else
        H_n(i) = normalise_Hi(H,i,m,P);
    end
end

end