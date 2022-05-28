function [phi_A] = compose_linear_flow_ssm(A,phi,type)
% This function computes the composition of the nonlinear flow phi with the
% linear map A
tol = 1e-5;

switch type
    case 'right'
        % Phi (A)
        % In this case, every map in phi will keep the same order
        
        phi_A = repmat(polynomial_initialisation(),numel(phi),1);
        N = size(A,1);
        H{1} = A;
        
        % Linear composition
        A_param.coeffs = A;
        A_param.ind = eye(N); % revlex order
        
        phi_A(1) = A_param;

        for k = 2:numel(phi)
            K = flip(sortrows(nsumk(N,k,'nonnegative')).',2);
            phi_k = phi(k)
            
            H_k  = linear_composition_coefficients(A,H,k); % Dependent on ordering, chooses conj, or revlex computation
            H{k} = H_k;
            if ~isempty(phi_k.coeffs)
                tmp.coeffs = phi_k.coeffs* compute_pi_hnf(phi_k.ind,K,H);
                tmp.ind    = K;
                
                phi_A(k)=tmp;
            end
        end
    case 'left'
        % compute A(phi)
        phi_A = repmat(polynomial_initialisation(),numel(phi),1);

        if ~isempty(phi(1)) && ~isempty(phi(1).coeffs)
        A_param.ind = phi(1).ind;
        A_param.coeffs = A* phi(1).coeffs;
        phi_A(1) = A_param;
        end
        for k = 2:numel(phi)
            phi_k = phi(k);


            if ~isempty(phi_k) && ~isempty(phi_k.coeffs)
                tmp.coeffs = A *phi_k.coeffs;
                tmp.ind    = phi_k.ind;
                
                phi_A(k)=tmp;
            end
        end

end
%{
for k = 2:numel(phi_A)
    if ~isempty(phi_A{k}) && ~isempty(phi_A{k}.coeffs)
    phi_A{k} = parametrisation_collapse(phi_A{k},tol);
    end
end
%}
end