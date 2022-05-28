function [phi12] = compose_flow(phi1,phi2,N,order,reduction,resModes)
% This function computes the composition of two nonlinear flows. This
% implies that both of them have unity linear flow and possibly nonlinear
% contributions.
%
% order  - highest order in which polynomials are to be kept
% N      - phase space size


% l      - multi-index size

phi12 = repmat(polynomial_initialisation(),order,1);

% Linear composition

phi12(1) = phi1(1);
H{1} = phi2(1).coeffs;

for k = 2:order % order of phi12
    %if reduction
     %   K = flip(sortrows(nsumk(2,k,'nonnegative')).',2);
    %else
        K = flip(sortrows(nsumk(N,k,'nonnegative')).',2);
    %end
    
    H_k  = coeffs_composition_hnf(phi2,H,k); % Dependent on ordering, chooses conj, or revlex computation
    H{k} = H_k;
    
    tmp.coeffs = zeros(N,size(K,2));
    
    for n = 2:k % orders of phi1
        phi_n = phi1(n);
        
        if ~isempty(phi_n.coeffs)

            tmp.coeffs = tmp.coeffs + phi_n.coeffs* compute_pi_hnf(phi_n.ind,K,H,N,resModes);
        end
    end
    tmp.ind    = K;
    
    % Add contributions due to the linear term in phi1
    
    if ~isempty(phi2(k))

        tmp.coeffs = [tmp.coeffs, phi2(k).coeffs];
        tmp.ind    = [tmp.ind, phi2(k).ind];
    end
    
    
    tmp = parametrisation_collapse(tmp); % add coefficients of same variate monomials
    tmp = flow_2_revlex(tmp); % In reverse lexicographic ordering
    
    phi12(k)=tmp;
    
    
end
phi12 = flow_2_revlex(phi12);


end