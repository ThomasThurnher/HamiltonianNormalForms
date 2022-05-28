function [K,Vh] = transform_Hamiltonian(H,N,varargin)

if isempty(varargin)
    tol = 1e-10;
else
    tol = varargin{1};
end

%% Move to complex coordinates, diagonalising the Hamiltonian
% With linear trafo $\mathcal{C}(p,q) = (x,y)$
[JA,~,~] = build_model_H(H,N);
[Vh, ~, ~] = linear_transformation(JA);
%% Transform Hamiltonian with a linear coordinate change
% K = H \circ V


K = repmat(polynomial_initialisation(),numel(H),1);
N = size(Vh,1);
H_comp{1} = Vh;


for m = 2:numel(H)
    M = flip(sortrows(nsumk(N,m,'nonnegative')).',2);
    phi_k = H(m);
    
    H_k  = linear_composition_coefficients(Vh,H_comp,m); % Dependent on ordering, chooses conj, or revlex computation
    H_comp{m} = H_k;
    
    if ~isempty(phi_k.coeffs)
        
        tmp.coeffs = phi_k.coeffs* compute_pi_H(phi_k.ind,M,H_comp,N);
        tmp.ind    = M;
        tmp        = polynomial_collapse(tmp,tol);
        
        K(m)=tmp;
    end
end
end