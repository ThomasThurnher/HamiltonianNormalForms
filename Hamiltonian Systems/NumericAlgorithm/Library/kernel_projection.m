function [Pm_out] = kernel_projection(Pm,m,phi,N)
% Choose the normalising polynomial such that the parametrisation has zero
% coefficients on the kernel elements
% Coefficients in phi are assumed to be stored in revlex ordering


restol = 1e-5;

% Check if lower orders contribute to parametrisation
if isempty(phi(m-1)) || isempty(phi(m-1).coeffs)
    Pm_out = Pm;
    return
else
    phi_m = phi(m-1)
end

% Current implementation assumes 1DOF system
if N > 2
    error('Kernel projection only for 1DOF Systems')
end

if ~isempty(Pm.coeffs)    
    ind = Pm.ind;
    coeffs = Pm.coeffs;
else 
    ind = [];
    coeffs = [];
end


%% Determine whether there are nonzero Kernel Elements in the parametrisation
K = flip(sortrows(nsumk(N,m-1,'nonnegative')).',2); % Multi-indices of parametrisation


% Find position of internal resonances
[ ~, I1] = find(K(1,:)-K(2,:) == 1);

if 0 ~= phi_m.coeffs(1,I1)
    phi_m.coeffs(1,I1)
    phi_m.coeffs(1,I1)/ceil(m-1/2)
    coeffs =  [coeffs, - phi_m.coeffs(1,I1)/ceil((m-1)/2)]
    ind = [ind, [ceil((m-1)/2); ceil((m-1)/2)]]
    
end


Pm_out.coeffs = coeffs;
Pm_out.ind    = ind;