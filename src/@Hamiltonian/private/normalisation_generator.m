function [P] = normalisation_generator( Lambda, Hm , resModes, reduction,restol)
% This function computes the generating polynomial of the canonical
% transformation normalising the ith term of a Hamiltonian which is in
% normal form up to order m-1, including the quadratic part of it
%
% Lambda - spectrum of linearised dynamical system, ordered like canoncial
%          coordinates
% res    - contains all resonances of linearised spectrum
% m      - order of Hamiltonian which is being normalised

tmp = polynomial_initialisation();


H_coeffs  = Hm.coeffs;
H_ind     = Hm.ind;
n         = size(H_ind,1)/2; % number of coordinates
N = 2*n;
if ~isempty(H_ind)
%% Find resonant multi-indices

if size(Lambda,1) > 1
    Lambda = Lambda.';
end
% Lambda contains linear coefficients of H in its columns
indicator     =  Lambda *( H_ind(n+1:end,:) - H_ind(1:n,:));

res_idx   =  abs(indicator) < restol;
%% Compute generating function
ii = 1;
if reduction
    %{
    %% Find reduction indicator
    target = [1;3];
    byproduct = [2;4];
    
    % Find multi-indices where byproduct indices at most contain 1 contribution
    byproductsum = sum(H_ind(byproduct,:),1);
    [~,byproductidx ] = find (byproductsum == 1);
    [~,targetidx] = find(byproductsum == 0);
    
    idcs = [byproductidx,targetidx];
    %}
    target = resModes;
    byproduct = find( ~ sum([1:N].' == target, 2)).';
    
    % Find multi-indices where byproduct indices at most contain 1 contribution
    targetsum = sum(H_ind(target,:),1);
    [~,targetidx ] = find (targetsum > 0);
    %[~,targetidx] = find(byproductsum == 0);
    
    idcs = targetidx;
else
    idcs = 1:numel(res_idx);
end

res_idx = res_idx(idcs);
for h = H_coeffs(idcs)
    if ~res_idx(ii)
        % nonResonant terms
        tmp(ii).ind = H_ind(:,idcs(ii));
        tmp(ii).coeffs = - h / indicator(idcs(ii));

    end
    ii = ii+1;
end
end

% Merge resulting monomials
P.ind = [tmp.ind];
P.coeffs = [tmp.coeffs];
end
