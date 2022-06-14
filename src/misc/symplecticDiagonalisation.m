function [V, D, W] = symplecticDiagonalisation(JA,varargin)

N = size(JA,1);
%{
if N == 4
[V, LAMBDA, W] = eig(JA);

[V, D, W] = sort_modes_canonical(V, LAMBDA, W);

J = complex_structure(2);

V1 = V(:,1); V2 = V(:,2); V3 = V(:,3);V4 = V(:,4);

% Normalise with condition R_1^T * J * S_1 = 1/4
scale1 = V1.' * J * V2;
V1 = V1 / sqrt(scale1); V2 = V2 / (sqrt(scale1));

scale2 = V3.' * J * V4;
V3 = V3 / sqrt(scale2); V4 = V4 / (sqrt(scale2));

%V = [V2, V4, V1, V3];
V = [V1, V3, V2, V4];
D = D([1,3,2,4]); % canonical ordering (x1, .. xn, y1, .. yn)
%% Check Symplecticity
if norm(V.' * J * V - J) > 1e-5
    warning('Linear Transformation is not symplectic')
end
W = W(:,[1,3,2,4]);
%% Normalise left eigenvectors accordingly
w_scale = diag((W'*V)).';
W = W./conj(w_scale);

elseif N == 5
   [V, LAMBDA, W] = eig(JA);
[V, D, W] = sort_modes_canonical(V, LAMBDA, W);

J = complex_structure(3);

V1 = V(:,1); V2 = V(:,2); V3 = V(:,3);V4 = V(:,4); V5 = V(:,5); V6 = V(:,6);

% Normalise with condition R_1^T * J * S_1 = 1/4
scale1 = V1.' * J * V2;
V1 = V1 / sqrt(scale1); V2 = V2 / (sqrt(scale1));

scale2 = V3.' * J * V4;
V3 = V3 / sqrt(scale2); V4 = V4 / (sqrt(scale2));

scale3 = V5.' * J * V6;
V5 = V5 / sqrt(scale3); V6 = V6 / (sqrt(scale3));

%V = [V2, V4, V1, V3];
V = [V1, V3, V5, V2, V4, V6];
D = D([1,3,5,2,4,6]); % canonical ordering (x1, .. xn, y1, .. yn)
%% Check Symplecticity
if norm(V.' * J * V - J) > 1e-5
    warning('Linear Transformation is not symplectic')
end
W = W(:,[1,3,5,2,4,6]);
%% Normalise left eigenvectors accordingly
w_scale = diag((W'*V)).';
W = W./conj(w_scale); 


elseif N == 6
    %}
[V, LAMBDA, W] = eig(JA);
[V, D, W] = sort_modes_canonical(V, LAMBDA, W);

J = complex_structure(N/2);

Vodd = V(:,1:2:end-1);
Veven = V(:,2:2:end);

scales = sum(Vodd .* (J * Veven),1);

Vodd = Vodd ./sqrt(scales);
Veven = Veven./sqrt(scales);


V = [Vodd, Veven];
maxidx = size(V,2);

D = D([1:2:maxidx-1,2:2:maxidx]); % canonical ordering (x1, .. xn, y1, .. yn)
%% Check Symplecticity
if norm(V.' * J * V - J) > 1e-5
    V.'*J*V
    warning('Linear Transformation is not symplectic')
end
W = W(:,[1:2:maxidx-1,2:2:maxidx]);
%% Normalise left eigenvectors accordingly
w_scale = diag((W'*V)).';
W = W./conj(w_scale); 
%{
elseif N == 2

    if ~isempty(varargin)
        V = varargin{1};
        D = inv(V)*JA*V
        W = eye(2)
        LAMBDA = diag(D);
    else
    [V, LAMBDA, W] = eig(JA);
    end
    [V, D, W] = sort_modes_canonical(V, LAMBDA, W);
    
    J = complex_structure(1);
    
    V1 = V(:,1); V2 = V(:,2);
    
    % Normalise with condition R_1^T * J * S_1 = 1/4
    scale1 = V1.' * J * V2;
    V1 = V1 / sqrt(scale1); V2 = V2 / (sqrt(scale1));

    
    %V = [V2, V4, V1, V3];
    V = [V1, V2];
    D = D([1,2]); % canonical ordering (x1, .. xn, y1, .. yn)
    %% Check Symplecticity
    if norm(V.' * J * V - J) > 1e-5
        warning('Linear Transformation is not symplectic')
    end
    W = W(:,[1,2]);
    %% Normalise left eigenvectors accordingly
    w_scale = diag((W'*V)).';
    W = W./conj(w_scale);
end#
%}
end

function [V, D, W] = sort_modes_canonical(V, Lambda, W)
%SORT_MODES: This function sorts the eigenvectors (V) in descending order of
%the real parts of the corresponding eigenvalues (D). The resulting
%eigenvectors are also normalized to unit magnitude. 

% obtain the eigenvalues as a vector instead of a diagonal matrix
D = diag(Lambda); 
if ~iscolumn(D)
    D = transpose(D);
end

% sort eigenvalues in canonical order for pairs of purely imaginary eigenvalue pairs in hamiltonian systems
[Lambda_sorted,I] = sortrows([abs(imag(D)) sign(imag(D))],[1 2],{'ascend' 'descend'});
D = 1i * Lambda_sorted(:,1) .*Lambda_sorted(:,2);
%D = diag(D);
% arrange eigenvectors accordingly
V = V(:,I);
W = W(:,I);
end

