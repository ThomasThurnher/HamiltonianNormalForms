function [V, D, W] = linear_spectral_analysis(obj)
if obj.CanonicalTrafo && isdiag(obj.A)
    % Do preserve canonical coordinate structure
    V = eye(obj.N);
    W = eye(obj.N);
    D = diag(obj.A);
    obj.spectrum.V = V;
    obj.spectrum.W = W;
    obj.spectrum.Lambda = D;
    return
end

if obj.N < obj.Options.Nmax
    if ~issparse(obj.A)
        if norm(obj.B-eye(size(obj.B)))<1e-8
            [V, LAMBDA, W] = eig(obj.A); % B=I, CONSISTENT WITH SSMTool 1.0
        else
            [V, LAMBDA, W] = eig(obj.A,obj.B);
        end
        
    else
        [V, LAMBDA, W] = eig(full(obj.A),full(obj.B));
    end
else
    
    E_max = obj.Options.Emax;
    n = obj.n;
    
    if obj.order == 2 && obj.Options.RayleighDamping
        
    disp(['Due to high-dimensionality, we compute only the first ' num2str(E_max) ' eigenvalues with the smallest magnitude. These would also be used to compute the spectral quotients'] )
    % Computing undamped eigenvalues and eigenvectors
    
    disp ('Assuming a proportional damping hypthesis with symmetric matrices') 
    [U, ~, NOT_CONVERGED] = eigs(sparse(obj.K),sparse(obj.M),E_max,'smallestabs');
    
    % Assuming proportional damping for estimating damped eigenvectors
    LAMBDA = zeros(2*E_max,2*E_max);
    V = zeros(2*n,2*E_max);
    
    % Rearranging eigenvalues
    for j = 1:E_max
        mu_j = U(:,j).'* obj.M * U(:,j);
        omega2_j = (U(:,j).'* obj.K * U(:,j))/mu_j;
        beta_j = (U(:,j).'* obj.C * U(:,j))/mu_j;
        
        fprintf('modal damping ratio for %d mode is %d\n', j, beta_j/(2*sqrt(omega2_j)));
        
        lambda1 = (-beta_j + sqrt(beta_j^2 - 4 * omega2_j) ) / 2;
        lambda2 = (-beta_j - sqrt(beta_j^2 - 4 * omega2_j) ) / 2;
        
                
        LAMBDA(2*j-1,2*j-1) = lambda1;
        LAMBDA(2*j,2*j) = lambda2;
        
        V(:,2*j-1) = [U(:,j); lambda1*U(:,j)];
        V(:,2*j) = [U(:,j); lambda2*U(:,j)];
    end
    W = conj(V);
    
    if NOT_CONVERGED
        error('The eigenvalue computation did not converge, please adjust the number of eigenvalues to be computed')
    end
    
    if ~issymmetric(obj.K) || ~issymmetric(obj.M)
        disp('the left eigenvectors may be incorrect in case of asymmetry of matrices')
%         [W, ~] = eigs(A',B',E_max,'smallestabs');
        % we assume here that both eigenvalue problems return eigenvalues in
        % the same order.
    end
    
    else
        % right eigenvectors
        [V, Dv] = eigs(obj.A,obj.B,E_max,'smallestabs');
        [Lambda_sorted,I] = sort(diag(Dv),'descend','ComparisonMethod','real');
        % further sort if real parts are equal (very close)
        [Lambda_sorted,II] = sort_close_real_different_imag(Lambda_sorted);
        LAMBDA = diag(Lambda_sorted);
        V = V(:,I(II));
        % rescale V w.r.t mass matrix
        if ~isempty(obj.M)
            for j=1:numel(I)
                vs = V(1:obj.n,j)'*obj.M*V(1:obj.n,j);
                V(:,j) = V(:,j)/sqrt(abs(vs));
            end
        end
        
        % left eigenvectors
        if issymmetric(obj.A) && issymmetric(obj.B)
            W = conj(V);
        else
            [W, Dw] = eigs(obj.A',obj.B',E_max,'smallestabs');
            [Lambda_sorted,I] = sort(diag(Dw),'descend','ComparisonMethod','real');
            % further sort if real parts are equal (very close)
            [Lambda_sorted,II] = sort_close_real_different_imag(Lambda_sorted);
            % make sure reordered Dw is consistent with reordered Dv
            assert(norm(Lambda_sorted-diag(LAMBDA))<1e-3*norm(diag(LAMBDA)),...
                'Orders for W and V are not consistent');
            W = W(:,I(II));
            W = conj(W);
        end
    end    
end
% Make Trafo canonical
[V,LAMBDA,W] = remove_stiff_modes(V,LAMBDA,W,obj.Options.lambdaThreshold);
[V,LAMBDA,W] = remove_zero_modes(V,LAMBDA,W);

if obj.CanonicalTrafo
    if obj.N == 4;
    [V, D, W] = sort_modes_canonical(V, LAMBDA, W);

    'Canonical Trafo'
    J = complex_structure(2);

    V1 = V(:,1); V2 = V(:,2); V3 = V(:,3);V4 = V(:,4); 

    % Normalise with condition R_1^T * J * S_1 = 1/4
    scale1 = V1.' * J * V2;
    V1 = V1 / sqrt(scale1); V2 = V2 / (sqrt(scale1));
    
    scale2 = V3.' * J * V4;
    V3 = V3 / sqrt(scale2); V4 = V4 / (sqrt(scale2));
    
    %V = [V2, V4, V1, V3];
    V = [V1, V3, V2, V4];
    D = D([1,3,2,4]);
    %% Check Symplecticity
    if norm(V.' * J * V - J) > 1e-5
        warning('Linear Transformation is not symplectic')
    end
    W = W(:,[1,3,2,4]);
    %% Normalise left eigenvectors accordingly
    w_scale = diag((W' * obj.B*V)).';    
    W = W./conj(w_scale);
    
    elseif obj.N == 6
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

    elseif obj.N == 2

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
end
else
    [V, D, W] = sort_modes(V, LAMBDA, W);
    [V,W] =  normalize_modes(V,W,obj.B);

end

% check the orthonoramlity of V and W with respect to B
numColV = size(V,2);
if norm(W'*obj.B*V-eye(numColV),'fro')/numColV>1e-4
    warning('V and W are not orthonormal');
end
obj.spectrum.V = V;
obj.spectrum.W = W;
obj.spectrum.Lambda = D;

fprintf('\n The first %d eigenvalues are given as \n',length(D))
disp(D)

end

function [V, Lambda, W] = sort_modes(V, D, W)
%SORT_MODES: This function sorts the eigenvectors (V) in descending order of
%the real parts of the corresponding eigenvalues (D). The resulting
%eigenvectors are also normalized to unit magnitude. 

% obtain the eigenvalues as a vector instead of a diagonal matrix
Lambda = diag(D); 
if ~iscolumn(Lambda)
    Lambda = transpose(Lambda);
end
% sort eigenvalues in the descending order of real parts, incase of tie by
% ascending order of magnitude of imaginary parts
[Lambda_sorted,I] = sortrows([real(Lambda), abs(imag(Lambda)) sign(imag(Lambda))],[1 2],{'descend' 'ascend'});
D = Lambda_sorted(:,1) + 1i * Lambda_sorted(:,2).*Lambda_sorted(:,3);
D = diag(D);
% arrange eigenvectors accordingly
V = V(:,I);
W = W(:,I);

% ensure positive imaginary part first in every complex pair
Lambda = diag(D);
skip = false;
for j = 1:length(Lambda)
    if skip 
        skip = false;
        continue;
    end    
    if ~isreal(Lambda(j))&& abs(Lambda(j)-conj(Lambda(j+1)))<1e-8*abs(Lambda(j))
        % extract complex eigenpair
        V0 = V(:,j:j+1);
        W0 = W(:,j:j+1);
        Lambda0 = Lambda(j:j+1);
        % sort eigenvalues in the descending order of imaginary parts
        [~,I] = sort(imag(Lambda0),'descend','ComparisonMethod','real');        
        % rearrange the ordre of the pair
        Lambda([j,j+1]) = Lambda0(I);
        V(:,[j,j+1]) = V0(:,I);
        W(:,[j,j+1]) = W0(:,I);
%         Lambda(j) = Lambda0(I(1));
%         V(:,j) = V0(:,I(1));
%         W(:,j) = W0(:,I(1));
        % ensure complex conjugate eigenvalues and eigenvectors - not true
        % if A and B are not real
%         Lambda(j+1) = conj(Lambda(j));                
%         V(:,j+1) = conj(V(:,j));
%         W(:,j+1) = conj(W(:,j));
        % move to the next pair of eigenvalues
        skip = true; 
    end
end
% D = diag(Lambda);
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


function [V,W] = normalize_modes(V,W,B)
%NORMALIZE_MODES: This function normalizes the right and left eigenvectors  
% W, V withrespect to the matrix B. 

% V = V*diag(1./vecnorm(V));
mu = diag(W'*B*V);

% V = V*diag(1./sqrt(mu));
% W = W*diag(1./(sqrt(mu)'));

W = W*diag(1./(mu'));
end

function [V,LAMBDA,W] = remove_stiff_modes(V,LAMBDA,W,lambdaThreshold)
%REMOVE_STIFF_MODES: This function removes modes with infinite eigenvalues
D = diag(LAMBDA);
D = abs(D);
n = numel(D);
idx0 = [find(isinf(D));find(isnan(D))]; % indices with inf and nan eigenvalues
D(idx0) = lambdaThreshold-1;
idx1 = find(D>lambdaThreshold);
idx2 = setdiff(1:n, [idx0;idx1]);
V = V(:,idx2);
W = W(:,idx2);
LAMBDA = LAMBDA(idx2,idx2);
if ~isempty(idx0)
    fprintf('%i nan/inf eigenvalues are removed\n',numel(idx0));
end
if ~isempty(idx1)
    fprintf('%i eigenvalues with mangnitude larger than %d are removed\n',...
        numel(idx1), lambdaThreshold);
end
end

function [V,LAMBDA,W] = remove_zero_modes(V,LAMBDA,W)
%REMOVE_STIFF_MODES: This function removes modes with zero eigenvalues
D = diag(LAMBDA);
D = abs(D);
n = numel(D);
idx1 = find(D<eps*max(D)); % indices with zero eigenvalues
idx2 = setdiff(1:n, idx1);
V = V(:,idx2);
W = W(:,idx2);
LAMBDA = LAMBDA(idx2,idx2);
if ~isempty(idx1)
    fprintf('%i zero eigenvalues are removed\n',numel(idx1));
end
end

function [y,idx] = sort_close_real_different_imag(x)
realx = real(x);
imagx = imag(x);
n = numel(realx);

drealx = realx(2:end)-realx(1:end-1);
drealx = drealx./(realx(1:end-1)+eps);
idgap  = find(abs(drealx)>1e-6);
idgap  = [1; idgap(:); n];
idx    = [];
for k=1:numel(idgap)-1
   ida = idgap(k);
   if k>1; ida=ida+1; end
   idb = idgap(k+1);
   tmp = ida:idb;
   [~,idab] = sort(imagx(tmp),'descend');
   idx = [idx tmp(idab)];
end
y = x(idx);

end