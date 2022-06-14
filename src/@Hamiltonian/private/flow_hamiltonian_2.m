function [phi] = flow_hamiltonian_2(P,m,order,n,varargin)
% This function computes the time 1 flow induced by the hamiltonian vector field
% X_P given as phi = exp(X_P)
% P   - Hamiltonian Function
% m     - Order of the homogenous polynomial P
% order - order of polynomials up to which the flow should be computed

% Flow Phi consists of nonlinear functions
phi = repmat(polynomial_initialisation(),order,1);
N = 2*n;

if isempty(varargin)
    reduction = false;
    
    % Near identity character of transformation

else
    reduction = varargin{1};
    resModes = varargin{2};
    
    %tmp.coeffs = double( [1:N].' == resModes);
    %tmp.ind    = eye(numel(resModes));
end

tmp.coeffs = eye(2*n);
tmp.ind    = eye(2*n);
phi(1) = tmp;


PBn = coordinate_polynomials(2*n); % Nfold Poisson Bracket 
PBn = coordinate_PB(PBn,P,2*n);
VF_order = m-1;
for i = 2:order

    tmp = vectorfield2parametrisation(PBn,n);
    
    W.coeffs = 1/factorial(i-1) * tmp.coeffs;
    W.ind    = tmp.ind;
    
    if reduction
    %W = reduce_W(W,VF_order,resModes);
    end
    
    phi(VF_order) = W;
    
    
    VF_order = 1 + i*(m-2);
    if VF_order > order
        for j = 1:numel(phi)
            [phi(j)] = parametrisation_collapse(phi(j));
        end
        %'output lex flow in hamiltonian_flow_2'
        phi = flow_2_revlex(phi);
        
        % Expansion has been obtained to desired order
        return
    else
        
    % Construct
    PBn = coordinate_PB(PBn,P,2*n);        
    end

end

phi = flow_2_revlex(phi);
[phi] = parametrisation_collapse(phi);
end


function [PBn] = coordinate_PB(PBn,P,N)

for i = 1:N
    PBn(i) = polynomial_Poisson_Bracket(PBn(i),P);
end
end

function [PBn] = coordinate_polynomials(N)

PBn = repmat(polynomial_initialisation(),N,1);
for i = 1:N
    PBn(i).coeffs = 1;
    PBn(i).ind    = (1:N == i).';
end
end

function [Wout] = reduce_W(Win,i,resModes)

if ~isempty(Win.ind)
ind = Win.ind(resModes,:);

%find all zero entries
idx = find( sum(ind,1) == i);

Wout.coeffs = Win.coeffs(:,idx);
Wout.ind = ind(:,idx);
else
    Wout = polynomial_initialisation();
end
end