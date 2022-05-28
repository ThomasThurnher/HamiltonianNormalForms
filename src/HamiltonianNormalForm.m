function [phi,phis,Hs,phicell] = HamiltonianNormalForm(H, order,n,varargin)
% Input H in quadratic complex normal form

if isempty(varargin)
    reduction = false;
    resModes = [];
else 
    reduction = varargin{1};
    if reduction
        fprintf('Performing reduction with modes \n')
        resModes = varargin{2}
    else
        resModes = [];
    end
end

Hs   = cell(order,1);
phis = cell(order,1);

Hs{2} = H;
phi = repmat(polynomial_initialisation(),order,1);
for i = 3:order+1
    fprintf('############## Normalisation order of Hamiltonian: %d ################ \n',i)

    
    Pi = compute_generator( H(2).coeffs, H(i), resModes , reduction, 1e-5);
    %Pi
    if ~isempty(Pi) && ~isempty(Pi.coeffs)
        % Normalise full Hamiltonian using the generating Polynomial
        H  = normalise_H(H,i,Pi);
        
        % This implmentation has a bug for H = H2 + x1^2y2 +x1x2^2
        %{
        % Get transformation
        XPi = polynomial_vectorfield(Pi);
        %phi_i_2 = hamiltonian_flow(XPi,i,order,n,reduction);
        %}
        if ~reduction
            phi_i = hamiltonian_flow_2(Pi,i,order,n);
        else
            phi_i = hamiltonian_flow_2(Pi,i,order,n,reduction, resModes);
        end
        %compare_flows(phi_i,phi_i_2,i)
        
        Hs{i} = H;
        phis{i-2} = phi_i;
        
        if i == 3
            phi = phi_i;
        else
            phi = compose_flow(phi,phi_i,2*n,order,reduction,resModes);
        end
        
        %' ----- Resulting flow -----'
        
        %phishow = flow_output(phi,order);
        %phishow.coeffs
    end
    phicell{i-2} = phi;
end

phi = flow_2_lex(phi);
phi = flow_sparsification(phi);

for i = 1:numel(phi)
    if ~isempty(phi(i))
    phi(i) = reduce_phi(phi(i),i,resModes);
    end
end

end


function [phi_out] = reduce_phi(phi_in,i,resModes)

if ~isempty(phi_in.ind)
ind = phi_in.ind(resModes,:);

%find all zero entries
idx = find( sum(ind,1) == i);

phi_out.coeffs = phi_in.coeffs(:,idx);
phi_out.ind = ind(:,idx);
else
    phi_out = polynomial_initialisation();
end
end
