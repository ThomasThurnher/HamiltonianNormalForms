function [phitilde,phi] = HamiltonianNormalForm(obj , order,varargin)
% Input H in quadratic complex normal form

reduction = obj.reduction;
n = obj.n;

if reduction
    if isempty(varargin)
        error('Choose reduction basis')
    end
    fprintf('Performing reduction with modes \n')
    resModes = varargin{1}
else
    resModes = [];
end


Hs   = cell(order,1);
phis = cell(order,1);

Hs{2} = obj.H;
phitilde = repmat(polynomial_initialisation(),order,1);
for i = 3:order+1
    fprintf('############## Normalisation order of Hamiltonian: %d ################ \n',i)

    
    Pi = normalisation_generator( obj.H(2).coeffs, obj.H(i), resModes , reduction, 1e-5);
    %Pi
    if ~isempty(Pi) && ~isempty(Pi.coeffs)
        % Normalise full Hamiltonian using the generating Polynomial
        obj.H  = normalise_H(obj.H,i,Pi);
        
        % This implmentation has a bug for H = H2 + x1^2y2 +x1x2^2
        %{
        % Get transformation
        XPi = polynomial_vectorfield(Pi);
        %phi_i_2 = hamiltonian_flow(XPi,i,order,n,reduction);
        %}
        if ~reduction
            phi_i = flow_hamiltonian_2(Pi,i,order,n);
        else
            phi_i = flow_hamiltonian_2(Pi,i,order,n,reduction, resModes);
        end
        %flow_compare(phi_i,phi_i_2,i)
        
        Hs{i} = obj.H;
        phis{i-2} = phi_i;

        if i == 3
            phitilde = phi_i;
        else
            phitilde = flow_compose(phitilde,phi_i,2*n,order,reduction,resModes);
        end
        
        %' ----- Resulting flow -----'
        
        %phishow = flow_output(phi,order);
        %phishow.coeffs
    end
    phicell{i-2} = phitilde;
end
phitilde = flow_2_lex(phitilde);
phitilde = flow_sparsification(phitilde);
if reduction
for i = 1:numel(phitilde)
    if ~isempty(phitilde(i))
    phitilde(i) = reduce_phi(phitilde(i),i,resModes);
    end
end
end
[phi] = compose_linear_flow(obj.V,phitilde);
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
