function [phi_out] = flow_output(phi,order)
%%
% Storing the coefficients in suitable format for function output.
% They are stored in lexicographic ordering of the multi-indices
%%
phi_out = repmat(struct('coeffs',[],'ind',[]),1,order);


l = size(phi(1).ind,1);
for i = 1:order
    %create multi-indices
    if l >1
    K = sparse(sortrows(nsumk(l,i,'nonnegative')));
    else
        K=i;
    end
    % SSM coefficients with multi-indices
    idx_W_0     = all(phi(i).coeffs==0);
    phi_out(i).coeffs = phi(i).coeffs(:,~idx_W_0);
    phi_out(i).ind    =  K(~idx_W_0,:).';
    
end
end