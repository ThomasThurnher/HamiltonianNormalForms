function [flow] =  flow_sparsification(phi) 

for i = 1:numel(phi)
    
    cols = find(~all(phi(i).coeffs==0));
    
    flow(i).coeffs = phi(i).coeffs(:,cols);
    flow(i).ind    = phi(i).ind(:,cols);
end

end