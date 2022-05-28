function [phi] = remove_slaves(phi,slaves)

for i = 1:numel(phi)
     phi(i) = remove_slaves_i(phi(i),slaves);
    
end
end

function [phii] = remove_slaves_i(phi,slaves)

% find columns which contain nonzero entries in slave rows
[row,col] = find(phi.ind);

slaveidx = (row == slaves(1) | row == slaves(2));

slavecols = unique(col(slaveidx));

numind = size(phi.ind,2);


mastercols = ~(repmat(1:numind,[numel(slavecols),1]) == repmat(slavecols,[1,numind]));
mastercols = find(sum(mastercols,1)== numel(slavecols)) ;

phii.ind = phi.ind(:,mastercols);
phii.coeffs = phi.coeffs(:,mastercols);
end