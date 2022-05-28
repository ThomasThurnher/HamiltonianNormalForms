function [Ft] = transform_nl(Fm)

for i = 1:numel(Fm)
    if ~isempty(Fm(i)) && ~isempty(Fm(i).coeffs)
    Ft{i} = multi_index_to_tensor(Fm(i).coeffs,Fm(i).ind);
    Ft{i} = sptensor(Ft{i});
    else
        sizei = 2*ones(1,i+1);
        Ft{i} = sptensor(sizei);
    end 
end
end