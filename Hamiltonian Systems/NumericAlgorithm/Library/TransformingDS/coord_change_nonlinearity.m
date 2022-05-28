function [Fx] = coord_change_nl(F,V)
F = num2cell(F);

for i = 1:numel(F)
    if ~isempty(F{i})
        F{i}.ind = F{i}.ind.';
    end
end

    % compute Fi \circ V   
    F = compose_linear_flow_ssm(V,F,'right');
    % compute   V^-1 \circ Fi
    F = compose_linear_flow_ssm(inv(V),F,'left');

for i = 1:numel(F)
    if ~isempty(F{i})
        F{i}.ind = F{i}.ind.';
    end
end

for i = 1:numel(F)
    if ~isempty(F{i})
        
        Fx(i).coeffs = F{i}.coeffs;
        Fx(i).ind = F{i}.ind;
    end
end
end