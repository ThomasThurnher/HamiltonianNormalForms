function [symExpr] = remove_zeros(symExprIn, vars,tol, varargin)

% {
if isempty(varargin)
    iter = 0;
else 
    iter = varargin{1};
end
% Remove small terms
childs = children(symExprIn);
childsnum = cellfun( @(a) double(subs(a,vars,ones(1,numel(vars)))), childs,'UniformOutput',false);
childsnum = cell2mat(childsnum);

I_imag = abs(imag(childsnum)) > tol;
I_real = abs(real(childsnum)) > tol;

I_out = I_imag | I_real;



symExpr = childs(I_out);


symExpr = sum(cell2sym(symExpr));

% If multiple children with same algebraic expression were present
%if iter < 2
%    [symExpr] = remove_zeros(symExpr, vars,tol, iter+1);
%end
%}

end