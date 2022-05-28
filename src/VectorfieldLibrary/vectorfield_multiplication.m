function [XY] = vectorfield_multiplication(X,Y,N)
% This function multiplies the vectorfields X and Y as XY
if isempty(X) || isempty(Y)
    XY = polynomial_initialisation();
    return
end



XY = repmat(polynomial_initialisation(),N,1);

for i = 1:N
    if ~isempty(X(i).coeffs)
        
        XYi = repmat(polynomial_initialisation(),1,N);
        for j = 1:N
            if ~isempty(Y(j).coeffs)
                
                
                dYidxj = polynomial_derivative(Y(i),j);
                
                XYi(j) = polynomial_multiplication(X(j),dYidxj);
            end
        end
        XY(i) = polynomial_collapse(XYi);
    end
end
end