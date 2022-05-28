function [J] = complex_structure(n)

J = [ zeros(n) , eye(n) ; ...
       - eye(n), zeros(n)];
   
end