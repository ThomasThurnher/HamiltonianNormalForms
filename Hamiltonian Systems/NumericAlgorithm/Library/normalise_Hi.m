function [Hi] = normalise_Hi( H, i , m ,P)
% Transforms order i >= m Hamiltonian for normalisation at order m
% P - generating function of normalising transformation

tol = 1e-4;

Hi(1) = H(i);

J = (m-2):(m-2):(i-2);

jj = 2;
for j = J
    
    jm = j/(m-2);
    % Compute n- fold poisson bracket
    PBj = Poisson_Bracket_j(jm, H(i-j),P);

    % Contributions to the order i Hamiltonian 
    Hi(jj).coeffs = 1/factorial(jm)   * [PBj.coeffs];
    Hi(jj).ind    = [PBj.ind];
    
    jj = jj +1;
end

tmp.coeffs = [Hi.coeffs];
tmp.ind    = [Hi.ind];

Hi = polynomial_collapse(tmp,tol);

end