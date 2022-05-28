function [PBj] = Poisson_Bracket_j(j,H,P)
% Computes the n-fold poisson bracket ad_P^j(H)

if j > 1
    PBj = Poisson_Bracket(Poisson_Bracket_j(j-1,H,P), P);
elseif j == 1
    PBj = Poisson_Bracket(H,P);
else
    error('Cannot compute j fold PoissonBracket for nonpositive j')
end
end