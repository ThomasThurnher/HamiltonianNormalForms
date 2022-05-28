% Dummy Polynomial for operations
N = 2;
n = 1;
P.coeffs = [1];

P.ind   = [2 ;1 ];


XP = polynomial_vectorfield(P);


%XPXP = vectorfield_multiplication(XP,XP,N)

phi = hamiltonian_flow(XP,3,3,n);

A = [2,1;1,0];


A_phi = compose_linear_flow(A,phi);