clc,clear all
%% Initialise Hamiltonian
omega1 = 1;
omega2 = 1;
H(2).coeffs = [omega1, omega2];
H(2).ind    = [ 1    0      ; ...
                   0    1       ; ...
                   1    0       ;...
                   0    1       ];
n = 2;

a = 1;
b = 1;

c = 0;
d = 0;
e = 0;
f = 0;
g = 0;
h = 0* pi/2;
k = 0*pi;

H(3).coeffs = [ a b c d e f g h k];
H(3).ind    = [ 2 1 1 1 0 0 1 0 0;...
                0 2 1 0 2 2 2 3 0
                0 0 1 1 0 1 0 0 0
                1 0 0 1 1 0 0 0 3];
               
order = 4;
reduction = false;
Pi = compute_generator( H(2).coeffs, H(3) , reduction, 1e-5);
%% Compute the poisson brackets

phi = hamiltonian_flow_2(Pi,3,order,2);
phi(4).coeffs
%slaves = [2,4];
%[phi] = remove_slaves(phi,slaves);

%phi(3).coeffs
%phi(4).coeffs


Hs   = cell(order,1);
phis = cell(order,1);

Hs{2} = H;


% Normalise full Hamiltonian using the generating Polynomial
H  = normalise_H(H,3,Pi);

% Get transformation
XPi = polynomial_vectorfield(Pi);

phi_i = hamiltonian_flow(XPi,3,order,n,reduction);

phi_i(4).coeffs
