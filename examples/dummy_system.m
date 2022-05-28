clear all;
% dummy hamiltonian in 1 DOf, complex normal form with cubic nonlinearity

%% Build Hamiltonian
H(2).coeffs = [1j,2j];
H(2).ind    = [1 0;  %q1
               0 1;  %q2
               1 0   %p1
               0 1]; %p2

H(3).coeffs = [1/2,-1  3];
H(3).ind    = [ 2  0   1 ;
                1  3   1 ; 
                0  0   1 ; 
                0  0   0 ];

H(3).coeffs = H(3).coeffs(1:2);
H(3).ind = H(3).ind(:,1:2);

H(4) = polynomial_initialisation();
%H(8) = polynomial_initialisation;

Lambda = H(2).coeffs;

%% Generating Function of 3rd order Normalisation
% Compute nonlinear coordinate change for cubic Hamiltonian
[P] = compute_generator( Lambda, H(3) , 1e-5);

% Compute normalised Hamiltonian
%N3 = normalise_Hi(H,3,3,P)
%H4 = normalise_Hi(H,4,3,P)

%% Normalise full Hamiltonian using the generating Polynomial
H  = normalise_H(H,3,P);