function [P] = build_model()
%% Makes dynamical system model of the quadratic Shaw Pierre

k1 = 1;
k2 = 1;
m1 = 1;
m2 = 1;
gamma = pi; % Cubic coefficient


% Hamiltonian
H2.coeffs = [m1/2, m2/2, 1/2*(k1+k2), 1/2*(k1+k2), -2*k2/2];
H2.ind    = [ 0    0     2           0          1     ; ...
                0    0     0           2          1     ; ...
                2    0     0           0          0     ;...
                0    2     0           0          0    ];

H3.coeffs = gamma;
H3.ind    = [3,0,0,0].';

P(2) = H2;
P(3) = H3;

end