%% Proof Theorem 1 for 2D Hamiltonian

% Quadratic Hamiltonian
omega1 = 1;
omega2 = 1;
H(2).coeffs = [omega1, omega2];
H(2).ind    = [ 1    0      ; ...
                   0    1       ; ...
                   1    0       ;...
                   0    1       ];
               
               
order = 3;
c_1 = 1;
c_2 = pi;

for alpha = 0:order -1
for beta = 0:order-alpha-1
theta = order - alpha - beta;

if alpha >= 0 && beta >=0 && theta >= 0
for delta = 0:order-1
for gamma = 0:order-delta-1
nu = order - delta - gamma;

if delta >= 0 && gamma >=0 && nu >= 0
if alpha+beta+delta+gamma == order

    parameters = [alpha, beta, theta, delta, gamma, nu].';
    names = {'alpha', 'beta', 'theta', 'delta', 'gamma', 'nu'}.';
    table(names,parameters)
    
    H(3).coeffs = [ c_1   c_2 ];
    H(3).ind    = [ alpha delta ;...
                   theta  0 
                   beta   gamma 
                   0      nu ];
    
     phi = hamiltonian_flow_2(H_in(3),3,5,2);


    slaves = [2,4];
    [phi] = remove_slaves(phi,slaves);
    phi(2).coeffs
    phi(3).coeffs
          
end
end

end
end
% delta and gamma loop finish
end
end
end
% alpha and beta loop finish