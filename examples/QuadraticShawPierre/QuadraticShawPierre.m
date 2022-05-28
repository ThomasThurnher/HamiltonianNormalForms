clear all
%% Quadratic Spring Shaw Pierre Example
k1 = 1;
k2 = 1;
m1 = 1;
m2 = 1;
gamma = pi; % Cubic coefficient
%% Initialise Hamiltonian

H_in(2).coeffs = [m1/2, m2/2, 1/2*(k1+k2), 1/2*(k1+k2), 2*k2/2];
H_in(2).ind    = [ 0    0     2           0          1     ; ...
                0    0     0           2          1     ; ...
                2    0     0           0          0     ;...
                0    2     0           0          0    ];

H_in(3).coeffs = gamma;
H_in(3).ind    = [3,0,0,0].';
            
            
J = [ 0  0 1 0 ;
      0  0 0 1 ;
     -1  0 0 0 ;
      0 -1 0 0 ];

A11 = [ k1 + k2 , -k2 ; -k2 , k1 + k2 ];
 
A = [A11 , zeros(2,2) ; zeros(2,2), eye(2)];
 
 
 %% Quadratic Normalisation
 [V,D] = eigs(J*A);
 Spectrum = diag(D);
 Lambda = imag(Spectrum([1,3]));
 omega1 = Lambda(1);
 omega2 = Lambda(2);
%% Real Normalisation
V1 = real(V(:,1)); V2 = imag(V(:,1)); V3 = real(V(:,3)); V4 = imag(V(:,3)); 

% Normalise with condition V_i^T * J * V_j = delta_ij
scale1 = V1.' * J * V2;
V1 = V1 / sqrt(scale1); V2 = V2 / sqrt(scale1);
scale2 = V3.' * J * V4;
V3 = V3 / sqrt(scale2); V4 = V4 / sqrt(scale2);

Trafo = [V1, V3, V2, V4];
%T.'*J*T
%% New Hamiltonian
% Quadratic real Hamiltonian
N2_real.coeffs = [ omega1 , omega1 , omega2 , omega2];
N2_real.ind    = [2         0          0       0    ;...
                 0         0          2       0    ; ...
                 0         2          0       0    ; ...
                 0         0          0       2    ];

          
H = polynomial_initialisation(); 

%% Move to complex coordinates, diagonalising the Hamiltonian

% Permutation matrix
P = [1,0,0,0;0,0,1,0;0,1,0,0;0,0,0,1];
Tj = 1/sqrt(2i)* [ 1i ,1 ; -1i, 1];
T = blkdiag(Tj,Tj);


Complexify = P * Trafo; %Permute to have pairs of coordinates in adjacent positions
Complexify = T * Complexify; % Move to complex coordinates
Complexify = P.' * Complexify; %Permute to regain canonical ordering of variables


%% Transform quadratic Hamiltonian
N2.coeffs = 1/1i * [omega1,omega2];
N2.ind    =        [  1   ,   0  ;...
                      0   ,   1  ;...
                      1   ,   0  ;...
                      0   ,   1  ];
H(2) = N2;
                
%% Transform cubic Hamiltonian to complex coordinates
% Transform x1 to normalised coordinates
x_can = Complexify* [1;0;0;0];

tmp1.coeffs = x_can.';
tmp1.ind    = eye(4);
tmp3 = polynomial_multiplication(polynomial_multiplication(tmp1,tmp1),tmp1);

H(3).ind  = tmp3.ind;
H(3).coeffs = gamma * tmp3.coeffs;

%H(3).ind = H(3).ind(:,1);
%H(3).coeffs = H(3).coeffs(1);
%% Normalisation of cubic Hamiltonian

%% Generating Function of 3rd order Normalisation
% Compute nonlinear coordinate change for cubic Hamiltonian
[P3] = compute_generator( N2.coeffs, H(3) , 1e-5);


%% Normalise full Hamiltonian using the generating Polynomial
H  = normalise_H(H,3,P3);

%% Transformation to real canonical coordinates
Tjinv = sqrt(1i)/sqrt(2)   * [-1i,1i;1,1];

Tinv = blkdiag(Tjinv,Tjinv);
Realise = P.' * Tinv * P;


%% Get transformation
XP3 = polynomial_vectorfield(P3);

n = size(H(2).ind,1)/2;
phi = hamiltonian_flow(XP3,3,3,n);


Transformation = compose_linear_flow(Complexify,phi,'right');
'complex'
%Transformation = compose_linear_flow(Realise,Transformation,'left');


%% Get SSM Transformation

[M,C,K,fnl] = build_model(m1,m2,k1,k2,gamma);

%% Dynamical system setup 

DS = DynamicalSystem();
set(DS,'M',M,'C',C,'K',K,'fnl',fnl);
set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
[V,D,~] = DS.linear_spectral_analysis();
%% 
% *Choose Master subspace (perform resonance analysis)*

S = SSM(DS);
set(S.Options, 'reltol', 0.1,'notation','multiindex')
% set(S.Options, 'reltol', 0.1,'notation','tensor')
masterModes = [1,2,3,4]; 
S.choose_E(masterModes);
[W,~] = S.compute_whisker(2);