clear all,clc
% {

%% Questions
% why does manual construction C and linear_transformation V not give same
% transformation x_can ?

% Why does construction C give negative eigenvalues in first coordinate
% directions in the complex normal coordinates ? 

% Does SSM computation work properly with the canonical linear trafo that
% is chosen?

J = complex_structure(2);
%% Setup example Hamiltonian in original coordinates
% $$\mathcal{H} = \frac{\omega_1}{2} (q_1^2 +p_1^2 ) + \frac{\omega_2}{2}(q_2^2 
% + p_2^2) + \gamma q_1^3$$

omega1 = 1;
omega2 = 8;
gamma  = 1;
Hpq(2).coeffs = [omega1/2 omega2/2 omega1/2 omega2/2];
Hpq(2).ind    =  2* eye(4);

Hpq(3).coeffs = [ gamma ];
Hpq(3).ind    = [ 3   ;...
                  0   ;
                  0;0];
              

JA = [ 0   0 omega1 0 ;...
      0   0 0 omega2 ;...
      -omega1  0 0 0 ;...
      0  -omega2 0 0 ];
%% Move to complex coordinates, diagonalising the Hamiltonian
% With linear trafo $\mathcal{C}(p,q) = (x,y)$

% Permutation matrix
P = [1,0,0,0;...
     0,0,1,0;...
     0,1,0,0;...
     0,0,0,1];

Tj =      1/sqrt(2)* [ 1/sqrt(1i) ,sqrt(1i) ;  -1/sqrt(1i) ,sqrt(1i)];
Tjinv = sqrt(1/(2i)) * [1i , -1i ; 1 , 1];
%Tj*Tjinv

Tinv = blkdiag(Tjinv,Tjinv);
T = blkdiag(Tj,Tj);


Ci    = P.' * Tinv * P;
C     = P.' * T * P;
Compl = C * JA * Ci;

[V, D, ~] = linear_transformation(JA); % JA V = V D - Why different from
% manual construction?


%% Transform the Hamiltonian
% Hamiltonian $H(x_1,x_2,y_1,y_2) = \omega_1 x_1 y_1 + \omega_2 x_2 y_2 +  \gamma 
% (P T P^T q_1)^3$

%% Transform quadratic Hamiltonian
N2.coeffs =  [D(1),D(2)];
N2.ind    =        [  1   ,   0  ;...
                      0   ,   1  ;...
                      1   ,   0  ;...
                      0   ,   1  ];
Hxy(2) = N2;                
% Transform cubic Hamiltonian to complex coordinates

x_can = inv(V) * [1;0;0;0];
%x_2 .* x_2
%x_can = C* [1;0;0;0]

tmp1.coeffs = x_can.';
tmp1.ind    = eye(4);

tmp1 = polynomial_collapse(tmp1);

tmp2 = polynomial_multiplication(tmp1,tmp1);

tmp3 = polynomial_multiplication(tmp2,tmp1);


Hxy(3).ind  = tmp3.ind;
Hxy(3).coeffs = gamma * tmp3.coeffs;

n = size(Hxy(2).ind,1)/2;
%% Cubic Normalisation
% Generating Function of 3rd order Normalisation compute nonlinear coordinate 
% change for cubic Hamiltonian
% 
% Find $P3$ such that $H_3 + \{ N_2 , P\} = 0$

[P] = compute_generator( Hxy(2).coeffs, Hxy(3) , 1e-5);
%% 
% Compute $H \circ \phi$ where $\phi = e^{X_{P}}$ is the time 1 flow induced 
% by polynomial $P$ by means of Taylor expansion as 
% 
% $$H \circ \phi = H + \{ H, P \} + \frac{1}{2!} \{ \{ H, P \} , P \} + ...$$

% Normalise full Hamiltonian using the generating Polynomial
Hxy_n  = normalise_H(Hxy,3,P);
%% 
% Compute the vector field induced by the polynomial $P$ given as  $X_P = J 
% \nabla P$

% Get transformation
XP3 = polynomial_vectorfield(P);
%% 
% Compute the flow induced by $X_P$ for polynomial of order $m$ on a $2n$ dimensional 
% phasespace. The coefficients are calculated up to polynomials of order $order$. 
% 
% $$\phi = 1 + X_P + \frac{1}{2!} X_P^2 + ...$$

m = 3;
order = 2;
phi3 = hamiltonian_flow(XP3,3,2,n);
%% 
% The linear transformation to complex normal coordinates and the nonlinear 
% flow then have to be composed as the normalised Hamiltonian is given as 
% 
% $$K = H \circ \phi \circ \mathcal{C}$$
[phi_A] = compose_linear_flow(V,phi3,'right');
%% Get SSM Transformation
% $$\dot{q}_1 = \frac{\partial \mathcal{H}}{\partial p_1} = \omega_1p_1$$
% 
% $$\dot{q}_2 = \frac{\partial \mathcal{H}}{\partial p_2} = \omega_2p_2$$
% 
% $$\dot{p}_1 = -\frac{\partial \mathcal{H}}{\partial q_1} = - \omega_1q_1 - 
% 3\gamma q_1^2$$
% 
% $$\dot{p}_2 = -\frac{\partial \mathcal{H}}{\partial q_2} = - \omega_2q_2$$

[A,B,F] = build_model_real(Hxy,4);
%% Dynamical system setup

DS = DynamicalSystem();
DS.CanonicalTrafo = true;
set(DS,'A',A,'B',B,'F',F);
set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
set(DS,'order',1);

[V,D,~] = DS.linear_spectral_analysis();
%% 
% *Choose Master subspace (perform resonance analysis)*


S = SSM(DS);
set(S.Options, 'reltol', 0.1,'notation','multiindex')
masterModes = [1,2,3,4]; 
S.choose_E(masterModes);
[W,R] = S.compute_whisker(2);

% {
for i = 1:2
'Order '
i
'SSM Coeffs and indices'

full(W(i).coeffs)
full(W(i).ind).'


'HNF Coeffs and indices'
full(phi3{i}.coeffs)
full(phi3{i}.ind)

full(phi_A{i}.coeffs)
full(phi_A{i}.ind)


end
%}