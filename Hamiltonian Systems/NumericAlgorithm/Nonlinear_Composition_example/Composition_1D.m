clear all
% {
%%
J = complex_structure(1);
%% Setup example Hamiltonian in original coordinates
% $$\mathcal{H} = \frac{\omega_1}{2} (q_1^2 +p_1^2 ) + \gamma q_1^3$$

omega1 = 1;
gamma  = 1;
order = 3;
Hpq(2).coeffs = [omega1/2 omega1/2];
Hpq(2).ind    =  2* eye(2);

Hpq(order).coeffs = [ gamma ];
Hpq(order).ind    = [ order   ;...
                  0   ];
              
Hpq(order+1).coeffs = [];
Hpq(order+1).ind = [];
JA = [ 0    omega1  ;...
      -omega1  0 ];
%% Move to complex coordinates, diagonalising the Hamiltonian
% With linear trafo $\mathcal{C}(p,q) = (x,y)$
%% Transform the Hamiltonian
% Hamiltonian $H(x_1,x_2,y_1,y_2) = -i\omega_1 x_1 y_1 +   \gamma  (V^{-1} q_1)^3$

[Hxy,Vh] = transform_Hamiltonian(JA,omega1,gamma,order);

n = size(Hxy(2).ind,1)/2;
%% Cubic Normalisation
% Generating Function of 3rd order Normalisation compute nonlinear coordinate 
% change for cubic Hamiltonian
% 
% Find $P3$ such that $H_3 + \{ N_2 , P\} = 0$

[P] = compute_generator( Hxy(2).coeffs, Hxy(order) , 1e-5);
%% 
% Compute $H \circ \phi$ where $\phi = e^{X_{P}}$ is the time 1 flow induced 
% by polynomial $P$ by means of Taylor expansion as 
% 
% $$H \circ \phi = H + \{ H, P \} + \frac{1}{2!} \{ \{ H, P \} , P \} + ...$$

% Normalise full Hamiltonian using the generating Polynomial
Hxy_3  = normalise_H(Hxy,order,P);
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

m = order;
phitilde3 = hamiltonian_flow(XP3,m,4,n);
%% Quartic Normalisation
% Generating Function of 3rd order Normalisation compute nonlinear coordinate 
% change for cubic Hamiltonian
% 
% Find $P3$ such that $H_3 + \{ N_2 , P\} = 0$

[P4] = compute_generator( Hxy_3(2).coeffs, Hxy_3(4) , 1e-5);
%% 
% Compute $H \circ \phi$ where $\phi = e^{X_{P}}$ is the time 1 flow induced 
% by polynomial $P$ by means of Taylor expansion as 
% 
% $$H \circ \phi = H + \{ H, P \} + \frac{1}{2!} \{ \{ H, P \} , P \} + ...$$

% Normalise full Hamiltonian using the generating Polynomial
Hxy_4  = normalise_H(Hxy_3,4,P4);
%% 
% Compute the vector field induced by the polynomial $P$ given as  $X_P = J  
% \nabla P$

% Get transformation
XP4 = polynomial_vectorfield(P4);
%% 
% Compute the flow induced by $X_P$ for polynomial of order $m$ on a $2n$ dimensional 
% phasespace. The coefficients are calculated up to polynomials of order $order$. 
% 
% $$\phi = 1 + X_P + \frac{1}{2!} X_P^2 + ...$$

m = 4;
phitilde4 = hamiltonian_flow(XP4,m,5,n);
%% 
% The linear transformation to complex normal coordinates and the nonlinear 
% flow then have to be composed as the normalised Hamiltonian is given as 
% 
% $$K = H \circ \mathcal{C} \circ \phi$$

N=2;
[phitilde] = compose_flow(phitilde3,phitilde4,N,4);
[phi] = compose_linear_flow(Vh,phitilde);
phi = flow_2_lex(phi);
phitilde = flow_2_lex(phitilde);
%% Get SSM Transformation
% $$\dot{q}_1 = \frac{\partial \mathcal{H}}{\partial p_1} = \omega_1p_1$$
% 
% $$\dot{p}_1 = -\frac{\partial \mathcal{H}}{\partial q_1} = - \omega_1q_1 - 
% 3\gamma q_1^2$$

[A,B,F] = build_model_H(Hpq,2);
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
masterModes = [1,2]; 
S.choose_E(masterModes);
[W,~] = S.compute_whisker(order);
%% Get SSM Transformation
% In complex coordinates $x,y$

[Ax,Bx,Fx] = build_model_H(Hxy_3,2);
%% Dynamical system setup

DSx = DynamicalSystem();
DSx.CanonicalTrafo = true;
set(DSx,'A',Ax,'B',Bx,'F',Fx);
set(DSx.Options,'Emax',5,'Nmax',10,'notation','multiindex')
set(DSx,'order',1);
[Vx,Dx,~] = DSx.linear_spectral_analysis();
%% 
% *Choose Master subspace (perform resonance analysis)*

Sx = SSM(DSx);
set(Sx.Options, 'reltol', 0.1,'notation','multiindex')
masterModes = [1,2]; 
Sx.choose_E(masterModes);
[Wtilde,~] = Sx.compute_whisker(order);
%%

% {
for i = order-1:order
' ############# Order #############'

i
'SSM Coeffs and indices'

'W'
full(W(i).coeffs)
'Wtilde'
full(Wtilde(i).coeffs)
full(Wtilde(i).ind.')
full(phitilde4(i).coeffs)
'HNF Coeffs and indices'
'Phi'
full(phi(i).coeffs)
'Phitilde'
full(phitilde(i).coeffs)
full(phitilde(i).ind)

'W'
%full(W(i).coeffs)-full(phi(i).coeffs)
'Wtilde'
%full(Wtilde(i).coeffs)-full(phitilde(i).coeffs)



end
%%
% Check changing property of SSM flow
[W_from_tilde] = compose_linear_flow(V,Wtilde);

'Compare'
'W'
full(W(order).coeffs)
'Wtilde'
full(W_from_tilde(order).coeffs)

% Check transformation property of HNF flow
% Trivially true as phi itself is obtained via the transformation