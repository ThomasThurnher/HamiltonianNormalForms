clear all, clc
%% Setup example Hamiltonian in normal quadratic coordinates
ifact = 2.1;
Hxy(2).coeffs = [1 ifact*1];
Hxy(2).ind    = [1 0 ; 0 1 ; 1 0 ;0 1];

Hxy(3).coeffs = [1/2/sqrt(2), -1];
Hxy(3).ind    = [ 2,1,0,0; 0,3,0,0].';

n = size(Hxy(2).ind,1)/2;
%% Cubic Normalisation
% Generating Function of 3rd order Normalisation
% Compute nonlinear coordinate change for cubic Hamiltonian
[P3] = compute_generator( Hxy(2).coeffs, Hxy(3) , 1e-5);
Hxy(3)

P3.ind
P3.coeffs
% Normalise full Hamiltonian using the generating Polynomial
Hxy  = normalise_H(Hxy,3,P3);

% Get transformation
XP3 = polynomial_vectorfield(P3);
phi3 = hamiltonian_flow(XP3,3,3,n);

%{
%% Quartic Normalisation

% Generating Function of 3rd order Normalisation
% Compute nonlinear coordinate change for cubic Hamiltonian
[P4] = compute_generator( Hxy(2).coeffs, Hxy(4) , 1e-5);

% Normalise full Hamiltonian using the generating Polynomial
Hxy  = normalise_H(Hxy,4,P4);


% Get transformation
XP4 = polynomial_vectorfield(P4);
phi4 = hamiltonian_flow(XP4,4,4,n);

%% Compose transformations

%}

%% Get SSM Transformation

[A,B,F] = build_model_complex(ifact);

%% Dynamical system setup 

DS = DynamicalSystem();
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
[W,R] = S.compute_whisker(3);

'Order 2'
'SSM Coeffs and indices'
full(W(2).coeffs)
full(W(2).ind).'

'HNF Coeffs and indices'
full(phi3{2}.coeffs)
full(phi3{2}.ind)

%{
'Order 3'
'SSM Coeffs and indices'
full(W(3).coeffs)
full(W(3).ind).'

'HNF Coeffs and indices'
full(phi{3}.coeffs)
full(phi{3}.ind)
%}