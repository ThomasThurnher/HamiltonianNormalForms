clear all
%% Quadratic Spring Shaw Pierre Example
% Hamiltonian mass spring system with a cubic mass
% 
% $$\mathcal{H} =  \frac{1}{2} (m_1\dot{x}_1^2 + m_2\dot{x}_2^2 ) + \frac{1}{2}k_1( 
% x_1 - x_2)^2 + \frac{1}{2} k_0 (x_1^2 + x_2^2) + g x_1^3$$
% 
% $$\dot{x}_1 = \frac{\partial \mathcal{H}}{\partial p_1} = \dot{x}_1$$
% 
% $$\dot{x}_2 = \frac{\partial \mathcal{H}}{\partial p_2} = \dot{x}_2$$
% 
% $$\dot{p}_1 = -\frac{\partial \mathcal{H}}{\partial x_1} = - x_1 (k_0 + k_1) 
% + k_1 x_2 - 3 g x_1^2$$
% 
% $$\dot{p}_2 = -\frac{\partial \mathcal{H}}{\partial x_2} = - x_2 (k_0 + k_1) 
% + k_1 x_2$$
% 
% 
% 
% In complex coordinates $(x, \dot{x}) = V (\xi , \eta)$  the Hamiltonian reads
% 
% $$\mathcal{H} = \omega_1 \xi_1 \eta_1 +  \omega_2 \xi_2 \eta_2 +g x_1^3(\xi_1, 
% \xi_2, \eta_1, \eta_2)$$ 

param_order = 4;
P = build_model(); % Polynomial representation of Hamiltonian

H = Hamiltonian();
set(H,'param_order',param_order); % needs to be provided before setting H
set(H,'H',P,'n',2); % set hamiltonian
set(H,'reduction',true);

[A,B,F] = H.dynamicalSystem();
%% Transform Hamiltonian to complex diagonal coordinates

[Hxy,Vh] = H.transform_Hamiltonian();
%% Normalisation

resModes = [1,3];
[phitilde,phi] = H.HamiltonianNormalForm(param_order,resModes);
%% Model for SSM Computation
% 
% Dynamical system setup

DS = DynamicalSystem();
DS.CanonicalTrafo = false;
set(DS,'A',A,'B',B,'F',F);
set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
set(DS,'order',1);

[V,D,~] = DS.linear_spectral_analysis();
%% 
% *Choose Master subspace (perform resonance analysis)*

S = SSM(DS);
set(S.Options, 'reltol', 0.1,'notation','multiindex')
masterModes = [3,4]; 
S.choose_E(masterModes);
[W,R] = S.compute_whisker(param_order);
%% Compare Backbone Curves
% SSM Theory

omegaRange = [0.8,1.2];
order      = param_order;

set(S.FRCOptions,'outdof',1);
%S.extract_backbone(masterModes, omegaRange, order);
figure(),figure()
figHandles = findobj('Type','figure');
figs = [figHandles(2),figHandles(1)];
% HNF Theory
% 

set(H.BBOptions,'outdof',1);

H.extract_backbone(1, omegaRange, order, phi,figs);