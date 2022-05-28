clear all,clc
%% Quadratic Spring Shaw Pierre Example

k1 = 1;
k2 = 1;
m1 = 1;
m2 = 1;
gamma = pi; % Cubic coefficient
param_order = 6;
N=4;
n = N/2;
%% Initialise Hamiltonian

Hpq(2).coeffs = [m1/2, m2/2, 1/2*(k1+k2), 1/2*(k1+k2), -2*k2/2];
Hpq(2).ind    = [ 0    0     2           0          1     ; ...
                0    0     0           2          1     ; ...
                2    0     0           0          0     ;...
                0    2     0           0          0    ];

Hpq(3).coeffs = [gamma, 2 * gamma];
Hpq(3).ind    = [[3,0,0,0].' , [0,2,0,1].' ];
if param_order >2
Hpq(param_order +1) = polynomial_initialisation();
end
%% Transform Hamiltonian to complex diagonal Coordinates

[Hxy,Vh] = transform_Hamiltonian(Hpq,N);
%% Normalisation

[phitilde,phis,Hs,phicell] = HamiltonianNormalForm(Hxy, param_order,n);
[phi] = compose_linear_flow(Vh,phitilde);
%% Model for SSM Computation
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


[A,B,F] = build_model_H(Hxy,N);
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
masterModes = [1,3]; 
S.choose_E(masterModes);
[W,R] = S.compute_whisker(param_order);
%% Compare Results

for i = 1:param_order
' ############# Order #############'

i
'SSM Coeffs and indices'

'W'
full(W(i).coeffs)
%'Wtilde'
%full(Wtilde(i).coeffs)
%full(Wtilde(i).ind.')
%full(phitilde4(i).coeffs)
% find indices
[tf, index]=ismember(phitilde(i).ind([1,3],:).',W(i).ind,'rows');
%'HNF Coeffs and indices'
%'Phi'
%full(phi(i).coeffs(:,tf))
%'Phitilde'
full(phitilde(i).coeffs(:,tf))
full(phitilde(i).ind(:,tf))
end