clc, clear all
%% Initialise Hamiltonian
% example with cubic hamiltonian and a distance 3 nonlinearity with respect
% to the first mode. 
order = 4;  % Order of parametrisation
omega1 = 1;
omega2 = 1;
H(2).coeffs = [omega1, omega2];
H(2).ind    = [ 1    0      ; ...
                   0    1       ; ...
                   1    0       ;...
                   0    1       ];

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
               
if order >2
H(order +1) = polynomial_initialisation();
end               
%% Compute HNF
n = 2;
reduction = false;
[phi,phis,Hs,phicell] = HamiltonianNormalForm(H, order,n,reduction);
phi = flow_sparsification(phi);

%% SSM Computation
N = 2*n;
[A,B,F] = build_model_H(H,N);

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
[W,R] = S.compute_whisker(order);


%% Compare results
for i = 1:order
' ############# Order #############'

i
'SSM Coeffs and indices'

full(W(i).coeffs)

% find indices

if ~isempty(phi(i).ind)
%[tf, index]=ismember(phi(i).ind([1,3],:).',W(i).ind,'rows');
'HNF Coeffs and indices'
full(phi(i).coeffs)
full(phi(i).ind)
end
end