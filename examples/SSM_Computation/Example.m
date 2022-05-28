clear all,clc, clear classes

%% Setup example $

omega1 = 1;
gamma  = 1;
order = 2;
A = [ 0    omega1  ;...
      -omega1  0 ];

B = eye(2);

F(order).coeffs = [0;-(order+1)*gamma];
F(order).ind    = [order,0];


[Aq,Bq,Fq] =deal(A,B,F);


%% Dynamical system setup

DSq = DynamicalSystem();
DSq.CanonicalTrafo = true;
set(DSq,'A',Aq,'B',Bq,'F',Fq);
set(DSq.Options,'Emax',5,'Nmax',10,'notation','multiindex')
set(DSq,'order',1);

[Vq,Dq,~] = DSq.linear_spectral_analysis();

%% 
% *Choose Master subspace (perform resonance analysis)*

Sq = SSM(DSq);
set(Sq.Options, 'reltol', 0.1,'notation','multiindex')
masterModes = [1,2]; 
Sq.choose_E(masterModes);%% Dynamical system setup

%%

[Ax,Bx,Fx] = transform_system(A,B,F,Vq);

[~,~,Finv] = transform_system(A,B,Fx,inv(Vq));


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

%%
[Wx,Rx] = Sx.compute_whisker(order);
[Wq,Rq] = Sq.compute_whisker(order);
% {



Wxcell = {Wx(1),Wx(order)};
Wqcell = {Wq(1),Wq(order)};

%Wq2x = compose_linear_flow(Trafo,Wqcell);
%Wq2x{2}.coeffs
Wx(order).coeffs
Wx(1).coeffs
Wx2q = compose_linear_flow(Vq,Wxcell);

full(Wx2q{1}.coeffs)
full(Wq(1).coeffs)

Wx2q{2}.coeffs
full(Wq(order).coeffs)