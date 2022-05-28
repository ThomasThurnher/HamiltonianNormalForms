clear all,clc, clear classes
% {
order =3;

%% Setup example Hamiltonian in original coordinates
% $$\mathcal{H} = \frac{\omega_1}{2} (q_1^2 +p_1^2 ) + \gamma q_1^3$$

omega1 = 1;
omega2 = 1.5;
gamma  = 1;
Hpq(2).coeffs = [omega1/2 omega1/2];
Hpq(2).ind    =  2* eye(2);

Hpq(order).coeffs = [ gamma ];
Hpq(order).ind    = [ order   ;...
                  0   ];
              

JA = [ 0    omega1  ;...
      -omega1  0 ];
[Hxy] = transform_Hamiltonian(JA,omega1,gamma,order);


[Ax,Bx,Fx] = build_model(Hxy,2);
[Aq,Bq,Fq] = build_model(Hpq,2);
full(Fx(2).coeffs)
Fx = transform_nl(Fx);
Fq = transform_nl(Fq);
%% Dynamical system setup

DSx = DynamicalSystem();
%DSx.CanonicalTrafo = true;
set(DSx,'A',Ax,'B',Bx,'F',Fx);
set(DSx.Options,'Emax',5,'Nmax',10,'notation','tensor')
set(DSx,'order',1);

[Vx,Dx,~] = DSx.linear_spectral_analysis();

%% 
% *Choose Master subspace (perform resonance analysis)*

Sx = SSM(DSx);
set(Sx.Options, 'reltol', 0.1,'notation','tensor')
masterModes = [1,2]; 
Sx.choose_E(masterModes);



%% Dynamical system setup

DSq = DynamicalSystem();
%DSq.CanonicalTrafo = true;
set(DSq,'A',Aq,'B',Bq,'F',Fq);
set(DSq.Options,'Emax',5,'Nmax',10,'notation','tensor')
set(DSq,'order',1);

[Vq,Dq,~] = DSq.linear_spectral_analysis();

%% 
% *Choose Master subspace (perform resonance analysis)*

Sq = SSM(DSq);
set(Sq.Options, 'reltol', 0.1,'notation','tensor')
masterModes = [1,2]; 
Sq.choose_E(masterModes);
[Wx,Rx] = Sx.compute_whisker(order-1);
[Wq,Rq] = Sq.compute_whisker(order-1);
% {



Wxcell = {Wx(1),Wx(order-1)};
Wqcell = {Wq(1),Wq(order-1)};

%Wq2x = compose_linear_flow(Trafo,Wqcell);
%Wq2x{2}.coeffs
Wx(order-1).coeffs
Wx2q = compose_linear_flow(Vq,Wxcell);

full(Wx2q{1}.coeffs)
full(Wq(1).coeffs)

Wx2q{2}.coeffs
full(Wq(order-1).coeffs)