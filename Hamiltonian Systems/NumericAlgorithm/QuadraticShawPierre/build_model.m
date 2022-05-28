function [M,C,K,fnl] = build_model(m1,m2,k1,k2,gamma)
%% Makes dynamical system model of the quadratic Shaw Pierre

M = diag([m1,m2]);
K = [k1+k2  -k2 ;
      -k2   k1+k2];
C = zeros(2);

fnl(1).coeffs = [gamma;0];
fnl(1).ind    = [3,0];
end