function FRC = compute_output_actionangle2D(I0,psi0,stability,epsilon,Omega,W0,eta,nt, saveIC, outdof)
% COMPUTE_OUTPUT_POLAR2D: This function computes the periodic response in
% full physical coordinates from reduced polar coordintes in 2D and 
% returns the output data structure. 
FRC = struct('rho', cell(numel(I0),1), 'psi', [], 'stability', [], ...
        'Omega', [] , 'epsilon',epsilon,'Aout', [], 'Zout', [], 'Znorm', ...
        [], 'Zic', [], 'ZoutNorm', []);
    
    for k = 1:numel(I0) % loop over each fixed point
        % periodic response in reduced complex coordinates
        p = actionangle2complex(I0(k),psi0(k),eta,nt);
        % periodic response in Physical Coordinates
        z = reduced_to_full_traj(linspace(0,2*pi,nt),p,W0,[],epsilon,1);
        
        % extract output
        [FRC(k).rho, FRC(k).psi, FRC(k).Omega, FRC(k).stability] = ...
            deal(I0(k), psi0(k), Omega(k), stability(k));
        [FRC(k).Zout, FRC(k).Aout, FRC(k).Znorm, ...
            FRC(k).ZoutNorm] = extract_output(z, outdof);
        
        if saveIC
            FRC(k).Zic = z(:,1)'; % record initial condition on periodic orbit
        end
    end
end


function p = actionangle2complex(I,psi,m,nt)
%%
% Complex coordinates: $\mathbf{p}(\phi)=\rho_{0}\left[\begin{array}{c}e^{\iota(\psi_{0}+m\phi)}\\e^{-\iota(\psi_{0}+m\phi)}\end{array}\right],$
phi = linspace(0,2*pi,nt);
p = sqrt(  I*[exp(2i*psi + m*phi);
             exp(-2i*psi + m*phi)] );
end