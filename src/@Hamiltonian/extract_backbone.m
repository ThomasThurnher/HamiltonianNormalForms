function BB = extract_backbone(obj, modepair, omegaRange, order,W0,figs)

nI   = 50;
Imax = 0.2;


I = (Imax/nI) * (1:nI);

lambda = obj.H(2).coeffs(modepair)

freq = compute_freq(obj.Hred)




%{
for i = 1:numel(W0)
    try
        W0(i).ind = W0(i).ind.';
    catch
    end
end
%  EXTRACT_BACKBONE This function extracts the *Backbone curves in Polar coordinates.* For two-dimensional
% SSMs, we use the normal form of paramaterization, where we choose the following
% form of autonomous reduced dynamics as
%
% $$\mathbf{R}_{0}(\mathbf{p})=\left[\begin{array}{c}\lambda p\\\bar{\lambda}\bar{p}\end{array}\right]+\sum_{j=1}^{M}\left[\begin{array}{c}\gamma_{j}p^{j+1}\bar{p}^{j}\\\bar{\gamma}_{j}p^{j}\bar{p}^{j+1}\end{array}\right],$$
%
% Subsitution of $p=\rho e^{\mathrm{i}\theta}$ and $\bar{p}=\rho e^{-\mathrm{i}\theta}$
% to $\dot{\mathbf{p}}=\mathbf{R}_0(\mathbf{p})$ yields
%
% $$\dot{\rho}=a(\rho), \dot{\theta}=b(\rho)$$
%
% where
%
% $$a(\rho)=\sum_{j=1}^{M}\Re(\gamma_{j})\rho^{2j+1}+\rho\Re(\lambda),$$
%
% $$b(\rho)=\sum_{j=1}^{M}\Im(\gamma_{j})\rho^{2j+1}+\rho\Im(\lambda).$$
%
% It follows that the _backbone curves_ in polar coordinates is given by $\Omega=\frac{b(\rho)}{\rho}$.

%f1 = figure('Name','Norm');
%f2 = figure('Name',['Amplitude at DOFs ' num2str(obj.FRCOptions.outdof(:)')]);
%figs = [f1, f2];
colors = get(0,'defaultaxescolororder');

% get options
[nt, nI, nOmega, IScale, outdof, saveIC]  = ...
    deal(obj.BBOptions.nt, obj.BBOptions.nI, ...
    obj.BBOptions.nPar, obj.BBOptions.IScale, ...
    obj.BBOptions.outdof, obj.BBOptions.saveIC);
%% setup
startBB = tic;

lambda = obj.H(2).coeffs(modepair);

% some checks
assert(~isreal(lambda),'The eigenvalues associated to the modal subspace must be complex for analytic backbone computation')
omega0 = abs(imag(lambda));
assert(prod([omega0-omegaRange(1),omega0-omegaRange(end)])<0,'The supplied omegaRange must contain the natural frequency associated to the modes')

%% loop over orders
norders = numel(order);
for k=1:norders
    %% compute autonomous SSM coefficients
    
    set(obj,'resModes',[modepair, obj.n+modepair]);
    freq = compute_freq(obj.Hred)

    %% compute backbone
    I   = compute_I_grid(omegaRange,nOmega,IScale,freq,lambda,nI);
    omega = omega_bb(I, freq, lambda);

    idx = [find(omega<omegaRange(1)) find(omega>omegaRange(2))];
    I(idx) = []; omega(idx) = [];

    %% Backbone curves in Physical Coordinates

   
    stability = true(size(I)); psi = zeros(size(I)); epsilon = 0;
    BB = compute_output_actionangle2D(I,psi,stability,epsilon,omega,W0,1,nt, saveIC, outdof);
 
    %% plotting
    plot_FRC(BB,outdof,order(k),'freq','lines',figs,colors(k+1,:));
    
    
    I
    omega    
    figure()
    plot(omega,I)
end
totalComputationTime = toc(startBB);
disp(['Total time spent on backbone curve computation = ' datestr(datenum(0,0,0,0,0,totalComputationTime),'HH:MM:SS')])

%}
end


function freq = compute_freq(H)
j = 4;

while j <= numel(H)
    Hj = H(j);
    Hj
    1/1i* j/2 * Hj.coeffs
    freq(j/2-1) = 1/1i* j/2 * Hj.coeffs;
    j = j+2;
end

end

function I = compute_I_grid(omegaRange,nOmega,IScale,freq,lambda,nI)
omega = linspace(omegaRange(1),omegaRange(end), nOmega);
%%
% *Explicit quadratic approximation of the backbone curve*
%
% $$\rho_{backbone} = \sqrt{\frac{\Omega-\Im(\lambda)}{\Im(\gamma_1)}}$$

I_bb = (abs((omega - lambda)/ freq(1)));
Imax = IScale * max(I_bb);
I = (Imax/nI) * (1:nI);

end