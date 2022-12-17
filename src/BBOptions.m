classdef BBOptions < matlab.mixin.SetGet
    %DSOptions Options for DynamicalSystems class
    properties
        markersize  = 10           % 'multiindex'
        
        nI = 20                 % number of discrete rho values in the range [0, rhomax] (relevant for method == 'level set') 
        nPar = 10                 % number of discrete parameter (Omega/epsilon) values in the parameter range 
        nt = 7                   % number of discrete time intervals over which the periodic orbit is discretized (relevant for post-processing only)
        IScale = 1               % factor for increasing rhomax polar FRC
        

        saveIC = true              % whether save initial conditions on periodic orbit or not 
        outdof = []                % output degree-of-freedom  
        
       
    end
    methods

    end
end

