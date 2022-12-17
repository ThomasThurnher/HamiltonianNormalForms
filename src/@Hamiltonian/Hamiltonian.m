classdef Hamiltonian < matlab.mixin.SetGetExactNames
    %SSM Spectral Submanifold Computation
    %   Detailed explanation goes here
    
    properties
        H          = [];       % struct for hamiltonian
        Hred       = [];       % reduced hamiltonian
        V          = [];       % symplectic transformation to diagonal coordinates
        degree     = [];       % degree of Hamiltonian
        n          = [];       % system size
        param_order= [];
        resModes   = [];       % resonant modes 
        
        tol        = 1e-10;    % tolerance for when a polynomial has zero coefficient
        reduction  = false;    % whether model reduction should be performed

        BBOptions = BBOptions();        
        Options = PolyOptions();
    end
    
    methods
        %% SET methods
        
        % set Hamiltonian
        function set.H(obj,P)
            order = numel(P);

            for i = 2:order
                if isempty(P(i))
                    H(i) = polynomial_initialisation();
                else
                    H(i) = P(i);
                end
            end

            if obj.param_order+1 > order
                 H(obj.param_order+1) = polynomial_initialisation();
            end 
            obj.H = H;
        end
        
        function set.n(obj,n)
            obj.n = n;
        end

       function set.resModes(obj,resModes)
            obj.resModes = resModes;
        end
        function set.param_order(obj,param_order)
            obj.param_order = param_order;
        end
        %% GET methods
        
        function degree = get.degree(obj)
            degree = numel(obj.H);
        end
        
        function n = get.n(obj)
            n = obj.n;
        end
        
        function Hred = get.Hred(obj)
            H = obj.H;
            
            for i = 1:numel(H)
                %try
                if ~isempty(H(i).ind)
                    ind = H(i).ind(obj.resModes,:);
                    
                    idx = find( sum(ind,1) == i );

                    Hred(i).coeffs = H(i).coeffs(idx);
                    Hred(i).ind    = ind(:,idx);
                %catch
                 %   Hred = []
                end
            end
        end
        %% other methods
           
        BB = extract_backbone(obj, modepair, omegaRange, order,W0,figs)               

        [H,V] = transform_Hamiltonian(obj);

        [A,B,F] = dynamicalSystem(obj);

        [phitilde,phi] = HamiltonianNormalForm(obj , order ,varargin)
    end
    
    methods (Access = protected)

    end
end

