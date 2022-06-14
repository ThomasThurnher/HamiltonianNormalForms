classdef Hamiltonian < matlab.mixin.SetGetExactNames
    %SSM Spectral Submanifold Computation
    %   Detailed explanation goes here
    
    properties
        H          = [];       % struct for hamiltonian
        V          = [];       % symplectic transformation to diagonal coordinates
        degree     = [];       % degree of Hamiltonian
        n          = [];       % system size
        param_order= [];

        tol        = 1e-10;    % tolerance for when a polynomial has zero coefficient
        reduction  = false;    % whether model reduction should be performed

        backboneOptions = backboneOptions();        
        Options = PolyOptions()
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
        %% other methods
           
        BB = extract_backbone(obj, parName, parRange, order)               

        [H,V] = transform_Hamiltonian(obj);

        [A,B,F] = dynamicalSystem(obj);

        [phitilde,phi] = HamiltonianNormalForm(obj , order ,varargin)
    end
    
    methods (Access = protected)

    end
end

