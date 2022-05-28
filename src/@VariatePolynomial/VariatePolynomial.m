classdef VariatePolynomial < matlab.mixin.SetGetExactNames
    % Construct a variate polynomial
    
    properties

        P                   % struct containing the polynomial
        N                   % dimension of space over which polynomial is defined
        degree              % degree of polynomial
        
        Options = PolyOptions()
        CanonicalTrafo = false;
    end
    
    methods
        %% SET methods
        
        % set polynomial
        function set.P(obj,P)
            obj.P = P;
        end
        
        
        %% GET methods
        function N = get.N(obj)
            
            N = [];
            i = 1;
            while isempty(N)
                try                   
                    N = size(obj.P(i).ind,1);                    
                catch
                    i = i+1;
                end
            end
            
        end
        
        
        %% other methods
        
        [p_revlex] = polynomial_2_revlex(obj,p);
        
        [PmF] = polynomial_addition(obj,P,F);
        
        [Pcoll] = polynomial_collapse(obj,P, varargin);
        
        [DjP] = polynomial_derivative(obj,P,j);
        
        [DP] = polynomial_gradient(obj,P);
        
        [P] = polynomial_initialisation(obj);
        
        [PF] = polynomial_multiplication(obj,P,F);

        [PmF] = polynomial_subtraction(obj,P,F);
        
        [JDP] = polynomial_vectorfield(obj,P);
    end
end

