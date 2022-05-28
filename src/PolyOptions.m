classdef PolyOptions < matlab.mixin.SetGet
    %DSOptions Options for DynamicalSystems class
    properties
        notation  = 'multiindex' % 'multiindex'
        
    end
    methods
        function set.notation(obj,notation)
            switch lower(notation)
                case 'tensor'
                    obj.notation = 'tensor';
                case 'multiindex'
                    obj.notation = 'multiindex';
                otherwise
                    error('Unknown notation type: set tensor or multiindex notation types')
            end
        end
    end
end

