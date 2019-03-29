%--------------------------------------------------------------------------
%
% Generic properties of all additional multiplicative operators class definitions
%
% Also to be used when an additional multiplicative operators is not available
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt's group
%
% see the README file for license details.


classdef generic < handle
    
    properties (Access = public) 
        label       % Label of this AMO function
        ind         % Index of this AMO function
        
        dvr         % Grid representation (in N dimensions)
        
    end
  
    methods (Access = public)
        
        % Constructor
        function obj = generic
        end
        
        % Initialize AMO
        function init_amo (obj)           
        end
        
        % Display AMO, overloading default disp method
        function disp(obj)
            log.disp ('Not available')
            log.disp ('***************************************************************')
        end
        
        % Default: No AMO function available
        function A = A(obj,r)
            A = [];
        end

        % Grid representation of additional multiplicative operators
        function grid_amo (obj)
            global space
            
            obj.dvr = A ( obj, space.dvr );
            
        end
        
    end
 
end

