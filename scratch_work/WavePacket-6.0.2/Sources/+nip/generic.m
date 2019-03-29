%--------------------------------------------------------------------------
%
% Generic properties of all negative imaginary potential class definitions
%
% Also to be used when a negative imaginary potential is not available
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
        
        ind         % Index of (coupled) channel
        
        dvr         % Grid representation (in N dimensions)

    end

    methods (Access = public)
        
        % Constructor
        function obj = generic
        end
        
        % Initialize negative imaginary potential
        function init_nip (obj)
        end
        
        % Display potential, overloading default disp method
        function disp(obj)
            log.disp ('Not available' )
            log.disp ('***************************************************************' )
        end
        
        % Default: No negative imaginary potential available
        function W = W(obj,r)
            W = [];            
        end

        % Grid representation of negative imaginary potential
        function grid_nip (obj)
            global space
            
            obj.dvr = W ( obj, space.dvr );
            
        end
        
    end

end

