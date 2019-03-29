%--------------------------------------------------------------------------
%
% Generic properties of all system-bath coupling class definitions
%
% Also to be used when a system-bath coupling is not available
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
        
        row         % row index (in diabatic potential matrix)
        col         % column index (in diabatic potential matrix)
        
        dvr         % Grid representation (in N dimensions)

    end
    
    methods (Access = public)

        % Constructor
        function obj = generic
        end
        
        % Initialize system-bath coupling
        function init_sbc (obj)
        end
        
        % Display system-bath coupling, overloading default disp method
        function disp(obj)
            log.disp ('Not available')
            log.disp ('***************************************************************')
        end

        % No system-bath coupling available
        function chi = chi(obj,r)
            chi = [];
        end

        % Grid representation of system-bath coupling
        function grid_sbc (obj)
            global space
            
            obj.dvr = chi ( obj, space.dvr );
            
        end
        
    end
    
end

