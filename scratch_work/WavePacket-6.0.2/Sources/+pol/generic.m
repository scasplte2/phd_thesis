%--------------------------------------------------------------------------
%
% Generic properties of all polarizability class definitions
%
% Also to be used when a polarizability is not available
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
        p_1         % polarization direction
        p_2         % polarization direction
        
        row         % row index (in diabatic representation)
        col         % column index (in diabatic representation)
        
        dvr         % Grid representation (in N dimensions)

    end
    
    methods (Access = public)
        
        % Constructor
        function obj = generic
        end
        
        % Initialize polarizability
        function init_pol (obj)
        end
        
        % Display polarizability, overloading default disp method
        function disp(obj)
            log.disp ('Not available')
            log.disp ('***************************************************************')
        end

        % No polarizability available
        function alpha = alpha(obj,r)
            alpha = [];
        end

        % Grid representation of polarizability
        function grid_pol (obj)
            global space
            
            obj.dvr = alpha ( obj, space.dvr );
            
        end
        
    end
    
end

