%--------------------------------------------------------------------------
%
% Generic properties of all dipole moment class definitions
%
% Also to be used when a dipole moment is not available
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
        pol         % polarization direction
        
        row         % row index (in diabatic representation)
        col         % column index (in diabatic representation)
        
        dvr         % Grid representation (in N dimensions)

    end
   
    methods (Access = public)
        
        % Constructor
        function obj = generic
        end
        
        % Initialize dipole moment
        function init_dip (obj)
        end
        
        % Display dipole moment, overloading default disp method
        function disp(obj)
            log.disp ('Not available')
            log.disp ('***************************************************************')
        end

        % Default: No dipole moment available
        function mu = mu(obj,r)
            mu = [];
        end

        % Grid representation of dipole moment
        function grid_dip (obj)
            global space
            
            obj.dvr = mu ( obj, space.dvr );
            
        end
        
    end

end

