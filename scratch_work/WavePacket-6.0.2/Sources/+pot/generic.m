%--------------------------------------------------------------------------
%
% Generic properties of all potential energy class definitions
%
% Also to be used when a potential energy function is not available
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
        dia         % Backup copy for dia=>adi transformation

    end
    
    methods (Access = public)

        % Constructor
        function obj = generic
        end
        
        % Initialize potential
        function init_pot (obj)
        end
        
        % Display potential, overloading default disp method
        function disp(obj)
            log.disp ('Not available')
            log.disp ('***************************************************************')
        end
   
        % Default: No potential available
        function V = V(obj,r)
            V = [];
        end
        
        % Default: No forces available
        function F = F(obj,r)
            F = cell(size(r));
            for d=1:length(r)
                F{d} = [];
            end
        end
        
        % Grid representation of potential energy function
        function grid_pot (obj)
            global space
            
            obj.dvr = V ( obj, space.dvr );
            
        end
        
    end
    
end

