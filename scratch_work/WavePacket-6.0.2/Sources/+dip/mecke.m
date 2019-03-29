%--------------------------------------------------------------------------
% Mecke model for (permanent) dipole moment
% see, e.g., Jakubetz, Manz, Mohan 
%            J. Chem. Phys. 90, 3686 (1989)    
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%
% see the README file for license details.

classdef mecke < dip.generic & handle
    
    properties (Access = public)
        
        q_0         % Charge parameter
        r_0         % Length parameter

    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = mecke
        end
        
        % Initialize dipole moment: Set/check parameters
        function init_dip (obj)
            global hamilt space
            
            if hamilt.coupling.n_eqs ~= 1
                log.error ('This polarizability is only for 1 state')
            end
            
            if space.n_dim > 1
                log.error ('This polarizability is only for 1 dof')
            end
        end
        
        % Display dipole moment, overloading default disp method
        function disp(obj)
            log.disp ('Mecke model function for permanent dipole moment     ')
            log.disp ('***************************************************************')
            log.disp ('                                                     ')
            log.disp ('  mu(r) = q_0 * r * exp ( - r / r_0 )                ')
            log.disp ('                                                     ')
            log.disp ('   see, e.g., Jakubetz, Manz, Mohan                  ')
            log.disp ('   J. Chem. Phys. 90, 3686 (1989)                    ')
            log.disp ('                                                     ')
            log.disp ( [ 'Charge parameter q_0: ' num2str(obj.q_0) ] )
            log.disp ( [ 'Length parameter r_0: ' num2str(obj.r_0) ] )
        end

        
        % Evaluate dipole moment
        function mu = mu(obj,r)

            mu = obj.q_0 * r{1}...
                .* exp ( - r{1} / obj.r_0 );

        end
    end
end
