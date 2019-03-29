%--------------------------------------------------------------------------
% Trigonometric (cosine-shaped) model for dipole moment
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2017 Burkhard Schmidt
%
% see the README file for license details.


classdef cosine < dip.generic & handle
    
    properties (Access = public)
        
        pre         % prefactor
        exp         % exponent
        
    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = cosine
            obj.pre = 1;
            obj.exp = 1;
        end
        
        % Initialize dipole moment: Set/check parameters
        function init_dip (obj)
            global hamilt space
            
            if hamilt.coupling.n_eqs ~= 1
                log.error ('This dipole function is only for 1 state')
            end
            
            if space.n_dim > 1
                log.error ('This dipole function is only for 1 dof')
            end
        end
        
        % Display dipole moment, overloading default disp method
        function disp(obj)
            log.disp (' mu(Theta) = f * cos^n  Theta                     ')
            log.disp ('***************************************************************')
            log.disp (' ')
            log.disp (['Prefactor f : ' num2str(obj.pre)])
            log.disp (['Exponent  n : ' num2str(obj.exp)])
        end
        
        % Evaluate dipole moment
        function mu = mu(obj,r)
            global hamilt space
            if isa (space.dof{hamilt.amo{1}.dof}, 'grid.legendre')
                mu = obj.pre *     r{1} .^obj.exp;
            else
                mu = obj.pre * cos(r{1}).^obj.exp;
            end
            
        end
    end
end
