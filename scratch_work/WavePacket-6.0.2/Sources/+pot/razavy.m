%--------------------------------------------------------------------------
% Potential energy for a Razayv single/double well
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2015 Burkhard Schmidt
%
% see the README file for license details.

classdef razavy < pot.generic & handle

    properties (Access = public)
        
        modified    % which version of Razavy potential
        
        s           % strength parameter
        n           % magic number: has to be integer for exact solvability
        
        eta         % prefactor of cosh
        zeta        % prefactor of cosh^2
        
    end
        
    methods (Access = public)
    
        % Constructor: Set default values
        function obj = razavy
            obj.modified = false;
        end
        
        % Initialize potential: Set/check parameters
        function init_pot (~)
            global space
            if space.n_dim ~= 1
                log.error ('This potential only for one dimension')
            end
         end
        

        % Display potential, overloading default disp method
        function disp(obj)
            log.disp ('Quasi-exactly solvable Razavy single/double well  ')
            log.disp ('***************************************************************')
            log.disp ('                                                  ')
            log.disp ('American Journal of Physics, 48(4), 285-288 (1980)')
            log.disp ('DOI:10.1119/1.12141, Eq. (2.2) with beta=1        ')
            log.disp ('                                                  ')
            log.disp (' V(R)= -s*(n+1)*cosh(2*x) +s^2/8*cosh(4*x) -s^2/8 ')
            log.disp ('                                                  ')
            log.disp (' or modified version (M. Mirahmadi/B. Schmidt)    ')
            log.disp (' Eur. Phys. J. D, 71 (6), 149 (2017)              ')
            log.disp ('                                                  ')
            log.disp (' V (R) = eta*cosh(x) + zeta*cosh^2(x)             ')
            log.disp ('                                                  ')
            if ~obj.modified
                log.disp ('Original version of Razayv single/double well potential')
                log.disp ( [ 'Strength parameter  s : ' num2str(obj.s) ] )
                log.disp ( [ 'Magic number n        : ' num2str(obj.n) ] )
            else
                beta  = sqrt(obj.zeta)*sign(obj.eta);
                kappa = abs(obj.eta/beta);
                e_min = -obj.eta^2/(4*obj.zeta);
                x_min = acosh(-obj.eta/(2*obj.zeta));
                e_max = obj.eta+obj.zeta;
                log.disp ('Modified version of Razayv single/double well potential')
                log.disp ( [ 'Prefactor of cosh    (eta) : ' num2str(obj.eta) ] )
                log.disp ( [ 'Prefactor of cosh^2 (zeta) : ' num2str(obj.zeta) ] )
                log.disp ( [ 'SUSY parameter      (beta) : ' num2str(    beta) ] )
                log.disp ( [ 'Topological index  (kappa) : ' num2str(   kappa) ] )
                log.disp ( [ 'Energy at minima           : ' num2str(   e_min) ] )
                log.disp ( [ 'Position of minima         : ' num2str(   x_min) ] )
                log.disp ( [ 'Energy at maximum          : ' num2str(   e_max) ] )
            end

        end

        % Evaluate potential energy function
        function V = V(obj,r)

            if ~obj.modified
                V = - obj.s*(obj.n+1) * cosh(2*r{1}) ...
                    + obj.s^2/8       * cosh(4*r{1}) ...
                    - obj.s^2/8;
            else   
                V = + obj.eta  * cosh(r{1}) ...
                    + obj.zeta * cosh(r{1}).^2; 
            end    
        end
        
        % Evaluate forces as negative gradients of potential
        function F = F(obj,r)
            log.error ('Code for calculation of forces still missing')
        end
        
    end
    
end
