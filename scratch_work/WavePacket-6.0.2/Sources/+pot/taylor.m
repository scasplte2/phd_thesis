%--------------------------------------------------------------------------
% Represent the potential energy function as a Taylor series
% Includes free particle (V=0) as a special case
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2016-2017 Burkhard Schmidt
%
% see the README file for license details.

classdef taylor < pot.generic & handle
    
    properties (Access = public)
        
        hshift      % Horizontal shift (scalar)
        vshift      % Vertical shift (row vector)
        coeffs      % Coefficients, i.e. derivatives
        
    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = taylor
            obj.hshift  = [];
            obj.vshift  = 0;
            obj.coeffs  = [];
        end
        
        % Initialize potential: Set/check parameters
        function init_pot (obj)
            global space
            if isempty(obj.hshift)
                obj.hshift = zeros(1,space.n_dim);
            end
        end
        
        % Display potential, overloading default disp method
        function disp(obj)
            if obj.vshift==0 && isempty (obj.coeffs)
                log.disp ('Not available ("free particle")                   ')
                log.disp ('***************************************************************')
            else
                log.disp ('Taylor series (diag. in N dimensions)')
                log.disp ('***************************************************************')
            end
        end
        
        % Evaluate potential energy function
        function V = V(obj,r)  
            V = math.taylor (...
                r, ...
                obj.hshift, ...
                obj.vshift, ...
                obj.coeffs, 0 );
        end
        
        % Evaluate forces as negative gradients of potential
        function F = F(obj,r)
            F = math.taylor_d (...
                r, ...
                obj.hshift, ...
                obj.coeffs );
        end
        
    end
end
