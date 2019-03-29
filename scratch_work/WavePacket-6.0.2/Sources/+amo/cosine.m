%--------------------------------------------------------------------------
%
% Additional multiplicative operator:
% Rotational characteristics using cos^n as projector
% 
% n=1: Orientation
% n=2: Alignment
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2008 Ulf Lorenz
%
% see the README file for license details.

classdef cosine < amo.generic & handle
    
    properties (Access = public)
        exp         % Exponent of cosine function
        dof         % Which degree of freedom        
    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = cosine
            obj.exp = 1;
            obj.dof = 1;
        end
        
        % Initialize AMO: Set/check parameters
        function init_amo (obj)
            
            if isempty(obj.label)
                switch obj.exp
                    case 1
                        obj.label = 'Orientation';
                    case 2
                        obj.label = 'Alignment';
                    otherwise
                        obj.label = ['<cos^' int2str(obj.exp) '\theta>'];
                end
            end
            
        end
        
        % Display AMO, overloading default disp method
        function disp(obj)
            global space
            log.disp ( 'Directional property (cos^n) in one dimension' )
            log.disp ( '***************************************************************' )
            log.disp (' ')
            if space.n_dim > 1
                log.disp ( [ 'Degree of freedom: ' int2str(obj.dof)] )
            end
            log.disp ( [ 'Exponent n       : ' num2str(obj.exp)])
            log.disp ( [ 'Label            : ' obj.label])            
        end
        
        % Evaluate AMO function
        function A = A(obj,r)
            global space
            
            if isa(space.dof{obj.dof}, 'grid.legendre')
                A =     r{obj.dof} .^obj.exp;
            else
                A = cos(r{obj.dof}).^obj.exp;
            end
            
        end
    end
end

