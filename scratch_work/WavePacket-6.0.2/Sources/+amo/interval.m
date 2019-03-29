%--------------------------------------------------------------------------
%
% Additional multiplicative operator:
% Characteristic function on an interval (in 1 dimension)
% Characteristic function on a rectangle (in 2 dimensions)
% etc ...
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%
% see the README file for license details.

classdef interval < amo.generic & handle
    
    properties (Access = public)
        min         % Beginning of interval
        max         % End of interval
    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = interval
        end
 
        % Initialize AMO: Set/check parameters
        function init_amo (obj)
            global space
            if length(obj.min)~=space.n_dim
                log.error ('Incompatible dimensionality for beginning of interval')
            end
            
            if length(obj.max)~=space.n_dim
                log.error ('Incompatible dimensionality for end of interval')
            end
            
            if any (obj.max <= obj.min)
                log.error ( 'Wrong ordering of interval min/max parameters' )
            end
            
            if isempty(obj.label)
                obj.label = 'Interval';
            end
            
        end
        
        % Display AMO, overloading default disp method
        function disp(obj)
            log.disp ( 'Interval / rectangle / cuboid / ...' )
            log.disp ( '***************************************************************' )
            log.disp (   ' ')
            log.disp ( [ 'Beginning of projection interval : ' num2str(obj.min) ] )
            log.disp ( [ 'End of projection interval       : ' num2str(obj.max) ] )
            log.disp ( [ 'Label                            : '         obj.label] )
        end
        
        % Evaluate AMO function
        function A = A(obj,r)
            global space
            
            A = ones(size(r{1}));
            
            % Tensor product of one-dimensional Gaussians
            for k = 1:space.n_dim
                A(r{k} < obj.min(k) ...
                    | r{k} > obj.max(k)) = 0;    
            end
            
        end
    end
end
