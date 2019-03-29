%--------------------------------------------------------------------------
%
% Negative imaginary potential (NIP) to be used  
% as a (smoothly) absorbing boundary condition
%
% Power function: (x-x0)^n
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2011 Ulf Lorenz
%
% see the README file for license details.

classdef power < nip.generic & handle
    
    properties (Access = public)
        exp         % Exponent of power function 
        min         % Beginning of inner grid region
        max         % End of inner grid region
        eps         % Small parameter

    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = power
            
            % Small parameter: exp(-NIP)=eps at grid boundaries
            obj.eps = eps('single');
            
        end
        
        % Initialize potential: Set/check parameters
        function init_nip (obj)
            global space
            
            if length(obj.exp)~=space.n_dim
                log.error ('Incompatible dimensionality for NIP exponenents')
            end
            
            if length(obj.min)~=space.n_dim
                log.error ('Incompatible dimensionality for beginning of inner grid')
            end
            
            if length(obj.max)~=space.n_dim
                log.error ('Incompatible dimensionality for end of inner grid')
            end
            
            if any (obj.exp < 1)
                log.error ( 'Exponent of negative imaginary potential must be at least 1' )
            end
            
            if any (obj.max <= obj.min)
                log.error ( 'Wrong ordering of "NIP" min/max parameters' )
            end
        end
        
        % Display potential, overloading default disp method
        function disp(obj)
            log.disp ( 'Power function' )
            log.disp ( '***************************************************************' )
            log.disp ( ' ' )
            log.disp ( [ 'Exponent of power function     : ' int2str(obj.exp) ] )
            log.disp ( ' ' )
            log.disp ( [ 'Beginning of inner grid region : ' num2str(obj.min) ] )
            log.disp ( [ 'End of inner grid region       : ' num2str(obj.max) ] )
        end
        
        
        % Evaluate potential energy function
        function W = W(obj,r)
            
            global space
            
            % Initialize grid representation
            W = zeros(size(r{1}));
            
            % Adding up contributions from each dimension
            for k = 1:space.n_dim
                
                % Left/lower "gobbler" region
                mask = find ( space.dvr{k} < obj.min(k) );
                if ~isempty ( mask )
                    factor = -log(obj.eps)/(obj.min(k)-space.dof{k}.dvr_min)^obj.exp(k);
                    W(mask) = W(mask) + factor * ...
                        (obj.min(k)-r{k}(mask)).^obj.exp(k);
                end
                
                % Right/upper "gobbler" region
                mask = find ( space.dvr{k} > obj.max(k) );
                if  ~isempty ( mask )
                    factor = -log(obj.eps)/(space.dof{k}.dvr_max-obj.max(k))^obj.exp(k);
                    W(mask) = W(mask) + factor * ...
                        (r{k}(mask)-obj.max(k)).^obj.exp(k);
                end
                
            end
            
        end
    end
end
