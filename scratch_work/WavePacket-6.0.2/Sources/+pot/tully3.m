%--------------------------------------------------------------------------
%
% Extended coupling example
%
% John C. Tully  
% Journal of Chemical Phyics 93, 1061 (1990)
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%
% see the README file for license details.

classdef tully3  < pot.generic & handle
    
    properties (Access = public)
        
        A
        B
        C
        
    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = tully3
            obj.A = 0.0006;
            obj.B = 0.10;
            obj.C = 0.90;
        end
        
        % Initialize potential: Set/check parameters
        function init_pot (~)
            
            global hamilt space
            
            if space.n_dim ~= 1
                log.error ('This potential is only for 1 dimension')
            end
            
            if hamilt.coupling.n_eqs ~= 2
                log.error ('This potential is only for 2 states')
            end
            
        end
        
        % Display potential, overloading default disp method
        function disp(obj)
            
            if obj.row==1 && obj.col==1
                log.disp ('Extended coupling example: J. of Chemical Phyics 93, 1061 (1990)')
                log.disp ('***************************************************************')
                log.disp (' ')
                log.disp ( [ 'Tully parameter A : ' num2str(obj.A) ] )
                log.disp ( [ 'Tully parameter B : ' num2str(obj.B) ] )
                log.disp ( [ 'Tully parameter C : ' num2str(obj.C) ] )
            else
                log.disp ('Same as above')
                log.disp ('***************************************************************')
            end
            
        end
       
        % Evaluate potential energy functions
        function V = V(obj,r)
            if obj.row==1 && obj.col==1
                V = + obj.A * ones(size(r{1}));
            elseif obj.row==2 && obj.col==2
                V = - obj.A * ones(size(r{1}));
            elseif obj.row==1 && obj.col==2
                left  = find ( r{1}< 0 );
                right = find ( r{1}>=0 );
                V(left)  = + obj.B *       exp (   obj.C * r{1}(left)  );
                V(right) = + obj.B * ( 2 - exp ( - obj.C * r{1}(right) ) );
                V = V';
            end
        end
        
        % Evaluate forces as negative gradients of potential
        function F = F(obj,r)
            if obj.row==1 && obj.col==1
                F{1} = zeros(size(r{1}));
            elseif obj.row==2 && obj.col==2
                F{1} = zeros(size(r{1}));
            elseif obj.row==1 && obj.col==2
                F{1}  = - obj.B * obj.C * exp ( - obj.C * abs(r{1})  );
            end
        end
        
    end

end


