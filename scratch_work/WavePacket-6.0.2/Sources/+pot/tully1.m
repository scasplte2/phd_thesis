%--------------------------------------------------------------------------
%
% Single crossing example
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

classdef tully1 < pot.generic & handle
    
    properties (Access = public)
        
        A           % Diabatic energy parameter
        B           % Diabatic range parameter
        C           % Coupling energy parameter
        D           % Coupling range parameter
                
    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = tully1
            obj.A = 0.01;
            obj.B = 1.6;
            obj.C = 0.005;
            obj.D = 1.0;  
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
                log.disp ('Single crossing example: J. of Chemical Phyics 93, 1061 (1990) ')
                log.disp ('***************************************************************')
                log.disp (' ')
                log.disp ( [ 'Tully energy parameter A : ' num2str(obj.A) ] )
                log.disp ( [ 'Tully range parameter  B : ' num2str(obj.B) ] )
                log.disp ( [ 'Tully energy parameter C : ' num2str(obj.C) ] )
                log.disp ( [ 'Tully range parameter  D : ' num2str(obj.D) ] )
            else
                log.disp ('Same as above')
                log.disp ('***************************************************************')
            end
        end
        
        % Evaluate potential energy functions
        function V = V(obj,r)
            if obj.row==1 && obj.col==1
                V = + obj.A * sign(r{1}) .* ( 1 - exp(-obj.B*sign(r{1}).*r{1}) );
            elseif obj.row==2 && obj.col==2
                V = - obj.A * sign(r{1}) .* ( 1 - exp(-obj.B*sign(r{1}).*r{1}) );
            elseif obj.row==1 && obj.col==2
                V = obj.C * exp ( - obj.D * r{1} .^ 2);
            end
        end
        
        % Evaluate forces as negative gradients of potential
        function F = F(obj,r)
            if obj.row==1 && obj.col==1
                F{1} = - obj.A * obj.B .* exp(-obj.B*sign(r{1}).*r{1});
            elseif obj.row==2 && obj.col==2
                F{1} = + obj.A * obj.B .* exp(-obj.B*sign(r{1}).*r{1});
            elseif obj.row==1 && obj.col==2
                F{1} = 2 * obj.C * obj.D .* r{1} .* exp ( - obj.D * r{1} .^ 2);
            end
        end

    end
    
end


