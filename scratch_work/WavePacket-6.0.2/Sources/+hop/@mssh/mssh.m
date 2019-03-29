%--------------------------------------------------------------------------
%
% MSSH = Multiple switches surface hopping (historic original)
%
% Determine transition probability from  |c(t)|^2
% from TDSEs attached to each of the trajectories.
% While in principle correct, this SH variant is 
% known to switch far too often, even when the
% trajectories are outside the regions of strong 
% non-adiabatic coupling
%
% see: J. C. Tully, R. K. Preston
%      J. Chem. Phys. 55(2), 562-572 (1971)
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt's group
%
% see the README file for license details.

classdef mssh < hop.generic & handle
    
    properties (Access = public)
         
        psi         % Associating quantum state vectors with each trajectory
        psi_old     % "Old" quantum state vectors for each trajectory
        psi_new     % "New" quantum state vectors for each trajectory

    end
    
    methods (Access = public)
        
        % Constructor: From superclass
        function obj = mssh
            obj = obj@hop.generic;
        end
        
        % Initialization
        function init_hop (obj,state)
            global hamilt

            % Inherit from initialization of superclass
            init_hop@hop.generic ( obj );
            
            % Preallocate quantum state vectors
            obj.psi     = cell(hamilt.coupling.n_eqs,1);
            obj.psi_new = cell(hamilt.coupling.n_eqs,1);
            obj.psi_old = cell(hamilt.coupling.n_eqs,1);
            
            % Set initial values for quantum state vectors
            for m=1:hamilt.coupling.n_eqs
                obj.psi{m}               = zeros (state.n_p, 1);
                obj.psi{m}(state.cha==m) = 1;
            end          

        end
        
        % Display surface hopping details, overloading default disp method
        function disp(obj)
            
            % Inherit from superclass
            disp@hop.generic ( obj );
            
            % Logfile output
            log.disp ('Multiple switches surface hopping')
            log.disp (' ')
        end
        
        % Preprocessing: before hopping
        function prep_hop ( obj,state,first_call )
            global hamilt
            
            % Inherit from superclass
            prep_hop@hop.generic ( obj,state,first_call );
            
            % Update Hamiltonians attached to Q/C trajectories
            eval_ham ( state, 1 )
            
            % Propagate quantum states
            if (strcmpi(hamilt.coupling.represent,'adi'))
                tdse_adi (obj,state) % adiabatic
            else
                tdse_dia (obj,state) % diabatic
            end
        end
            
        % Get probabilities of hopping from state "m" to state "n"
        function probable = prob_hop (obj,state,m,n,ind_m)
            
            % Multiple switches surface hopping (historic original)
            probable = abs(obj.psi{n}(ind_m)).^2;
            
        end
        
        % Tidying up after hopping
        function tidy_hop ( obj,state,m,n,allowed )
            % deliberately left empty
        end
                
        % see separate files for the following public methods
        tdse_adi (obj,state)
        tdse_dia (obj,state)

        
    end
    
end

