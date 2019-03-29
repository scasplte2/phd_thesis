%--------------------------------------------------------------------------
%
% LZ_1: Landau-Zener formula, involving couplings
%
% see: A. K. Belyaev, O. V. Lebedev
%      Phys. Rev. A 84, 014701 (2011)
%      DOI:10.1103/PhysRevA.84.014701
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt's group
%
% see the README file for license details.

classdef lz_1 < hop.generic & handle
    
    properties (Access = public)
        gap_m1      % Energy gaps at previous time step
        gap_m2      % Energy gaps at pre-previous time step
    end
    
    methods (Access = public)
        
        % Constructor: From superclass
        function obj = lz_1
            obj = obj@hop.generic;
        end
        
        % Initialization
        function init_hop (obj,state)
            global hamilt
            
            % Inherit from initialization of superclass
            init_hop@hop.generic ( obj );
                       
            % Preallocate gap widths, criticality indices
            obj.gap_m1 = cell(hamilt.coupling.n_eqs,1);
            obj.gap_m2 = cell(hamilt.coupling.n_eqs,1);
            for m=1:hamilt.coupling.n_eqs
                for n=1:hamilt.coupling.n_eqs
                    obj.gap_m1{n,m} = zeros(state.n_p,1);
                    obj.gap_m2{n,m} = zeros(state.n_p,1);
                end
            end
        end
        
        % Display surface hopping details, overloading default disp method
        function disp(obj)
            
            % Inherit from superclass
            disp@hop.generic ( obj );
            
            % Logfile output
            log.disp ('Landau-Zener formula: 1st variant, involving couplings')
            log.disp (' ')
        end
        
        % Preprocessing: before hopping
        function prep_hop ( obj,state,first_call )
            
            % Inherit from superclass
            prep_hop@hop.generic ( obj,state,first_call );
            
            % Update Hamiltonians attached to Q/C trajectories
            eval_ham ( state, 1 )
            
        end
        
        % Get probabilities of hopping from state "m" to state "n"
        function probable = prob_hop (obj,state,m,n,ind_m)
            global hamilt space time
            
            % Preallocate
            probable = zeros(size(ind_m));

            % Eq. (1) from doi:10.1103/PhysRevA.84.014701
            if strcmpi(hamilt.coupling.represent,'dia')
                
                % Diabatic energy gap with first time derivatives
                gap    = state.ham {n,n}(ind_m) - state.ham {m,m}(ind_m);
                gap_d1 = (gap - obj.gap_m1{n,m}(ind_m) )/time.steps.s_delta;
                
                % Single switch citerion: Hopping can only occur at critical phase space points
                % Critical phase space points: Detect sign change of gap
                ind_c = find ((obj.gap_m1{n,m}(ind_m) > 0 & gap <0) | (obj.gap_m1{n,m}(ind_m) < 0 & gap > 0)) ;
                
                % Diabatic Landau Zener probability: only for "critical" trajectories
                coupling = state.ham {n,m}(ind_m);
                probable (ind_c) = exp ( - 2*pi * coupling(ind_c).^2 ./ abs(gap_d1(ind_c)) );
                
                % Getting ready for next iteration
                obj.gap_m1{n,m}(ind_m) = gap;
                
            % Eq. (3) from doi:10.1103/PhysRevA.84.014701 or 5 lines above
            % Eq. (3) from doi:10.1063/1.4882073 but why is gap squared there?!?
            elseif strcmpi(hamilt.coupling.represent,'adi')
                
                % Adiabatic energy gap
                gap = abs ( state.ham {n,n}(ind_m) - state.ham {m,m}(ind_m) );
                
                % Single switch citerion: Hopping can only occur at critical phase space points
                % Detect sign change of first derivative: from negative to positive
                ind_c = find (obj.gap_m2{n,m}(ind_m) > obj.gap_m1{n,m}(ind_m) & gap > obj.gap_m1{n,m}(ind_m));
                
                % Non-adiabatic coupling
                coupling = zeros (size(probable));
                for d = 1:space.n_dim
                    coupling = coupling + state.nac_1 {d}{n,m}(ind_m) .* state.mom{d}(ind_m) / space.dof{d}.mass;
                end
                
                % Adiabatic Landau Zener probability: only for "critical" trajectories
                probable (ind_c) = exp ( - pi/4 * gap(ind_c) ./ abs(coupling(ind_c)) );
                
                % Getting ready for next iteration
                obj.gap_m2{n,m}(ind_m)  = obj.gap_m1{n,m}(ind_m);
                obj.gap_m1{n,m}(ind_m)  =     gap        ;
                
            end
 
        end
        
        % Tidying up after hopping
        function tidy_hop ( obj,state,m,n,allowed )
            global hamilt
            
            % Reset gaps from previous (two!) time steps
            for k=1:hamilt.coupling.n_eqs
                obj.gap_m2{k,m}(allowed)  = 0;
                obj.gap_m1{k,m}(allowed)  = 0;
            end
 

        end
        
    end
end
 
