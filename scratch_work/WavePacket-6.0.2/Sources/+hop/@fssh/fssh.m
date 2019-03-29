%--------------------------------------------------------------------------
%
% FSSH = Fewest switching surface hopping algorithm
% 
% Determine transition probability from d/dt |c(t)|^2
% from TDSEs attached to each of the trajectories.
% According to the proof given by Tully, this algorithm
% really results in the lowest number of switches.
% 
% see: J. C. Tully
%      J. Chem. Phys. 93(2), 1061-1071 (1990)
%      DOI:10.1063/1.459170
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt's group
%
% see the README file for license details.

classdef fssh < hop.mssh & handle
    
    methods (Access = public)
        
        % Display surface hopping details, overloading default disp method
        function disp(obj)
            
            % Inherit from superclass
            disp@hop.generic ( obj );
            
            % Logfile output
            log.disp ('Fewest switches surface hopping')
            log.disp (' ')
        end
        
        % Get probabilities of hopping from state "m" to state "n"
        % Eqs. (14,15) from doi:10.1063/1.459170
        function probable = prob_hop (obj,state,m,n,ind_m)
            global hamilt space time
            
            % Preallocate
            probable = zeros(size(ind_m));
            
            % Quantum coherence
            coherence = conj ( obj.psi{n}(ind_m) ) .* obj.psi{m}(ind_m);
            
            % Diabatic picture
            if strcmpi(hamilt.coupling.represent,'dia')
                coupling = state.ham {n,m}(ind_m);
                probable = +2 * imag (coupling .* coherence);
                
            % Adiabatic picture
            elseif strcmpi(hamilt.coupling.represent,'adi')
                coupling = zeros (size(probable));
                for d = 1:space.n_dim
                    coupling = coupling + ...
                        state.nac_1 {d}{n,m}(ind_m) .* state.mom{d}(ind_m) / space.dof{d}.mass;
                end
                probable = 2 * real (coupling .* coherence);
            end
            probable = time.steps.s_delta * probable ./ abs(obj.psi{m}(ind_m)).^2;
            
        end
        
    end
end

