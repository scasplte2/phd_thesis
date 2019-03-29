%--------------------------------------------------------------------------
%
% LZ_2: Landau-Zener formula, involving gaps only
%
% see: A. K. Belyaev, O. V. Lebedev
%      Phys. Rev. A 84, 014701 (2011)
%      DOI:10.1103/PhysRevA.84.014701
%
% see: A. K. Belyaev, C. Lasser, G. Trigila
%      J. Chem. Phys. 140, 224108 (2014)
%      DOI:10.1063/1.4882073
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt's group
%
% see the README file for license details.

classdef lz_2 < hop.lz_1 & handle
    
    methods (Access = public)
        
        % Initialization
        function init_hop (obj,state)
            global hamilt
            
            % Inherit from initialization of superclass
            init_hop@hop.lz_1 ( obj,state );
            
            % LZ_2 not available in diabatic representation
            if strcmpi(hamilt.coupling.represent,'dia')
                log.error ('LZ_2 not available in diabatic representation')
            end
            
        end
        
        % Display surface hopping details, overloading default disp method
        function disp(obj)
            
            % Inherit from superclass
            disp@hop.generic ( obj );
            
            % Logfile output
            log.disp ('Landau-Zener formula: 2nd variant, involving gaps only!')
            log.disp (' ')
        end
        
        % Get probabilities of hopping from state "m" to state "n"
        % Eq. (2) from doi:10.1103/PhysRevA.84.014701
        % Eq. (4) from doi:10.1063/1.4882073
        function probable = prob_hop (obj,state,m,n,ind_m)
            global time
            
            % Preallocate
            probable = zeros(size(ind_m));
            
            % Adiabatic energy gap
            gap    = abs ( state.ham {n,n}(ind_m) - state.ham {m,m}(ind_m) );
            
            % Single switch citerion: Hopping can only occur at critical phase space points
            % Detect sign change of first derivative: from negative to positive
            ind_c = find (obj.gap_m2{n,m}(ind_m) > obj.gap_m1{n,m}(ind_m) & gap > obj.gap_m1{n,m}(ind_m));
            
            % Local minimum adiabatic energy gap (previous step) of the trajectories
            % with second time derivatives
            % https://en.wikipedia.org/wiki/Finite_difference_coefficient
            gap_min     = obj.gap_m1{n,m}(ind_m(ind_c));
            gap_min_d2  = (gap(ind_c) - 2* gap_min + obj.gap_m2{n,m}(ind_m(ind_c))) / time.steps.s_delta^2;
            
            % Adiabatic Landau Zener probability: only for "critical" trajectories
            probable (ind_c) = exp ( - pi/2 * sqrt(gap_min.^3 ./ gap_min_d2 ) );
            
            % Getting ready for next iteration
            obj.gap_m2{n,m}(ind_m)  = obj.gap_m1{n,m}(ind_m);
            obj.gap_m1{n,m}(ind_m)  =     gap        ;
            
        end
    end
    
end

