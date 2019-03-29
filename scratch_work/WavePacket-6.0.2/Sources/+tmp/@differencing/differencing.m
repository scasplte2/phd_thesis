%--------------------------------------------------------------------------
%
% Propagate objects of class wave (wavefunctions on grids) or traj (trajectory 
% bundles) by time.steps.s_number substeps of constant size time.steps.s_delta
% by repeated action of symmetrized finite differences 
%
% Classical dynamics 
% ==================
%
% Velocity Verlet (second order):
%                        
%    q(t+dt) = q(t) + dt * p(t)/m + dt^2 * F(q(t)) / (2*m) + O(dt^3)
% 
%    p(t+dt) = p(t) + dt * ( F(q(t)) + F(q(t+dt)) ) / 2 + O(dt^3)
%
% Stoermer-Verlet (third order):
%                                 F(q(t))   2        4
%    q(t+dt) = - q(t-dt) + 2 q(t) + ------- dt + O ( dt ) 
%                                    m
%
%    L. Verlet, Phys. Rev. 159, 98-103 (1967)  
%    http://dx.doi.org/10.1103/PhysRev.159.98
%
% Quantum dynamics (only second order available):
% ===============================================
%
%                               i      ^                3
%    psi(t+dt) = psi(t-dt) - 2 ---- dt H psi(t) + O ( dt ) 
%                              hbar              
%
%    A. Askar, A. S. Cakmak, J. Chem. Phys. 68, 2794 (1978)
%    http://dx.doi.org/10.1063/1.436072
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2008 Ulf Lorenz
%
% see the README file for license details.

classdef differencing < handle
    
    properties (Access = public)
        order       % order of method
    end
    
    properties (Access = private)
        pos_old     % Used for Stoermer Verlet algorithm
        pos_new     % Used for Stoermer Verlet algorithm
        frc_old     % Used for velocity Verlet algorithm
        % For future versions: add these private properties
        % psi_new
        % psi_old
    end
    
    methods (Access = public)
        
        % Construct propagator: Set default values
        function obj = differencing
            obj.order = 2;
        end
        
        % Display propagator, overloading default disp method
        function disp(obj)
            
            log.disp ('Finite differences, symmetrized in time           ')
            log.disp ('***************************************************************')
            log.disp (' ')
            
            switch obj.order
                case 2
                    log.disp ( 'Second order differencing aka Velocity-Verlet' )
                case 3
                    log.disp ( 'Third order differencing aka Stoermer-Verlet' )
                    log.disp ( 'Classical propagation only' )
            end
            
        end

        % Initialize QC propagator: Verlet, one time step backward
        function traj_init (obj, state)
            global space time

            % Check parameters
            switch obj.order
                case (2) % Velocity-Verlet propagator
                case (3) % Stoermer-Verlet propagator
                otherwise
                    log.error ( 'Wrong choice of splitting order for classical propagation' )
            end
            
            % Get initial energies
            apply_ham ( state );
            
            % Propagate positions one step backward in time (first order)
            for d = 1:space.n_dim
                obj.pos_old{d} = state.pos{d} - time.steps.s_delta * state.mom{d} / space.dof{d}.mass;
            end
            
            % Initialize hopping, if desired
            if isfield (time,'hop')
                init_hop (time.hop,state)
                disp     (time.hop)
            end
            
        end
        
        % Perform QC propagation: Verlet, one time step forward
        function traj_propa (obj, state)
            global space time
            
            % Loop over substeps
            for k = 1:time.steps.s_number
                
                % Save previous eigenvectors
                state.U_old = state.U_new;
                
                switch obj.order
                    case 2 % Velocity Verlet
                        
                        % Get (potentials and) forces
                        eval_V_F ( state )
    
                        % Propagate positions by one time step
                        for d = 1:space.n_dim
                            state.pos{d} = state.pos{d} + time.steps.s_delta     * state.mom{d} / space.dof{d}.mass;
                            state.pos{d} = state.pos{d} + time.steps.s_delta^2/2 * state.frc{d} / space.dof{d}.mass;
                        end
                        
                        % Get (potentials and) forces at new positions
                        obj.frc_old = state.frc;
                        eval_V_F ( state )
                        
                        % Propagate positions by one time step
                        for d = 1:space.n_dim
                            state.mom{d} = state.mom{d} + time.steps.s_delta/2 * (state.frc{d}+obj.frc_old{d}) / space.dof{d}.mass;
                        end
               
                        
                    case 3 % Stoermer Verlet
                        
                        % Get (potentials and) forces
                        eval_V_F ( state )
                        
                        % Propagate positions by one time step
                        % Get momenta by simple differencing
                        for d = 1:space.n_dim
                            obj.pos_new{d} = 2*state.pos{d} - obj.pos_old{d} + time.steps.s_delta^2 * state.frc{d} / space.dof{d}.mass;
                            state.mom{d} = ( obj.pos_new{d} - obj.pos_old{d} ) / (2*time.steps.s_delta) * space.dof{d}.mass;
                            obj.pos_old{d} = state.pos{d};
                            state.pos{d} = obj.pos_new{d};
                        end
                        
                end
                
                % Q/C surface hopping trajectories only
                if isfield (time,'hop')

                    % Optionally perform surface hopping
                    prep_hop (time.hop,state,k==1)            % Preprocessing
                    traj_hop (time.hop,state)                 % Trajectory hopping
                    
                end
                
            end % Loop over substeps
            
            % Get new energies
            apply_ham ( state );
            
        end
        
        % Initialize QM propagator
        function wave_init (obj, psi)
            
            % Check parameters
            switch obj.order
                case (2) % Second order differencing
                otherwise
                    log.error ( 'Wrong choice of splitting order for quantum propagation' )
            end
            
            % First step only: Initialize "old" wavefunction; no propagation yet
            psi.old = psi.dvr;
            
        end


        % Perform QM propagation
        function wave_propa (~, psi)
            global time hamilt
            
            % Loop over substeps
            for k = 1:time.steps.s_number
                
                % Perform propagation for one substep
                if isfield(time,'pulse')
                    apply_ham(psi,[ ...
                        time.efield.grid{1}(k+time.steps.offset) ...
                        time.efield.grid{2}(k+time.steps.offset)], 0);
                else
                    apply_ham(psi,[0 0], 0);
                end
                
                for m = 1:hamilt.coupling.n_eqs
                    psi.new{m} = psi.old{m} - 2 * 1i * time.steps.s_delta * psi.new{m};
                end
                
                % Get ready for next time step
                psi.old = psi.dvr;
                psi.dvr = psi.new;
                
                % Autocorrelation function (using position representation)
                time.steps.acf(k+time.steps.offset) = wave.braket(psi.ini,psi.dvr);
                
            end
        end
    end
end

