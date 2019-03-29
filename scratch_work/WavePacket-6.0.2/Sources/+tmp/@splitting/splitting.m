%-------------------------------------------------------------------------
%
% Operator splitting
% ==================
%
% Propagate objects of class wave (wavefunctions on grids) or traj (trajectory 
% bundles) by time.steps.s_number substeps of constant size time.steps.s_delta
% by (repeated action of) Trotter or Strang splitting 
%
% Classical dynamics 
% ===================
%
% (1) March frog (first order):
%
%    q(t+dt) = q(t) + dt * p(t)/m +
%    p(t+dt) = p(t) + dt * F(q(t+dt))
%
% (2) Leap frog (second order):
%
%    p(t+dt/2) = p(t     ) + dt * F(q(t)    )/2 
%    q(t+dt  ) = q(t     ) + dt *   p(t+dt/2)/m
%    p(t+dt  ) = p(t+dt/2) + dt * F(q(t+dt  ))/2 
%
% (3) Beeman (third order):
%
%     Reference: D. Beeman, J. Comp. Phys. 20, 130 (1976)
%     DOI:10.1016/0021-9991(76)90059-0
%
% (4)Yoshida (fourth order):
%
%     Reference: H. Yoshida, Phys. Lett. A 150, 262 (1990)
%     DOI:10.1016/0375-9601(90)90092-3
%
% Quantum dynamics 
% ================
%
%                                          F^2(t)          (       d)
%   H = V (R) - F(t) * ( D (R) + D (R) ) - ------ P(R) + T (R, -i --)
%                         p       t          2             (      dR)
%
%
% with kinetic operator T and potential V (R) energy (may be matrix-valued) 
% and where the interaction between an (optional) external electric 
% field F(t) and the dipole moment D(R) and/or the polarizability P(R)
% of the system is treated semiclassically. Here, the dipole matrix D 
% has been split into its diagonal (D_p) and offdiagonal (D_t) parts
% containing permanent and transition dipoles, respectively. In most
% cases, only one of the two is expected to play a role, depending on 
% the light frequencies under consideration. Moreover, note that F(t) 
% as well D_p(R), D_t(R), and P(R) can have two cartesian components 
% (along x,y) corresponding to different polarization directions.
%
% (1) Exact (for time-independent Hamiltonian)
%
% psi(t+tau) = exp(-i*H*tau  ) * psi(t) 
%
% (2) Lie-Trotter splitting
%
% psi(t+tau) = exp( -i*         V *tau ) 
%            * exp( +i*F(t)    *Dp*tau ) 
%            * exp( +i*F(t)    *Dt*tau ) 
%            * exp( +i*F^2(t)/2*P *tau ) 
%            * exp( -i*         T *tau )
%            * psi(t)                    + O(tau^2)
%
% (3) Strang-Marchuk splitting
%
% psi(t+tau) = exp( -i*         V *tau/2 ) 
%            * exp( +i*F(t)    *Dp*tau/2 ) 
%            * exp( +i*F(t)    *Dt*tau/2 ) 
%            * exp( +i*F^2(t)/2*P *tau/2 ) 
%            * exp( -i*         T *tau   ) 
%            * exp( +i*F^2(t)/2*P *tau/2 ) 
%            * exp( +i*F(t)    *Dt*tau/2 )
%            * exp( +i*F(t)    *Dp*tau/2 )
%            * exp( -i*         V *tau/2 ) 
%            * psi(t)                    + O(tau^3)
%
%
% Implementation for a single Schroedinger equation
% ------------------------------------------------
%
% The above Trotter and Strang formulae are straightforward to evaluate for
% FFT grids (plane wave DVR):
%   Operator V, and hence exp(V), are diagonal in position representation
%   Operator D, and hence exp(D), are diagonal in position representation
%   Operator T, and hence exp(T), are diagonal in momentum representation
% Hence, two FFTs have to be performed for every timestep to transform
% from position to momentum representation and back. Similar arguments
% apply for other DVR methods.
%
% Implementation for coupled Schroedinger equation
% -----------------------------------------------
%
% In case of coupled Schroeinger equations, the potential energy V as well
% as the dipole moment operator D may be represented by (real, symmetric) 
% matrices. Hence, for each spatial discretization point, the exponential 
% of these matrices has to be calculated. This can be achieved either by
% Matlab's expm function or by the method of eigenvectors and eigenvalues, 
% where the accuracy is determined by the condition of the matrix.
% 
% For the case of two equations, the matrix exponential can be calculated 
% analytically, see  Eq. (11.204) on page 316 of David Tannor's book. 
% 
%     ( alpha   beta )
% V = (              )
%     ( beta^* gamma )
%
% delta = (alfa-gamma)/2
%   eta = (alfa+gamma)/2
%   rho = sqrt(delta^2+beta^2)
%
% expm(-iV*tau) 
%   = exp(-i*eta*tau) *
%   * [              cos(tau*rho)*I 
%       -i*delta/rho*sin(tau*rho)*S3
%       -i* beta/rho*sin(tau*rho)*S1  ]
%
% where I=(1 0; 0 1), S1 = (0 1; 1 0), S3 = (1 0; 0 -1)
%
% 
% References
% ==========
%
% Original version: 
%     J. A. Fleck et al.,                  Appl. Phys.,   10,  129 (1976)
%     M. D. Feit, J. A. Fleck, A. Steiger, J. Comp. Phys. 47,  412 (1982) 
% Coupled equations: 
%     J. Alvarellos and H. Metiu,          J. Chem. Phys. 88, 4957 (1988)
%
%-------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%               2011 Ulf Lorenz
%
% see the README file for license details.

classdef splitting < handle
    
    properties (Access = public)
        order       % order of method
    end
    
    properties (Access = private)
        frac        % fraction: One (Trotter) or one half (Strang)
        frc_old     % Used for Beeman's algorithm
    end
    
    % For future versions: add these private properties (from global hamilt)
        % pot_couple
        % pot_expo{m,n}
        % dip_prod{p}{m}
        % pol_prod{p,q}{m}
        % dip_trans{p}
        % dip_eig_vals{p}{m}
        % dip_eig_vecs{p}{m,n}
        % dip_diagonal{p}{m}
    % end
    
    
    methods (Access = public)
        
        % Construct propagator: Set default values
        function obj = splitting
            obj.order = 2;  % Strang-Marchuk splitting
        end
              
        % Display propagator, overloading default disp method
        function disp(obj)

            log.disp ('Operator splitting        ')
            log.disp ('***************************************************************')
            log.disp (' ')
            
            switch obj.order
                case (1)
                    log.disp ( 'Lie-Trotter splitting (first order)' )
                    log.disp ( 'In classical mechanics a.k.a. march frog' )
                case (2)
                    log.disp ( 'Strang-Marchuk splitting (second order)' )
                    log.disp ( 'In classical mechanics a.k.a. leap frog' )
                case (3)
                    log.disp ( 'Beeman''s algorithm (third order)' )
                    log.disp ( 'For trajectories only' )
                case (4)
                    log.disp ( 'Yoshida''s algorithm (fourth order)' )
                    log.disp ( 'For trajectories only' )
            end
            
        end
        
        % Initialize QC propagator
        function traj_init (obj, state)
            global time
            
            % Check parameters
            switch obj.order
                case (1) % Lie-Trotter splitting a.k.a. "march frog"
                case (2) % Strang-Marchuk splitting a.k.a. "leap frog"
                case (3) % Beeman's algorithm
                case (4) % Yoshida integrator
                otherwise
                    log.error ( 'Wrong choice of splitting order for classical propagation' )
            end
            
            % Get initial energies
            apply_ham ( state );
            
            % Initialize hopping, if desired
            if isfield (time,'hop')
                init_hop (time.hop,state)
                disp     (time.hop)
            end
            
        end

        % Perform QC propagation
        function traj_propa (obj, state)
            global space time
            
            % Loop over substeps
            for k = 1:time.steps.s_number
                
                % Save previous eigenvectors
                state.U_old = state.U_new;
                
                switch obj.order
                    case 1 % Lie-Trotter splitting  a.k.a. "march frog"
                        
                        % Propagate positions by full time step
                        for d = 1:space.n_dim
                            state.pos{d} = state.pos{d} + time.steps.s_delta * state.mom{d} / space.dof{d}.mass;
                        end
                        
                        % Get potentials and forces at new positions
                        eval_V_F ( state )
                        
                        % Propagate momenta by full time step
                        for d = 1:space.n_dim
                            state.mom{d} = state.mom{d} + time.steps.s_delta * state.frc{d};
                        end
                        
                    case 2 % Strang-Marchuk splitting a.k.a. "leap frog"
                        
                        % Propagate momenta by half time step
                        for d = 1:space.n_dim
                            state.mom{d} = state.mom{d} + 0.5 * time.steps.s_delta * state.frc{d};
                        end
                        
                        % Propagate positions by full time step
                        for d = 1:space.n_dim
                            state.pos{d} = state.pos{d} + time.steps.s_delta * state.mom{d} / space.dof{d}.mass;
                        end
                        
                        % Get potentials and forces at new positions
                        eval_V_F ( state )
                        
                        % Propagate momenta by half time step
                        for d = 1:space.n_dim
                            state.mom{d} = state.mom{d} + 0.5 * time.steps.s_delta * state.frc{d};
                        end
                        
                    case 3 % Beeman's algorithm
                        
                        if(isempty(obj.frc_old)) % Initialization by Strang splitting
                            
                            % Propagate momenta by half time step
                            for d = 1:space.n_dim
                                state.mom{d} = state.mom{d} + 0.5 * time.steps.s_delta * state.frc{d};
                            end
                            
                            % Propagate positions by full time step
                            for d = 1:space.n_dim
                                state.pos{d} = state.pos{d} + time.steps.s_delta * state.mom{d} / space.dof{d}.mass;
                            end
                            
                            % Get potentials and forces at new positions
                            obj.frc_old = state.frc;
                            eval_V_F ( state )
                            
                            % Propagate momenta by half time step
                            for d = 1:space.n_dim
                                state.mom{d} = state.mom{d} + 0.5 * time.steps.s_delta * state.frc{d};
                            end
                            
                        else
                            
                            for d = 1:space.n_dim
                                state.pos{d} = state.pos{d} + time.steps.s_delta * state.mom{d} / space.dof{d}.mass ...
                                    + 1/6 * time.steps.s_delta^2 * (4 * state.frc{d} - obj.frc_old{d}) / space.dof{d}.mass;
                            end
                            
                            frc_current = state.frc;
                            eval_V_F ( state )
                            
                            for d = 1:space.n_dim
                                state.mom{d} = state.mom{d} + 1/6 * time.steps.s_delta * ...
                                    (2 * state.frc{d} + 5 * frc_current{d} - obj.frc_old{d});
                            end
                            
                            obj.frc_old = frc_current;
                            
                        end
                        
                    case 4 % 4th order Yoshida integrator
                        
                        coeff_c = zeros(4,1);
                        coeff_d = zeros(4,1);
                        
                        coeff_c([1,4]) = 1 / (2 * ( 2 - 2^(1/3) ) );
                        coeff_c([2,3]) = ( 1 - 2^(1/3) ) / (2 * ( 2 - 2^(1/3) ) );
                        
                        coeff_d([1,3])  = 1 / ( 2 - 2^(1/3) ) ;
                        coeff_d(2)      = - 2^(1/3) / ( 2 - 2^(1/3) ) ;
                        coeff_d(4)      = 0;
                        
                        for indices = 1:4
                            
                            % Propagate momenta
                            for d = 1:space.n_dim
                                state.mom{d} = state.mom{d} + coeff_c(indices) .* time.steps.s_delta * state.frc{d};
                            end
                            
                            % Propagate positions
                            for d = 1:space.n_dim
                                state.pos{d} = state.pos{d} + coeff_d(indices) .* time.steps.s_delta * state.mom{d} / space.dof{d}.mass;
                            end
                            
                            if( indices < 4)
                                % Get potentials and forces
                                eval_V_F ( state )
                            end
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
        function wave_init (obj, ~)
            global space time hamilt
            
            % Check/set parameters
            switch obj.order
                case (1) % Lie-Trotter splitting
                    obj.frac = 1;
                case (2) % Strang-Marchuk splitting
                    obj.frac = 1/2;
                otherwise
                    log.error ( 'Wrong choice of splitting order for quantum propagation' )
            end
            
            % Propagator for kinetic energy (all but the last one, see below)
            for k = 1:space.n_dim-1
                init_kin(space.dof{k}, obj.frac, false);
            end
            
            % External kinetic energies; Since they are typically the most
            % time-consuming operators, we put them in the middle of the split operator
            % If none are specified, take the last grid-internal kinetic energy
            % in the middle.
            if isfield(hamilt, 'kin')
                for n = 1:length(hamilt.kin)-1
                    init_kin(hamilt.kin{n}, obj.frac, false);
                end
                init_kin(hamilt.kin{end}, 1, false);
                init_kin(space.dof{end}, obj.frac, false);
            else
                init_kin(space.dof{end}, 1, false);
            end
            
            % Propagator for the potential energy
            pot_init (obj)
            
            % Propagator for (permanent and/or transition) dipoles, polarizabilities
            if isfield(time,'pulse')
                perm_init  (obj);
                trans_init (obj);
                pol_init   (obj);
            end
            
        end
        
        % Perform QM propagation
        function wave_propa (obj, psi)
            global space time hamilt
            
            % Loop over N sub-steps:
            for k = 1:time.steps.s_number
                
                % Potential energy: Full (Trotter) or half (Strang) sub-step
                pot_propa (obj, psi);
                
                % Permanent/transition dipoles, polarizabilities: Full (Trotter) or half (Strang) sub-step
                if isfield(time,'pulse')
                    e = [ ...
                        time.efield.grid{1}(k+time.steps.offset) ...
                        time.efield.grid{2}(k+time.steps.offset)];
                    perm_propa  (obj, psi, e);
                    trans_propa (obj, psi, e, k==1);
                    pol_propa   (obj, psi, e);
                end
                
                
                % Kinetic energy: Propagate each DOF
                for l = 1:space.n_dim
                    kinetic_exp(space.dof{l}, psi);
                end
                
                % Propagate each external kinetic energy
                if isfield(hamilt, 'kin')
                    for n = 1:length(hamilt.kin)
                        kinetic_exp(hamilt.kin{n}, psi);
                    end
                end
                
                % Strang-Marchuk splitting only
                if time.propa.order==2
                    
                    % Inverse order of the external kinetic energies
                    % And apply the last grid kinetic operator as well.
                    if isfield(hamilt, 'kin')
                        for n = length(hamilt.kin)-1:-1:1
                            kinetic_exp(hamilt.kin{n}, psi);
                        end
                        kinetic_exp(space.dof{end}, psi);
                    end
                    
                    % Apply all but the last grid kinetic energy operator
                    % (it has either been applied before, or it was the
                    % innermost split operator).
                    for l = space.n_dim-1:-1:1
                        kinetic_exp(space.dof{l},psi);
                    end
                    
                    % Transition/permanent dipoles, polarizabilities
                    if isfield(time,'pulse')
                        
                        % If possible, use field at next time step
                        if k+1+time.steps.offset < length (time.efield.grid{1})
                            e = [ ...
                                time.efield.grid{1}(k+1+time.steps.offset) ...
                                time.efield.grid{2}(k+1+time.steps.offset)];
                        else
                            e = [ ...
                                time.efield.grid{1}(k+time.steps.offset) ...
                                time.efield.grid{2}(k+time.steps.offset)];
                        end
                        pol_propa   (obj, psi, e);
                        trans_propa (obj, psi, e, true);
                        perm_propa  (obj, psi, e);
                    end
                    
                    % Potential energy
                    pot_propa (obj, psi);
                    
                end
                
                % Autocorrelation function (using position representation)
                time.steps.acf(k+time.steps.offset) = wave.braket(psi.ini,psi.dvr);
                
            end
            
        end
        
    end
    
    methods (Access = private)
        
        pot_init   (obj)
        pol_init   (obj)
        trans_init (obj)
        perm_init  (obj)

        pot_propa   (obj, psi)
        pol_propa   (obj, psi, e);
        trans_propa (obj, psi, e, recalc);
        perm_propa  (obj, psi, e);
                
    end
end
        
    
