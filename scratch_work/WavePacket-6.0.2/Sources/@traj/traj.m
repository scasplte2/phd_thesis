%--------------------------------------------------------------------------
%
% Sampling densities by bundles of trajectories, i.e. swarms of points 
% in phase space, the evolution of which is governed by classical dynamics. 
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2017-.... Burkhard Schmidt
%
% see the README file for license details.

classdef traj < handle
    
    properties (Access = public)
        
        n_p         % number of phase space points/particles
        q_c         % variants of quantum-classical dynamics
        
        pos         % positions of particles (cell vector)
        mom         % momenta of particles (cell vector)
        
        frc         % forces acting on "particles" (cell vector)
        
        pot         % potential energy
        kin         % kinetic energy
        
        pot_mat     % diabatic potential matrices
        frc_mat     % diabatic forces
        
        cha         % Assigning channel indices to each trajectory
        ham         % Associating Hamiltonian matrices with each trajectory
        nac_1       % Associated first-order non-adiabatic coupling tensor

        D_new       % Adiabatic energies (eigenvalues)
        U_new       % Adiabatic states (eigenvectors) 
        U_old       % Adiabatic states (eigenvectors) of the previous step
        
        sav_export  % Toggle export of wavefunction
        sav_dir     % Directory where to store the wave function
        sav_file    % File name template for the stored wave function
        sav_step    % Start a new file for saving wave functions every step main time steps
        sav_mem     % Start a new file for saving wave functions as soon as mem bytes have been filled
        
        rnd_seed    % Toggle use of a predictable sequence of random numbers

        
    end
        
    methods (Access = public)
        
        % Constructor: Console/logfile output and default values
        function obj = traj (n,seed)
            
            % Number of phase space points/particles
            if nargin<1
                obj.n_p = 1000;          % number of phase space points/particles
            else
                obj.n_p = n;
            end
            
            % Predictable sequence of random numbers
            if nargin<2
                obj.rnd_seed = [];        
            else
                obj.rnd_seed = seed;
            end
                        
            obj.sav_export = false;      % Toggle saving to binary
            obj.sav_dir    = pwd;        % Default directory
            obj.sav_file   = 'traj';     % Default file name template
            obj.sav_step   = [];         % Start a new file for saving wave functions
            obj.sav_mem    = 500 * 2^20; % Max. file size: 500 MB
            
            log.disp (' ')
            log.disp ('***************************************************************')
            log.disp ('Solving the (quantum-)classical Liouville equation using ')
            log.disp ('sampling of phase-space densities by trajectory ensembles')
            log.disp ('***************************************************************')
            log.disp (' ')
            log.disp (['Number of particles : ' int2str(obj.n_p)])
            if isempty(obj.rnd_seed)
                log.disp ('Arbitrary sequence of random numbers')
            else
                log.disp (['Predictable sequence of random numbers; seed : ' int2str(obj.rnd_seed)])
                rng(obj.rnd_seed);
            end
            log.disp (' ')
            
        end
        
        % see separate files for the following public methods
        init_obj  ( obj )                % Initial conditions
        propagate ( obj, step )          % Propagation
        observe   ( obj, step )          % Expectation values and uncertainties
        init_ham  ( obj )                % Initialization of Hamiltonian
        apply_ham ( obj )                % Application of Hamiltonian
        eval_V_F  ( obj )                % Evaluate potential energy and forces 
        eval_ham  ( obj, nac )           % Evaluate Hamiltonian operators
        adiabatic ( obj, step, direction) % Adiabatic<=>diabatic transformation

    end
    
end

