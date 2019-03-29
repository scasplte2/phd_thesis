%--------------------------------------------------------------------------
%
% Wave functions represented on DVR / FBR grids
% for use in fully quantum-mechanical dynamics
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2017-.... Burkhard Schmidt
%
% see the README file for license details.

classdef wave < handle
    
    properties (Access = public)
        
        dvr         % Discrete variable representation (cell vector for coupled states)
        adi         % Backup copy for adi=>dia transformation
        dia         % Backup copy for dia=>adi transformation
        fbr         % Finite basis representation (cell vector for coupled states)
        ini         % Initial wavefunction (DVR, cell vector)
        new         % New wavefunction (DVR, cell vector)
        old         % Old wavefunction (DVR, cell vector)
        sum         % Sum of wavefunctions (DVR, Chebychev only, cell vector)
        mom         % Momentum operator acting on wavefunction (DVR, single channel)
        
        redu        % Reduced densities (for plots only)
        wig         % Wigner transform (for plots only, cell vector)
        wig_max     % Maximum of Wigner transform
            
        sav_export  % Toggle export of wavefunction
        sav_dir     % Directory where to store the wave function
        sav_file    % File name template for the stored wave function
        sav_step    % Start a new file for saving wave functions every step main time steps
        sav_mem     % Start a new file for saving wave functions as soon as mem bytes have been filled
        
    end
    
    methods (Access = public)
        
        % Constructor: Console/logfile output and default values
        function obj = wave

            obj.sav_export = false;      % Toggle saving to binary
            obj.sav_dir    = pwd;        % Default directory 
            obj.sav_file   = 'wave';     % Default file name template
            obj.sav_step   = [];         % Start a new file for saving wave functions 
            obj.sav_mem    = 500 * 2^20; % Max. file size: 500 MB
            
            log.disp (' ')
            log.disp ('***************************************************************')
            log.disp ('Solving coupled Schroedinger equations using DVR/FBR techniques')
            log.disp ('***************************************************************')
            log.disp (' ')
            log.disp ('For quantum systems interacting with electrical fields              ')
            log.disp ('(using semiclassical dipole approximation)                          ')
            log.disp ('using partial differential equations (pde) solvers                  ')
            log.disp ('Atomic units (hbar = m_e = e = 1) are used throughout.              ')
            log.disp ('https://sf.net/p/wavepacket/matlab/wiki/Reference.Programs.qm_propa ')
            log.disp ('                                                                    ')
            log.disp ('H = T + V -iW -F*\mu -F^2/2*\alpha                                  ')
            log.disp ('                                                                    ')
            log.disp ('with T = T(-i d/dR)       Kinetic energy    (scalar)                ')
            log.disp ('with V = V(R)             Potential energy  (matrix)                ')
            log.disp ('                                                                    ')
            log.disp ('TDSE (qm_propa) only:                                               ')
            log.disp ('with W = W(R)             Negative imaginary potential (vector)     ')
            log.disp ('with F = F(t)             Electrical field  (scalar, time-dependent)')
            log.disp ('with \mu = \mu(R)         Dipole moment     (matrix)                ')
            log.disp ('with \alpha = \alpha(R)   Polarizability    (vector)                ')
            log.disp (' ')
            
        end
        
        % More methods: see separate files
        init_obj  ( obj )                % Initial conditions
        propagate ( obj, step )          % Propagation
        eigen     ( obj, step )          % Eigenfunctions of Hamiltonian
        observe   ( obj, step )          % Expectation values / uncertainties
        init_ham  ( obj, e, norm )       % Application of Hamiltonian
        apply_ham ( obj, e, norm )       % Application of Hamiltonian
        adiabatic ( obj, step, direction ) % Adiabatic<=>diabatic transformation
        diabatic  ( obj )                % Adiabatic=>diabatic (initial only)
        
    end
    
    methods (Static)
        
        wig = wigner (dvr)
        
        retval =   braket ( bra,           ket )
        retval = sandwich ( bra, operator, ket )
        
    end
    
end

