%--------------------------------------------------------------------------
%
% Representing quantum states by ket vectors 
% based on an eigen|energy representation
%
% Should we introduce a super class ABNCD ?!?
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt
%
% see the README file for license details.

classdef ket < handle
    
    properties (Access = public)
        
        x           % input (state vector)
        y           % output (observables)
        
        A           % matrix;   from energies 
        B           % vectors;  from transition dipole moments
        N           % matrices; from transition dipole moments
        D           % matrices; from (quadratic) observables
        
    end
    
    methods (Access = public)
        
        % Constructor: Console/logfile output and default values
        function obj = ket ( )
            
            log.disp ('                                                               ')
            log.disp ('***************************************************************')
            log.disp ('Matrices A, N, B and C, D and input x and output y             ')
            log.disp ('***************************************************************')
            log.disp ('                                                               ')
            log.disp ('for use in a bilinear control problem                          ')
            log.disp ('https://sourceforge.net/p/wavepacket/matlab/wiki/Reference.Programs.qm_abncd/')
            log.disp ('                                                               ')
            log.disp ('         d                                                     ')
            log.disp ('control: -- x(t) = ( A + iu(t)N ) x(t) + iu(t)B                ')
            log.disp ('         dt                                                    ')
            log.disp ('                 T                                             ')
            log.disp ('observe: y(t) = x (t) D x(t)                                   ')
            log.disp ('                                                               ')
            log.disp ('see B. Schaefer-Bung, C. Hartmann, B. Schmidt, Ch. Schuette    ')
            log.disp ('J. Chem. Phys. 135, 014112-1-13 (2011)                         ')            
            log.disp ('***************************************************************')
            log.disp ('                                                       ')
            log.disp (' coming from the time-dependent Schoedinger equation   ')
            log.disp (' using atomic units throughout, i.e. hbar = m_e = e = 1')
            log.disp (' for quantum systems interacting with electrical fields')
            log.disp ('                                                       ')
            log.disp ('    d                                                  ')
            log.disp (' i -- |psi(t)> = ( H  - F(t) mu ) |psi(t)>             ')
            log.disp ('   dt               0                                  ')
            log.disp ('                                                       ')
            log.disp (' H_0 is a diagonal matrix for the unperturbed system,  ')
            log.disp (' mu is the Hermitian matrix with dipole moments        ')
            log.disp (' and F(t) is the electric field.                       ')
            log.disp ('                                                       ')
            log.disp ('  <D>(t)= <psi(t)| D |psi(t)>                          ')
            log.disp ('                                                       ')
            
        end
        
        % see separate files for the following public methods
        init_obj  ( obj )                % Initial conditions
        progagate ( obj, step )          % Solving the TDSE using ODE methods
        observe   ( obj, step )          % Expectation values from CD
        init_ham  ( obj )                % Initialization of ABNCD
        apply_ham ( obj )                % Application of ABN

    end
end

