% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init (~)
global hamilt plots space time

log.disp ( '*************************************' )
log.disp ( 'Gaussian packet in a Morse oscillator' )
log.disp ( '*************************************' )

% Spatial discretization
space.dof{1} = grid.fft;                 % using fft grid
space.dof{1}.mass = 1728.539;            % Reduced mass (OH molecule)
space.dof{1}.n_pts = 128;                % Number of grid points
space.dof{1}.x_min = 0.7;                % Lower bound of grid 
space.dof{1}.x_max = 7.0;                % Upper bound of grid

% Temporal discretization
time.steps.m_start  = 000;               % Index of initial time step
time.steps.m_stop   = 200;               % Index of final time step
time.steps.m_delta  = 005;               % Size of time steps 
time.steps.s_number = 020;               % Number of sub steps per time step

% Propagator
time.propa = tmp.splitting;              % Split operator method
time.propa.order = 2;                    % Strang splitting

% Initial wave function
time.dof{1}       = init.gauss;          % Gaussian-shaped wavepacket
time.dof{1}.width =  0.1;                % Width 
time.dof{1}.pos_0 =  1.44;               % Center in position representation
time.dof{1}.mom_0 =  0.0;                % Center in momentum representation

% Hamiltonian operator 
hamilt.truncate.e_min  =  0.0;           % Lower truncation of energy
hamilt.truncate.e_max  =  1.0;           % Upper truncation of energy

hamilt.pot{1,1}      = pot.morse;        % Harmonic oscillator
hamilt.pot{1,1}.d_e  = 0.1994;           % Dissociation energy
hamilt.pot{1,1}.r_e  = 1.821;            % Equilibrium length
hamilt.pot{1,1}.alf  = 1.189;            % Range parameter

% Absorbing boundary conditions
hamilt.nip{1}      = nip.power;          % Negativ imaginary potential
hamilt.nip{1}.exp  = 4;                  % Exponent
hamilt.nip{1}.min  = 1.0;                % Beginning of non-gobbler region
hamilt.nip{1}.max  = 6.0;                % End of non-gobbler region

% Plot time evolution of the density
plots.density       = vis.surface;       % Contour plot of Wigner transform

% Plot time evolution of expectation values
plots.expect          = vis.expect;
plots.expect.e_max    = 0.1;             % Manually set range for energy plot