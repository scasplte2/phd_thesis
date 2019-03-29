% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init (~)
global hamilt plots space time

log.disp ( '****************************************************************' )
log.disp ( 'Single crossing example, k=10 (Tully 1990)' )
log.disp ( '****************************************************************' )

% Number of (coupled) Schrödinger equations
hamilt.coupling.n_eqs      = 2;
hamilt.coupling.represent  = 'adi';
hamilt.coupling.ini_rep    = 'adi';
hamilt.coupling.ini_coeffs = [1 0];      % Initially only lower adiabatic state 

% Spatial discretization
space.dof{1} = grid.fft;                 % using FFT grid
space.dof{1}.mass = 2000;                % mass
space.dof{1}.n_pts = 256;                % Number of grid points
space.dof{1}.x_min = -10.0;              % Lower bound of grid 
space.dof{1}.x_max =  10.0;              % Upper bound of grid

% Temporal discretization
time.steps.m_start  = 000;               % Index of initial time step
time.steps.m_stop   = 020;               % Index of final time step
time.steps.m_delta  = 100.0;             % Size of time steps 
time.steps.s_number = 050;               % Number of sub steps per time step

% Propagator
time.propa = tmp.splitting;              % Split operator method
time.propa.order = 2;                    % Strang splitting

% Hamiltonian operator 
hamilt.truncate.e_min    = -0.01;        % Lower truncation of energy
hamilt.truncate.e_max    =  0.10;        % Upper truncation of energy

hamilt.pot{1,1} = pot.tully1;            % Single Crossing
hamilt.pot{2,2} = pot.tully1;            % Single Crossing
hamilt.pot{1,2} = pot.tully1;            % Single Crossing

% Absorbing boundary conditions
hamilt.nip{1} = nip.power;               % Negative imaginary potential
hamilt.nip{1}.exp  = 3;                  % Exponent
hamilt.nip{1}.min = -8.0;                % Beginning of inner grid region
hamilt.nip{1}.max = +8.0;                % End of inner grid region
hamilt.nip{2} = nip.power;               % Negative imaginary potential
hamilt.nip{2}.exp  = 3;                  % Exponent
hamilt.nip{2}.min = -8.0;                % Beginning of inner grid region
hamilt.nip{2}.max = +8.0;                % End of inner grid region

% Initial wave function
time.dof{1}       = init.gauss;          % Gaussian-shaped wavepacket
time.dof{1}.width =  0.75;               % Width 
time.dof{1}.pos_0 = -6.0;                % Center in position representation
time.dof{1}.mom_0 = 10.0;                % Center in momentum representation

% Plot densities
plots.density         = vis.contour;
plots.density.expect  = false;
plots.density.range   = true;            % manual setting of plotting range
plots.density.x_min   = -10;
plots.density.x_max   = 10;
plots.density.y_min   = 0;
plots.density.y_max   = 15;
plots.density.pot_max = 0.025;

% Plot expectation values
plots.expect       = vis.expect;
plots.expect.e_max = 0.04;
