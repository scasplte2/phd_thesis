% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init (~)
global hamilt plots space time

log.disp ( '*******************************************************************************' )
log.disp ( 'Free particle with momentum: Wavepacket translation and dispersion' )
log.disp ( '*******************************************************************************' )

% Spatial discretization
space.dof{1} = grid.fft;                 % using fft grid
space.dof{1}.mass = 1;                   % mass
space.dof{1}.n_pts = 096;                % Number of grid points
space.dof{1}.x_min = -10.0;              % Lower bound of grid 
space.dof{1}.x_max =  20.0;              % Upper bound of grid

% Temporal discretization and propagator
time.steps.m_start = 000;                % Index of initial time step
time.steps.m_stop  = 020;                % Index of final time step
time.steps.m_delta  = 0.5;               % Size of time steps 
time.steps.s_number       = 050;         % Number of sub steps per time step

time.propa = tmp.splitting;              % Split operator method
time.propa.order = 2;                    % Strang splitting

% Initial wave function
time.dof{1}       = init.gauss;          % Gaussian-shaped wavepacket
time.dof{1}.width = sqrt(1/2);           % Width 
time.dof{1}.pos_0 =  -5.0;               % Center in position representation
time.dof{1}.mom_0 =   1.0;               % Center in momentum representation

% Hamiltonian operator 
hamilt.truncate.e_min  =   0.0;          % Lower truncation of energy
hamilt.truncate.e_max  =  10.0;          % Upper truncation of energy

% Absorbing boundary conditions
hamilt.nip{1}      = nip.power;          % Negative imaginary potential
hamilt.nip{1}.exp  = 2;                  % Exponent
hamilt.nip{1}.min = -12;                 % Beginning of inner grid region
hamilt.nip{1}.max =  12;                 % End of inner grid region


% Plot densities
plots.density       = vis.contour;       % contour plot
plots.density.range = true;              % manual setting of plotting range
plots.density.x_min = -10;
plots.density.x_max = 10;
plots.density.y_min = -3;
plots.density.y_max = 3;

% Plot expectation values
plots.expect      = vis.expect;
plots.expect.e_max = 0.5;
