% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008 Ulf Lorenz
%               2008-2009 Burkhard Schmidt
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init (~)
global hamilt plots space time

log.disp ( '*********************************************************************' )
log.disp ( 'Evolution of pendulum with squeezed initial function    ' )
log.disp ( '*********************************************************************' )

% Spatial discretization
space.dof{1}       = grid.fft;           % using fft grid
space.dof{1}.mass  = 1;                  % (Reduced) moment of inertia
space.dof{1}.n_pts = 64;                 % Number of grid points
space.dof{1}.x_min = 0;                  % Lower bound of grid 
space.dof{1}.x_max = 2*pi;               % Upper bound of grid

% Additional multiplicative operator
hamilt.amo{1}     = amo.cosine;          % Use cos-function as AMO
hamilt.amo{1}.exp = 1;                   % Degree of orientation

% Temporal discretization
time.steps.m_start = 000;                % Index of initial time step
time.steps.m_stop  = 100;                % Index of final time step

time.steps.m_delta = pi/50;              % Size of time steps
time.steps.s_number =  05;               % Number of sub steps per time step

% Propagator
time.propa = tmp.splitting;              % Split operator method
time.propa.order = 2;                    % Strang splitting

% Hamiltonian operator 
hamilt.truncate.e_min  = 000;            % Lower truncation of energy
hamilt.truncate.e_max  = 200;            % Upper truncation of energy

hamilt.pot{1,1}     = pot.pendulum;      % Plane pendulum
hamilt.pot{1,1}.eta = -50;               % Prefactor of cos-potential
hamilt.pot{1,1}.v_0 = +50;               % Energy offset

% Initial wave function
time.dof{1}          = init.pendulum1;   % Mathieu type wavefunction
time.dof{1}.parity   = 'c';              % cosine elliptic  
time.dof{1}.order    = 0;                % Order of Mathieu function 
time.dof{1}.multiple = 1;                % Potential multiplicity
time.dof{1}.barrier  = 100;              % Potential barrier
time.dof{1}.shift    = pi/1;             % Potential shift (horizontal)

% Plot densities
plots.density      = vis.polar;          % Contour plot of Wigner transform

% Plot expectation values
plots.expect       = vis.expect;
plots.expect.e_max = 100;                % Range for energy plot

