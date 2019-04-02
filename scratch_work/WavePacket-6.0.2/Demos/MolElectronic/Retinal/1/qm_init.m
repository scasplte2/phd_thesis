% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init (~)
global atomic hamilt plots space time

log.disp ( '*****************************************************************' )
log.disp ( 'Retinal isomerization: model in 1 dimension' )
log.disp ( 'S.Hahn, G.Stock, J. Phys. Chem. B 104(6),1149 (2000)' )
log.disp ( '*****************************************************************' )

% Number of (coupled) Schr�dinger equations
hamilt.coupling.n_eqs      = 2;
hamilt.coupling.represent  = 'adi';
hamilt.coupling.ini_rep    = 'adi';
hamilt.coupling.ini_coeffs = [0 1];      % Initially only upper adiabatic state populated

% Spatial discretization
space.dof{1}       = grid.fft;           % using FFT grid
space.dof{1}.mass  = 56198.347;          % effective mass
space.dof{1}.n_pts = 384;                % Number of grid points
space.dof{1}.x_min = -pi;                % Lower bound of grid 
space.dof{1}.x_max =  pi;                % Upper bound of grid

% Temporal discretization
time.steps.m_start  = 000;               % Index of initial time step
time.steps.m_stop   = 050;               % Index of final time step
time.steps.m_delta  = 100;               % Size of time steps 
time.steps.s_number = 100;               % Number of sub steps per time step

% Propagator
time.propa       = tmp.splitting;        % Operator splitting
time.propa.order = 2;                    % Strang splitting

% Hamiltonian operator 
hamilt.truncate.e_min    =  -0.1;        % Lower truncation of energy
hamilt.truncate.e_max    =  +0.5;        % Upper truncation of energy

hamilt.pot{1,1}         = pot.retinal;         % Hahn&Stock: Retinal potential
hamilt.pot{1,1}.shift   = 0.0  / atomic.E.eV;  % Vertical offset
hamilt.pot{1,1}.barrier = 3.6  / atomic.E.eV;  % Rotational barrier
hamilt.pot{1,1}.omega   = 0.19 / atomic.E.eV;  % Vibrational frequency
hamilt.pot{1,1}.kappa   = 0.0  / atomic.E.eV;  % Linear parameter
hamilt.pot{1,1}.xc      = 0.5;                 % Constant coupling coordinate

hamilt.pot{2,2}         = pot.retinal;         % Hahn&Stock: Retinal potential
hamilt.pot{2,2}.shift   =  2.48 / atomic.E.eV; % Vertical offset
hamilt.pot{2,2}.barrier = -1.09 / atomic.E.eV; % Rotational barrier
hamilt.pot{2,2}.omega   =  0.19 / atomic.E.eV; % Vibrational frequency
hamilt.pot{2,2}.kappa   = -0.10 / atomic.E.eV; % Linear parameter
hamilt.pot{2,2}.xc      =  0.5;                % Constant coupling coordinate

hamilt.pot{1,2}         = pot.retinal;         % Hahn&Stock: Retinal potential
hamilt.pot{1,2}.lambda  = 0.19 / atomic.E.eV;  % Linear coupling parameter
hamilt.pot{1,2}.xc      = 0.5;                 % Constant coupling coordinate

% Initial wave function
time.dof{1}        = init.gauss;         % Gaussian-shaped wavepacket
time.dof{1}.width  =  0.25;              % Width 
time.dof{1}.pos_0  =  0;                 % Center in position representation
time.dof{1}.mom_0  =  0;                 % Center in momentum representation

% Plot densities
plots.density        = vis.contour;      % Contour plot of Wigner transform
plots.density.energy = true;
plots.density.expect = false;

% Plot expectation values
plots.expect       = vis.expect;
plots.expect.e_max = 0.10;               % Range for energy plot