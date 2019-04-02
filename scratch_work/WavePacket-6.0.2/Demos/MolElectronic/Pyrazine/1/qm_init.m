% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init (~)
global atomic hamilt plots space time

log.disp ( '***************************************' )
log.disp ( 'Two state / three mode model of        ' )
log.disp ( 'S2->S1 internal conversion of pyrazine ' )
log.disp ( '***************************************' )

omega = [0.126 0.074 0.118]/atomic.E.eV; % Harmonic frequencies in eV

% Number of (coupled) Schr�dinger equations
hamilt.coupling.n_eqs     = 2;
hamilt.coupling.represent = 'dia';
hamilt.coupling.ini_rep   = 'dia';       
hamilt.coupling.ini_coeffs = [0 1];      % Initially upper adiabatic state populated

% Grid definition
space.dof{1}       = grid.fft;           % Using a FFT grid
space.dof{1}.x_min = -10;                % Lower bound of grid
space.dof{1}.x_max = 10;                 % Upper bound of grid
space.dof{1}.n_pts = 96;                 % Number of grid points
space.dof{1}.mass  = 1/omega(1);         % "Mass" for kinetic energy

space.dof{2}       = grid.fft;           % Using a FFT grid
space.dof{2}.x_min = -10;                % Lower bound of grid
space.dof{2}.x_max = 10;                 % Upper bound of grid
space.dof{2}.n_pts = 96;                 % Number of grid points
space.dof{2}.mass  = 1/omega(2);         % "Mass" for kinetic energy

space.dof{3}       = grid.fft;           % Using a FFT grid
space.dof{3}.x_min = -10;                % Lower bound of grid
space.dof{3}.x_max = 10;                 % Upper bound of grid
space.dof{3}.n_pts = 96;                 % Number of grid points
space.dof{3}.mass  = 1/omega(3);         % "Mass" for kinetic energy

% Temporal discretization
time.steps.m_start  = 000;               % Index of initial time step
time.steps.m_stop   = 300;               % Index of final time step
time.steps.m_delta  = 1/atomic.t.fs;     % Size of time steps: 1 femtosecond 
time.steps.s_number = 100;               % Number of sub steps per time step

% Propagator
time.propa = tmp.splitting;              % Split operator method
time.propa.order = 2;                    % Strang splitting

% Hamiltonian operator 
hamilt.truncate.e_min    =    0.0;       % Lower truncation of energy
hamilt.truncate.e_max    =   +1.0;       % Upper truncation of energy

hamilt.pot{1,1} = pot.pyrazine;          % Pyrazine model potential: S1
hamilt.pot{1,1}.omega  = omega;          % Force constants
hamilt.pot{1,1}.kappa = [+0.037 -0.105]/atomic.E.eV;  % Linear parameter (tuning)
hamilt.pot{1,1}.energ =   3.940/atomic.E.eV; % Vertical excitation energy

hamilt.pot{2,2} = pot.pyrazine;          % Pyrazine model potential: S2
hamilt.pot{2,2}.omega  = omega;          % Force constants
hamilt.pot{2,2}.kappa = [-0.254 +0.149]/atomic.E.eV;  % Linear parameter (tuning)
hamilt.pot{2,2}.energ =   4.840/atomic.E.eV; % Vertical excitation energy

hamilt.pot{1,2} = pot.pyrazine;          % Pyrazine model potential: S1-S2
hamilt.pot{1,2}.lambda =  +0.262/atomic.E.eV; % Interstate coupling

% Initial wave function
time.corr        = init.gauss_corr;      % Gaussian-shaped wavepacket
time.corr.width  = [1 1 1]/sqrt(2);      % Width of coherent state
time.corr.pos_0  = [0 0 0];              % Center in position representation
time.corr.mom_0  = [0 0 0];              % Center in momentum representation

% Modify settings for appearance of plots (if desired)
plots.density.type  = 'reduced';         % View point for surface plot: [az el]
plots.density.contour.nlev  = [60 15];   % Number of contours: density/energy

% Plot densities
plots.density          = vis.reduced_1d;
plots.density.cnt_nlev = 50;             % Number of contours: density

% Plot expectation values
plots.expect = vis.expect;