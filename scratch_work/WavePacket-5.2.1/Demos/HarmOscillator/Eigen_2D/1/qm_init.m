% Copyright (C) 2008-2009 Ulf Lorenz
%               2008 Burkhard Schmidt
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots psi space time

util.disp ( '********************************************' )
util.disp ( 'Eigenstates of the harmonic oscillator      ' )
util.disp ( 'We propagate an eigenstate of a 2D harmonic ' )
util.disp ( 'oscillator in time and check whether the    ' )
util.disp ( 'autocorrelation looks good.                 ' )
util.disp ( '********************************************' )

% Spatial discretization
space.dof{1}       = grid_fft;           % Using FFT grid 
space.dof{1}.n_pts = 32;                 % Number of grid points
space.dof{1}.x_min = -5.0;               % Lower bound of grid 
space.dof{1}.x_max = +5.0;               % Upper bound of grid
space.dof{1}.mass  = 1.0;                % Mass

space.dof{2}       = grid_fft;           % Using FFT grid 
space.dof{2}.n_pts = 64;                 % Number of grid points
space.dof{2}.x_min = -8.0;               % Lower bound of grid 
space.dof{2}.x_max = +8.0;               % Upper bound of grid
space.dof{2}.mass  = 1.0;                % Mass

% Temporal discretization
time.main.start = 000;                   % Index of initial time step
time.main.stop  = 30;                    % Index of final time step

time.main.delta = 1;                     % Size of time steps 
time.sub.n      = 0100;                  % Number of sub steps per time step

% Propagator
time.propa.handle = @ket.splitting;      % Split operator method
time.propa.params.order = 3;             % Strang splitting

% Hamiltonian
hamilt.truncate.min    =   0.0;          % Lower truncation of energy
hamilt.truncate.max    = 100.0;          % Upper truncation of energy

hamilt.pot.handle = @pot.taylor;         % Taylor series
hamilt.pot.params.v{1,1} =  [0 0;1 0.25];% Force constants

% Initial wave function
psi.init.dof{1}.handle   = @wav.harmonic;
psi.init.dof{1}.omega    = 1;
psi.init.dof{1}.r_e      = 0;
psi.init.dof{1}.n_q      = 0;

psi.init.dof{2}.handle   = @wav.harmonic;
psi.init.dof{2}.omega    = 0.5;
psi.init.dof{2}.r_e      = 0;
psi.init.dof{2}.n_q      = 3;

% Plots do not really tell us anything interesting here, so turn
% them off for speed.
plots.density.on = false;
plots.expect.on  = false;
