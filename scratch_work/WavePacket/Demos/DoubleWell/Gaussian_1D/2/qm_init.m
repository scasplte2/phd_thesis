% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots psi space time

util.disp ( '*******************************************' )
util.disp ( 'Double well: Initial energy above threshold' )
util.disp ( '*******************************************' )

% Spatial discretization
space.dof{1} = grid_fft;                 % using FFT grid
space.dof{1}.mass = 1;                   % mass
space.dof{1}.n_pts = 128;                % Number of grid points
space.dof{1}.x_min = -8.2;               % Lower bound of grid 
space.dof{1}.x_max =  8.2;               % Upper bound of grid

% Temporal discretization and propagator
time.main.start = 000;                   % Index of initial time step
time.main.stop  = 025;                   % Index of final time step

time.main.delta  = 0.2;                  % Size of time steps 
time.sub.n   =  100;                     % Number of sub steps per time step

% Propagator
time.propa.handle = @ket.splitting;      % Split operator method
time.propa.params.order = 3;             % Strang splitting

% Hamiltonian operator 
hamilt.truncate.min    = -15.0;          % Lower truncation of energy
hamilt.truncate.max    =  85.0;          % Upper truncation of energy

hamilt.pot.handle      = @pot.taylor;    % Taylor series: Double well potential
hamilt.pot.params.v{1,1} = [0;-3;0;+1];  % Quadratic and quartic constant

% Initial wave function
psi.init.dof{1}.handle= @wav.gauss;      % Gaussian-shaped wavepacket
psi.init.dof{1}.width = sqrt(1/2);       % Width 
psi.init.dof{1}.pos_0 = -4.0;            % Center in position representation
psi.init.dof{1}.mom_0 =  6.0;            % Center in momentum representation

% Modify settings for appearance of plots (if desired)
plots.density.type = 'contour';

plots.expect.energies.max = 20;
