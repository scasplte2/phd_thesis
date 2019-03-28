% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots psi space time

util.disp ( '*******************************************************' )
util.disp ( 'Free particle with zero momentum: Wavepacket dispersion' )
util.disp ( '*******************************************************' )

% Spatial discretization
space.dof{1} = grid_fft;                 % using fft grid
space.dof{1}.mass = 1;                   % mass
space.dof{1}.n_pts = 128;                % Number of grid points
space.dof{1}.x_min = -20.0;              % Lower bound of grid 
space.dof{1}.x_max = +20.0;              % Upper bound of grid

% Temporal discretization and propagator
time.main.start = 000;                   % Index of initial time step
time.main.stop  = 010;                   % Index of final time step

time.main.delta  = 0.5;                  % Size of time steps 
time.sub.n       = 050;                  % Number of sub steps per time step

time.propa.handle = @ket.splitting;      % Split operator method
time.propa.params.order = 3;             % Strang splitting

% Hamiltonian operator 
hamilt.truncate.min    =   0.0;          % Lower truncation of energy
hamilt.truncate.max    =  10.0;          % Upper truncation of energy

% Absorbing boundary conditions
hamilt.nip.handle      = @nip.power;     % Negative imaginary potential
hamilt.nip.params.exp  = 2;              % Exponent
hamilt.nip.params.min = -12;             % Beginning of inner grid region
hamilt.nip.params.max =  12;             % End of inner grid region

% Initial wave function
psi.init.dof{1}.handle= @wav.gauss; % Gaussian-shaped wavepacket
psi.init.dof{1}.width = sqrt(1/2);       % Width 
psi.init.dof{1}.pos_0 =  0.0;            % Center in position representation
psi.init.dof{1}.mom_0 =  0.0;            % Center in momentum representation

% Modify settings for appearance of plots (if desired)
plots.density.type          = 'contour'; % contour plot

plots.density.range.on    = true;        % manual setting of plotting range
plots.density.range.x_min = -10;
plots.density.range.x_max = 10;
plots.density.range.y_min = -3;
plots.density.range.y_max = 3;
plots.expect.energies.max = 0.5;
