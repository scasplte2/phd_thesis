% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%               2008-2009 Burkhard Schmidt
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots psi space time

util.disp ( '************************************************' )
util.disp ( 'Fulvene torsionial dynamics (electronic B state)' )
util.disp ( '************************************************' )

% Spatial discretization
space.dof{1} = grid_fft;                 % using fft grid
space.dof{1}.mass = 11152;               % Reduced moment of inertia: 1/2.44 meV
space.dof{1}.n_pts = 128;                % Number of grid points
space.dof{1}.x_min = -pi;                % Lower bound of grid 
space.dof{1}.x_max =  pi;                % Upper bound of grid

% Temporal discretization
time.main.start = 000;                   % Index of initial time step
time.main.stop  = 075;                   % Index of final time step

time.main.delta  = 206.7;                % Size of time steps: 05 fs
time.sub.n   =  50;                      % Number of sub steps per time step

% Propagator
time.propa.handle = @ket.splitting;      % Split operator method
time.propa.params.order = 3;             % Strang splitting

% Hamiltonian operator 
hamilt.truncate.min    = -0.015;          % Lower truncation of energy
hamilt.truncate.max    = +0.050;          % Upper truncation of energy

hamilt.pot.handle  = @pot.pendulum;       % Intramolecular torsion
hamilt.pot.params.zeta = + 0.0095;        % Barrier height: 106*2.44 meV

% Initial wave function
psi.init.dof{1}.handle   = @wav.pendulum;% Eigenstate of plane pendulum
psi.init.dof{1}.parity   = 'c';          % Parity: cosine=even
psi.init.dof{1}.order    = 0;            % Order (quantum number)
psi.init.dof{1}.multiple = 2;            % Multiplicity: Double well
psi.init.dof{1}.barrier  = 0.076757*space.dof{1}.mass; % 856*2.44 meV
psi.init.dof{1}.shift    = 0;            % Potential shift: None

% Modify settings for appearance of plots (if desired)
plots.density.type          = 'contour'; % Contour plot of the Wigner function

plots.pot.min = -10e-3;                  % adjust plot energy range
plots.pot.max = 4e-3;