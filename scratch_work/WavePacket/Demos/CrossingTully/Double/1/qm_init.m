% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots psi space time

util.disp ( '***************************************************' )
util.disp ( 'Dual crossing example, k=13 (Tully 1990)' )
util.disp ( '***************************************************' )

% Number of (coupled) Schrödinger equations
hamilt.coupling.n_eqs = 2;
hamilt.coupling.representation = 'adi';

% Spatial discretization
space.dof{1} = grid_fft;             % using FFT grid
space.dof{1}.mass = 2000;            % mass
space.dof{1}.n_pts = 512;            % Number of grid points
space.dof{1}.x_min = -12.0;          % Lower bound of grid 
space.dof{1}.x_max =  12.0;          % Upper bound of grid

% Temporal discretization
time.main.start = 000;               % Index of initial time step
time.main.stop  = 040;               % Index of final time step

time.main.delta  = 050.0;            % Size of time steps 
time.sub.n       = 050;              % Number of sub steps per time step

% Propagator
time.propa.handle = @ket.splitting;  % Split operator method
time.propa.params.order = 3;         % Strang splitting

% Hamiltonian operator 
hamilt.truncate.min    = -0.05;      % Lower truncation of energy
hamilt.truncate.max    =  0.30;      % Upper truncation of energy

hamilt.pot.handle      = @pot.tully2;% Single Crossing
hamilt.pot.params.A    = 0.10;       % Tully's parameter A
hamilt.pot.params.B    = 0.28;       % Tully's parameter B
hamilt.pot.params.C    = 0.015;      % Tully's parameter C
hamilt.pot.params.D    = 0.06;       % Tully's parameter D
hamilt.pot.params.E    = 0.05;       % Tully's parameter E

% Absorbing boundary conditions
hamilt.nip.handle      = @nip.power; % Negative imaginary potential
hamilt.nip.params.exp  = 4;          % Exponent
hamilt.nip.params.min = -10;         % Beginning of inner grid region
hamilt.nip.params.max = +10;         % End of inner grid region

% Initial wave function
psi.init.dof{1}.handle= @wav.gauss;  % Gaussian-shaped wavepacket
psi.init.dof{1}.width =  0.75;       % Width 
psi.init.dof{1}.pos_0 = -7.0;        % Center in position representation
psi.init.dof{1}.mom_0 = 13.0;        % Center in momentum representation

psi.init.representation = 'adi';
psi.init.coeffs       = [1 0];       % Initially only lower adiabatic state populated

% Modify settings for appearance of plots (if desired)
plots.density.type = 'contour';

plots.density.surface.view  = [80 06];   % View point for surface plot: [az el]

if strcmp(plots.density.type, 'contour')
    plots.density.range.on = true;       % manual setting of plotting range
    plots.density.range.x_min = -9;
    plots.density.range.x_max = 9;
    plots.density.range.y_min = 0;
    plots.density.range.y_max = 30;
end

plots.expect.energies.max = 0.1;
