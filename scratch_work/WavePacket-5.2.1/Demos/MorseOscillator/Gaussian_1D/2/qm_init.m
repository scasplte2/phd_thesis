% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots psi space time

util.disp ( '*************************************' )
util.disp ( 'Gaussian packet in a Morse oscillator' )
util.disp ( '*************************************' )

% Spatial discretization
space.dof{1} = grid_fft;                 % using fft grid
space.dof{1}.mass = 1728.539;            % Reduced mass (OH molecule)
space.dof{1}.n_pts = 128;                % Number of grid points
space.dof{1}.x_min = 0.7;                % Lower bound of grid 
space.dof{1}.x_max = 7.0;                % Upper bound of grid

% Temporal discretization
time.main.start =  000;                  % Index of initial time step
time.main.stop  = 200;                  % Index of final time step

time.main.delta  = 50;                   % Size of time steps 
time.sub.n   =    100;                   % Number of sub steps per time step

% Propagator
time.propa.handle = @ket.cheby_real;     % Chebychev polynomials
time.propa.params.order = 0;             % Automatically detecting order
time.propa.params.precision = 10^-8;     % Threshold for truncation

% Hamiltonian operator 
hamilt.truncate.min    =  0.0;           % Lower truncation of energy
hamilt.truncate.max    =  1.0;           % Upper truncation of energy

hamilt.pot.handle      = @pot.morse;     % Harmonic oscillator
hamilt.pot.params.d_e  = 0.1994;         % Dissociation energy
hamilt.pot.params.r_e  = 1.821;          % Equilibrium length
hamilt.pot.params.alf  = 1.189;          % Range parameter
hamilt.pot.params.t_e  = 0.0;            % Energetic shift

% Absorbing boundary conditions
hamilt.nip.handle = @nip.power;          % Negativ imaginary potential
hamilt.nip.params.exp  = 4;              % Exponent
hamilt.nip.params.min  = 1.0;            % Beginning of non-gobbler region
hamilt.nip.params.max  = 6.0;            % End of non-gobbler region

% Initial wave function
psi.init.dof{1}.handle= @wav.gauss; % Gaussian-shaped wavepacket
psi.init.dof{1}.width =  0.1;            % Width 
psi.init.dof{1}.pos_0 =  1.44;           % Center in position representation
psi.init.dof{1}.mom_0 =  0.0;            % Center in momentum representation

% Modify settings for appearance of plots (if desired)
plots.density.type        = 'contour';

plots.expect.energies.max = 0.1;
