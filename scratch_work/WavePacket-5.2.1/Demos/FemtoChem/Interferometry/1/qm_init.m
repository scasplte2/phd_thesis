% Copyright (C) 2009,2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots psi space time

util.disp ('**************************************')
util.disp ('Relaxation of an wavefunction in the  ')
util.disp ('electronic groundstate of Na2         ')
util.disp ('**************************************')

hamilt.coupling.labels = {'X^1\Sigma_g^+'};

% Grid definition
space.dof{1}       = grid_fft;  % using FFT grid
space.dof{1}.label = 'R';       % label for axis descriptions
space.dof{1}.n_pts = 256;       % number of points
space.dof{1}.x_min = 4;         % lower boundary of the grid
space.dof{1}.x_max = 14;        % upper boundary of the grid
space.dof{1}.mass  = 23 * 1836/2; % reduced mass of Na2

% Temporal discretisation
time.main.start = 0;            % index of first time step
time.main.stop  = 10;           % index of last time step

time.main.delta = 410;          % approx 10 fs per time step
time.sub.n      = 1;            % propagation time step not relevant (Chebychev)

% Propagator
time.propa.handle = @ket.cheby_imag;% Chebychev in imaginary time
time.propa.params.order = 0;        % auto-determine
time.propa.params.precision = 1e-10;% requested precision

% Hamiltonian operator
hamilt.truncate.min = -0.03;    % lower truncation of energy
hamilt.truncate.max = 0.05;     % upper truncation of energy

hamilt.pot.handle = @pot.interp;       % potential energy
hamilt.pot.params.pos_conv = 0.5292;   % conversion factor for coordinates
hamilt.pot.params.pot_conv = 27.212;   % conversion factor for energies
hamilt.pot.params.method   = 'spline'; % spline interpolation

% Initial wave function; Gaussian around the equilibrium
psi.init.dof{1}.handle = @wav.gauss;    % Gaussian initial wave function
psi.init.dof{1}.width  = 0.5;           % width of Gaussian
psi.init.dof{1}.pos_0  = 6;             % center in position representation
psi.init.dof{1}.mom_0  = 0;             % center in momentum representation

% Modify appearance of plots
plots.density.type        = 'contour';  % Draw contour lines

plots.density.range.on    = true;       % manually set ranges
plots.density.range.x_min = 4;
plots.density.range.x_max = 8;
plots.density.range.y_min = -10;
plots.density.range.y_max = 10;

plots.expect.energies.max = 4e-4;

plots.dvr.rho_max = 2;                  % manual height of wave function
