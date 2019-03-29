% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%               2011 Ulf Lorenz
%                    Boris Schaefer-Bung
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global control hamilt plots psi space time dim_red

util.disp ( '**********************************' )
util.disp ( 'Asymmetric double well potential  ' )
util.disp ( ' see example 3 in:                ' )
util.disp ( ' B. Schaefer-Bung, C. Hartmann,   ' )
util.disp ( ' B. Schmidt, and Ch. Schuette,    ' )
util.disp ( ' J. Chem. Phys. 135, 014112 (2011)' )
util.disp ( '**********************************' )

% Spatial discretization
space.dof{1} = grid_fft;                 % using FFT grid
space.dof{1}.mass = 162;                 % I=162 [hbar^2/D], see section IV.D
space.dof{1}.n_pts = 256;                % Number of grid points
space.dof{1}.x_min = -2.0;               % Lower bound of grid 
space.dof{1}.x_max =  2.0;               % Upper bound of grid

% Temporal discretization
time.main.start = 0;                     % Index of initial time step
time.main.stop  = 100;                   % Index of final time step
time.main.delta = 7.2;                   % Size of time steps 

% Propagator
time.propa.handle = @ket.splitting;      % Operator splitting
time.propa.params.order = 3;             % Symmetrized (Strang) scheme

% Electric field as half-cycle pulse
time.efield.dressed = false;             % Not using dressed state picture
time.efield.shape   = 'sin^2';           % sin^2 shape of envelope
time.efield.polar   = 0;                 % Polarization angle [rad]
time.efield.delay   = 180;               % Time delay of pulse center
time.efield.fwhm    = 180;               % Full width at half maximum
time.efield.ampli   = 0.2;               % u_0=0.2 [D/(e*a_0)] see section IV.D
time.efield.frequ   = 0.196;             % Omega = 0.196 [D/hbar], see section IV.D
time.efield.phase   = 35.28;             % Phase = frequ * delay = 0.196* 180

% Hamiltonian operator 
hamilt.truncate.min    =  -1.0;          % Lower truncation of energy
hamilt.truncate.max    =   10.0;         % Upper truncation of energy

hamilt.pot.handle = @pot.taylor;         % Taylor series: Double well potential
hamilt.pot.params.v{1,1} = [0.055;-4;0;+24]; % Linear, quadratic, quartic constant
hamilt.pot.params.c{1,1} = 1;            % Constant energy shift

hamilt.dip.handle = @dip.taylor;         % Dipole moment: Taylor series
hamilt.dip.params.x.v{1,1} = 1;          % Linear dipole moment, slope 1

hamilt.sbc.handle = @sbc.taylor;         % System-bath coupling: Taylor series
hamilt.sbc.params.v{1,1} = 1;            % Linear coupling, slope 1

% Calculate and save (bound) eigen states (==> qm_bound)
psi.eigen.start        = 0;              % Lower index
psi.eigen.stop         = 20;             % Upper index
psi.save.export        = true;

% Select ODE solver, parameters (==> qm_control)
control.solvers.handle2 = @ode45;        % Runge-Kutta
control.solvers.reltol = 1e-6;           % Relative error tolerance for ode45

% Initial density
control.initial.choice = 'thermal';      % thermal = Boltzmann density
control.initial.temperature = 0.1;       % choice of temperature

% Temperature and system-bath coupling
control.lvne.order = 'df';               % ordering vectorized density matrices: diagonals first
control.lvne.temperature = 0.1;          % temperature in atomic units: 315,923.5 K
control.lvne.relax.rate =  3.73779924e-4;% Gamma_{2->0} should be = 1e-3 [D/hbar] 
control.lvne.relax.lower = 0;            % Lower state for reference transition
control.lvne.relax.upper = 2;            % Upper state for reference transition
control.lvne.relax.model = 'omega';      % Omega dependent (Andrianov&Saalfrank)

% Define output observables and choose control targets
control.observe.choices={0:2:10 1:2:9 11:20}; 
control.observe.labels={'left well' 'right well' 'delocalized'};
control.observe.targets=1:3; 

% Parameters for balanced truncation
dim_red.balance.A_stable = 'ssu';        % SSU or EVS method for stabilizing A  
dim_red.balance.A_shift = 1e-6;          % magnitude of eigenvalue shift
dim_red.balance.A_split = 1;             % dimension of instable part of A
dim_red.balance.BN_scale = 250;          % Scaling factor for control fields
dim_red.balance.acf_couple = true;       % additional control field coupling A and rho/x
dim_red.balance.method = 'iter';         % ITER or BICG solver for GLYAP
dim_red.balance.transform = 'srbt';      % SRBT or MRMR balancing transform;
dim_red.balance.truncate = 'normal';     % normal or confinement truncation

% Parameters for H2 model reduction (Benner/Breiten @ MPI Magdeburg)
dim_red.reduce.A_stable = 'ssu';         % SSU or EVS method for stabilizing A  
dim_red.reduce.A_shift = 1e-6;           % magnitude of eigenvalue shift
dim_red.reduce.A_split = 1;              % dimension of instable part of A
dim_red.reduce.BN_scale = 150;           % Scaling factor for control fields
dim_red.reduce.acf_couple = true;        % additional control field coupling A and rho/x
dim_red.reduce.method = 'iter';          % ITER or BICG solver for BIRKA

% Parameters for calculation of H2 error
dim_red.error.A_shift = 1e-2;            % magnitude of eigenvalue shift
dim_red.error.BN_scale = 10;             % Scaling factor for control fields
dim_red.error.method = 'iter';           % ITER or BICG solver for GSYLV

% Modify settings for appearance of plots (if desired)
plots.density.type = 'contour';          % Contour plot of Wigner transform
plots.density.export.on = true;
plots.expect.export.on = true;
% plots.expect.population.min = -0.101;
% plots.expect.population.max =  0.062;

plots.expect.energies.max = 2;           % Range of energy plot
