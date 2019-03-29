% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%               2011 Ulf Lorenz
%                    Boris Schaefer-Bung
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init (state)
global hamilt plots space time control reduce

log.disp ( '**********************************' )
log.disp ( 'Asymmetric double well potential  ' )
log.disp ( ' see example 1 in:                ' )
log.disp ( ' B. Schaefer-Bung, C. Hartmann,   ' )
log.disp ( ' B. Schmidt, and Ch. Schuette,    ' )
log.disp ( ' J. Chem. Phys. 135, 014112 (2011)' )
log.disp ( '**********************************' )

% Spatial discretization
space.dof{1} = grid.fft;                 % using FFT grid
space.dof{1}.mass = 162;                 % I=162 [hbar^2/D], see section IV.B
space.dof{1}.n_pts = 256;                % Number of grid points
space.dof{1}.x_min = -2.0;               % Lower bound of grid 
space.dof{1}.x_max =  2.0;               % Upper bound of grid

% Temporal discretization
time.steps.m_start = 0;                  % Index of initial time step
time.steps.m_stop  = 100;                % Index of final time step
time.steps.m_delta = 3*pi/(100*0.196);   % Size of time steps 

% Propagator
time.propa = tmp.splitting;              % Operator splitting
time.propa.order = 2;                    % Strang splitting

% Electric field as half-cycle pulse
time.pulse{1}       = efi.recta;         % Rectangular shape of envelope
time.pulse{1}.delay = pi/(2*0.196);      % Time delay of pulse center
time.pulse{1}.fwhm  = pi/0.196;          % Full width at half maximum
time.pulse{1}.ampli = 7.5;               % u_0=7.5 [D/(e*a_0)] see section IV.B
time.pulse{1}.frequ = 0.196;             % Omega = 0.196 [D/hbar], see section IV.B

% Hamiltonian operator 
hamilt.truncate.e_min  =  -1.0;          % Lower truncation of energy
hamilt.truncate.e_max  =   10.0;         % Upper truncation of energy

hamilt.pot{1,1}        = pot.taylor;     % Taylor series: Double well potential
hamilt.pot{1,1}.coeffs = [0.055;-4;0;+24];  % Linear, quadratic, quartic constant
hamilt.pot{1,1}.vshift = 1;              % Constant energy shift

hamilt.dip{1}{1,1}        = dip.taylor;  % Dipole moment: Taylor series
hamilt.dip{1}{1,1}.coeffs = 1;           % Linear dipole moment, slope 1

hamilt.sbc{1,1}        = sbc.taylor;     % System-bath coupling: Taylor series
hamilt.sbc{1,1}.coeffs = 1;              % Linear coupling, slope 1

% Calculate (bound) eigen states (==> qm_bound)
hamilt.eigen.start = 00;                 % Lower index
hamilt.eigen.stop  = 20;                 % Upper index

% Save (bound) eigen states (==> qm_bound)
state.sav_export = true;                 % Toggle saving
state.sav_file   = 'bound';              % Filename for saving
state.sav_dir    = pwd;                  % Directory for saving

% Select ODE solver, parameters (==> qm_control)
control.solvers.handle2 = @ode45;        % Runge-Kutta
control.solvers.reltol = 1e-6;           % Relative error tolerance for ode45
          
% Initial density
control.initial.choice = 'thermal';      % thermal = Boltzmann density
control.initial.temperature = 1.0;       % choice of temperature

% Temperature and system-bath coupling
control.lvne.order = 'df';               % ordering vectorized density matrices: diagonals first
control.lvne.temperature = 1.0;          % temperature in atomic units: 315,923.5 K
control.relax.rate =  0.264386355;       % Gamma_{2->0} should be = 1 [D/hbar]
                                         % set to 0.2644 due to error in ancestor code 
control.relax.lower = 0;                 % Lower state for reference transition
control.relax.upper = 2;                 % Upper state for reference transition
control.relax.model = 'fermi';           % Omega dependent (Andrianov&Saalfrank)

% Define output observables and choose control targets
control.observe.types = 'prj';           % choose populations of observables
control.observe.choices = {0:2:10 1:2:9 11:20}; 
control.observe.labels = {'left well' 'right well' 'delocalized'};
control.observe.targets = 1:3;           % choose control targets 

% Parameters for balancing transformation
reduce.balance.A_stable = 'ssu';         % SSU or EVS method for stabilizing A  
reduce.balance.A_shift = 1e-6;           % magnitude of eigenvalue shift
reduce.balance.A_split = 1;              % dimension of instable part of A
reduce.balance.BN_scale = 3;             % Scaling factor for control fields
reduce.balance.acf_couple = true;        % additional control field coupling A and rho/x
reduce.balance.method = 'iter';          % ITER or BICG solver for GLYAP
reduce.balance.transform = 'srbt';       % SRBT or MRMR balancing transform;

% Parameters for H2 optimal model reduction (Benner/Breiten @ MPI Magdeburg)
reduce.H2model.A_stable = 'ssu';         % SSU or EVS method for stabilizing A  
reduce.H2model.A_shift = 1e-6;           % magnitude of eigenvalue shift
reduce.H2model.A_split = 1;              % dimension of instable part of A
reduce.H2model.BN_scale = 3;             % Scaling factor for control fields
reduce.H2model.acf_couple = true;        % additional control field coupling A and rho/x
reduce.H2model.method = 'iter';          % ITER or BICG solver for BIRKA

% Parameters for calculation of H2 error
reduce.H2error.A_shift = 1e-2;           % magnitude of eigenvalue shift
reduce.H2error.BN_scale = 10;            % Scaling factor for control fields
reduce.H2error.method = 'iter';          % ITER or BICG solver for GSYLV

% Plot of densities
plots.density        = vis.contour;
plots.density.expect = false;

% Plot of expectation values as time series
plots.expect       = vis.expect;        
plots.expect.p_min = 0.26;               % Range of population plot
plots.expect.p_max = 0.42;               % Range of population plot
plots.expect.e_max = 2;                  % Range of energy plot

% Plot of input/output
plots.control      = vis.uxy;


