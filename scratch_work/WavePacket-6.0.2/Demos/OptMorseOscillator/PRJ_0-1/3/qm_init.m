% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global atomic control hamilt plots space time

log.disp ( '*****************************************************************' )
log.disp ( 'Vibrationally state selective excitation of a Morse ' )
log.disp ( 'oscillator resembling one of the OH bonds in water  ' )
log.disp ( 'Laser pulse optimized for 1-photon transition: 0->1 ' )
log.disp ( '                                                    ' )
log.disp ( 'see M.V.Korolkov, G.K.Paramonov, and B. Schmidt     ' )
log.disp ( 'J. Chem. Phys. 105(5), 1862-1879  (1996)            ' )
log.disp ( 'see W. Zhu, J. Botina, H. Rabitz                    ' )
log.disp ( 'J. Chem. Phys. 108(5), 1953-1963  (1998)            ' )
log.disp ( '*****************************************************************' )

% Isotopic masses
m_O = atomic.mass.O16;
m_H = atomic.mass.H1;

% Spatial discretization ==> qm_bound
space.dof{1}       = grid.fft;           % Using fft grid
space.dof{1}.mass  = m_O*m_H/(m_O+m_H);  % Reduced mass (OH radical)
space.dof{1}.n_pts = 512;                % Number of grid points
space.dof{1}.x_min =  0.7;               % Lower bound of grid 
space.dof{1}.x_max = 10.0;               % Upper bound of grid     

% Hamiltonian operator 
hamilt.truncate.e_min = -0.1;              % Lower truncation of energy
hamilt.truncate.e_max = +1.0;              % Upper truncation of energy

hamilt.pot{1,1}      = pot.morse;        % Harmonic oscillator
hamilt.pot{1,1}.d_e  = 0.1994;           % Dissociation energy
hamilt.pot{1,1}.r_e  = 1.821;            % Equilibrium length
hamilt.pot{1,1}.alf  = 1.189;            % Range parameter

hamilt.dip{1}{1,1}     = dip.mecke;      % Mecke dipole function
hamilt.dip{1}{1,1}.r_0 = 0.6/atomic.l.A; % Length parameter: 0.6 A
hamilt.dip{1}{1,1}.q_0 = 7.85/atomic.d.D*atomic.l.A; % Charge parameter: 7.85 D/A

hamilt.sbc{1,1}        = sbc.taylor;     % System-bath coupling: Taylor series
hamilt.sbc{1,1}.coeffs = 1;              % Linear coupling, slope 1

% Temporal discretization
time.steps.m_delta = 10/atomic.t.fs;        % Size of time steps: 10 fs 
time.steps.m_start = 0000;                  % Index of initial time step
time.steps.m_stop  = 0100;                  % Index of final time step
time.steps.s_number      = 200;                   % Number of sub steps per time step

time.frog.choice = 'none';               % Frequency resolved optical gating
time.frog.zoom = 15;                     % Zoom factor for frequency axes 

% Calculate and save (bound) eigen states (==> qm_bound)
hamilt.eigen.stop = 00;                     % Lower index
hamilt.eigen.stop = 21;                     % Upper index

% Electric field pulses: Initial guess
time.pulse{1}       = efi.sin_2;         % Shape of envelope
time.pulse{1}.delay = 0500/atomic.t.fs;  % Time delay of pulse center
time.pulse{1}.fwhm  = 0500/atomic.t.fs;  % Pulse length
time.pulse{1}.ampli = 0.005;             % Amplitude of electric field         

% Select ODE solver, parameters
control.solvers.handle1 = @ode.RuKu4;    % Runge-Kutta or other integrators
control.solvers.handle2 = @ode45;        % Runge-Kutta: Dormand-Prince (built-in)
control.solvers.reltol = 1e-6;           % Relative error tolerance of solution

% Initial state (==> qm_abncd)
control.initial.choice = 'pure';         % starting from a pure state
control.initial.pure = 0;                % vibrational ground state

% Define outputs: projections as observables
control.observe.types = 'prj';           % types of observables
control.observe.choices = {[0] [1] [2] [3] [4] [5:10] [11:21]};
control.observe.labels  = {'|0>' '|1>' '|2>' '|3>' '|4>' '|5:10>' '|11:21>'};
control.observe.targets = 1:7;           % choose control targets 

% Optimal control 
control.optimal.terminal = 2;            % Which observable to be optimized
control.optimal.tolerance = 1e-10;       % Threshold terminating iteration
control.optimal.max_iter = 050;          % Max. number of iteration steps
control.optimal.alpha = 1.00;            % Penalty factor for laser fluence
control.optimal.eta  = 1.00;             % Calculate optimal backward field
control.optimal.zeta = 1.00;             % Calculate optimal forward field
control.optimal.order = 2;               % Error order for optimal fields
control.optimal.prefactor = 'current';   % How to calculate prefactor (C2)
control.optimal.fb_test = false;         % Only testing forward/backward

% Plotting options
control.plot.uxy = true;                 % Plot u(t), x(t), y(t)
control.plot.j12 = true;                 % Plot j_1(t), j_2(t), and total
control.plot.psd = true;                 % Plot power spectral density
control.plot.mov = true;                 % Animation of u(t), x(t), y(t)

% Open quantum system: system-bath coupling, temperature
% These parameters yield 0.5 population of |1> after 1000 fs (LvNE)
control.lvne.temperature = 0.00;         % Temperature in atomic units: 315,923.5 K
control.lvne.order = 'df';               % Split into densities and coherences
control.relax.model = 'fermi';           % Relaxation model
control.relax.rate  = 1e-4;              % Relaxation rate (inverse of 1/e time)
control.relax.lower = 0;                 % Lower state for reference transition
control.relax.upper = 1;                 % Upper state for reference transition

% Modify settings for appearance of plots (if desired)
% plots.density.type        = 'contour';   % Contour plot of Wigner transform
% plots.density.range.on    = true;        % Ranges of contour plot
% plots.density.range.x_min = 0.7;
% plots.density.range.x_max = 10.0;
% plots.density.range.y_min = -25;
% plots.density.range.y_max = 25;

% Plot u(t), x(t), y(t)
plots.control = vis.uxy;
