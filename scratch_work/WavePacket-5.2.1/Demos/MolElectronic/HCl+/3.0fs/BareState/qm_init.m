% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots psi space time

util.disp ( '*************************************' )
util.disp ( 'HCl+ cation ( X==>A excitation )' )
util.disp ( 'A.D.Pradhan, K.P.Kirby, A.Dalgarno  )' )
util.disp ( 'J. Chem. Phys. 95, 9009 (1991)      )' )
util.disp ( '*************************************' )

% Number of (coupled) Schr�dinger equations
hamilt.coupling.n_eqs = 2;
hamilt.coupling.representation='dia';
hamilt.coupling.labels = {'X ^2\Pi', 'A ^2\Sigma^+'};

% Spatial discretization
space.dof{1} = grid_fft;             % using FFT grid
space.dof{1}.mass = 0.97989*1822.9;  % Reduced mass
space.dof{1}.n_pts = 128;            % Number of grid points
space.dof{1}.x_min =  1.5;           % Lower bound of grid 
space.dof{1}.x_max =  7.5;           % Upper bound of grid

% Temporal discretization
time.main.start = 000;               % Index of initial time step
time.main.stop  = 100;               % Index of final time step

time.main.delta  = 0.03*41.341;      % Size of time steps: 30 as 
time.sub.n   =   030;                % Number of sub steps per time step

% Propagator
time.propa.handle = @ket.splitting;  % Split operator method
time.propa.params.order = 3;         % Strang splitting

% Electric field as sequence of pulses
time.efield.shape   = 'sin^2';       % Shape of envelope
time.efield.polar   = 0.0;           % Polarization angle [rad]
time.efield.delay   = 1.5*41.341;    % Time delay of pulse center: 1.5 fs
time.efield.fwhm    = 1.5*41.341;    % Full width at half maximum: 1.5 fs
time.efield.ampli   = 0.5;           % Field amplitude
time.efield.frequ   = 0.1295;        % Carrier frequency: 3.52 eV
time.efield.phase   = +time.efield.delay*time.efield.frequ+pi/2;    % Phase

% Hamiltonian operator 
hamilt.truncate.min    =  -1.9;      % Lower truncation of energy
hamilt.truncate.max    =  -1.2;      % Upper truncation of energy

hamilt.pot.handle      = @pot.interp;% Use tabulated values
hamilt.pot.params.pos_conv = 1;      % conversion factor for coordinates
hamilt.pot.params.pot_conv = 1;      % conversion factor for energies
hamilt.pot.params.method = 'spline'; % Spline interpolaation

hamilt.dip.handle      = @dip.interp;% Use tabulated values
hamilt.dip.params.pos_conv = 1;      % conversion factor for coordinates
hamilt.dip.params.dip_conv = 1;      % conversion factor for dipole moments
hamilt.dip.params.method = 'spline'; % Spline interpolaation

% Initial wave function
psi.init.dof{1}.handle = @wav.gauss; % Gaussian-shaped wavepacket
psi.init.dof{1}.width  =  0.153688;  % Width 
psi.init.dof{1}.pos_0  =  2.519868;  % Center in position representation
psi.init.dof{1}.mom_0  =  0.0;       % Center in momentum representation

psi.init.coeffs         = [1 0];     % Initially: ground state only
psi.init.representation = 'dia';

% Modify settings for appearance of plots (if desired)
plots.density.type = 'contour';

plots.expect.energies.min = -2;
plots.expect.energies.max = 0.5;