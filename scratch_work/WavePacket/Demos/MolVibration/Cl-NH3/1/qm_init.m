% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots psi space time

util.disp ( '*********************************************' )
util.disp ( 'Tabulated potential (neutral state of Cl-NH3)' )
util.disp ( 'Physical Review Letters 93 (4), 048301 (2004)' )
util.disp ( '*********************************************' )

% Spatial discretization
space.dof{1}       = grid_fft;           % Fourier grid
space.dof{1}.label = 'NH_2-H';           % Axis label
space.dof{1}.x_min = 0.25/0.5292;        % Lower bound of grid
space.dof{1}.x_max = 3/0.5292;           % Upper boundof grid
space.dof{1}.n_pts = 128;                % Number of grid points

space.dof{2}       = grid_fft;           % Fourier grid
space.dof{2}.label = 'H-Cl';             % Axis label
space.dof{2}.x_min = 0.70/0.5292;        % Lower bound of grid
space.dof{2}.x_max = 3.1/0.5292;         % Upper boundof grid
space.dof{2}.n_pts = 128;                % Number of grid points

% Temporal discretization
time.main.start = 000;                   % Index of initial time step
time.main.stop  = 020;                   % Index of final time step

time.main.delta = 050;                   % Size of time steps 
time.sub.n      = 100;                   % Number of sub steps per time step

% Propagator
time.propa.handle = @ket.splitting;      % Split operator method
time.propa.params.order = 3;             % Strang splitting

% Hamiltonian operator 
hamilt.truncate.min    =  -0.05;         % Lower truncation of energy
hamilt.truncate.max    =  +0.20;         % Upper truncation of energy

hamilt.kin{1}      = kin_triatomic;      % Jacobi coordinates with fixed angle
hamilt.kin{1}.dof  = [1 2];              % on which degrees of freedom we act
hamilt.kin{1}.mass = [16 1 35] * 1822.9; % masses of NH2,H,Cl
hamilt.kin{1}.theta= pi;                 % NH2-H-Cl bending angle

hamilt.pot.handle      = @pot.interp;    % Interpolation of potential fct.
hamilt.pot.params.n_pts  = [56 49];      % Number of tabulated values
hamilt.pot.params.pos_conv = 0.5292;     % Conversion of coordinates
hamilt.pot.params.pot_conv = 219473.6;   % Conversion of pot. energies
hamilt.pot.params.method  = 'spline';    % Spline method

% Initial wave function
psi.init.dof{1}.handle = @wav.gauss;
psi.init.dof{1}.width  = 0.28 / 2;
psi.init.dof{1}.pos_0  = 1.95;
psi.init.dof{1}.mom_0  = 0;

psi.init.dof{2}.handle = @wav.gauss;
psi.init.dof{2}.width  = 0.47 / 2;
psi.init.dof{2}.pos_0  = 4.52;
psi.init.dof{2}.mom_0  = 0;

% Modify settings for appearance of plots (if desired)
plots.density.surface.view  = [60 65];   % View point for surface plot: [az el]
plots.density.type = 'contour';          % Choose plot type
plots.expect.energies.max = 0.015;       % Manually set the range for the energy plot
