% Copyright (C) 2008-2009 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots psi space time

util.disp ( '*******************************************' )
util.disp ( 'Double well: Initial energy above threshold' )
util.disp ( 'Calculated using Gauss-Hermite quadrature. ' )
util.disp ( '*******************************************' )

% Spatial discretization
space.dof{1}       = grid_hermite;       % use Gauss-Hermite grid
space.dof{1}.mass  = 1;                  % mass
space.dof{1}.v_2   = 3;                  % frequency of the oscillator
space.dof{1}.n_pts = 40;                 % number of basis functions

% Temporal discretization and propagator
time.main.start = 000;                   % Index of initial time step
time.main.stop  = 025;                   % Index of final time step

time.main.delta  = 0.2;                  % Size of time steps 
time.sub.n   =  100;                     % Number of sub steps per time step

% Propagator
time.propa.handle = @ket.splitting;      % Split operator method
time.propa.params.order = 3;             % Strang splitting

% Hamiltonian operator 
hamilt.truncate.min    = -15.0;          % Lower truncation of energy
hamilt.truncate.max    =  85.0;          % Upper truncation of energy

hamilt.pot.handle      = @pot.taylor;    % Taylor series: Double well potential
hamilt.pot.params.v{1,1} = [0;-3;0;+1];  % Quadratic and quartic constant

% Initial wave function
psi.init.dof{1}.handle= @wav.gauss;      % Gaussian-shaped wavepacket
psi.init.dof{1}.width = sqrt(1/2);       % Width 
psi.init.dof{1}.pos_0 = -4.0;            % Center in position representation
psi.init.dof{1}.mom_0 =    0;            % Center in momentum representation

% Modify settings for appearance of plots (if desired)
plots.density.type        = 'curve';     % simple plot of wavefunction

plots.expect.energies.max = 10;

plots.dvr.rho_max = 1;                   % without manual intervention, the plot
                                         % looks bad
