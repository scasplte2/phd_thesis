% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008 Ulf Lorenz
%               2008 Burkhard Schmidt
%               2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots psi space

util.disp ( '*********************************************' )
util.disp ( 'Fulvene torsion: Stationary torsion (A state)' )
util.disp ( '*********************************************' )

% Spatial discretization
space.dof{1} = grid_fft;                 % using fft grid
space.dof{1}.mass = 11152;               % Reduced moment of inertia: 1/2.44 meV
space.dof{1}.n_pts = 128;                % Number of grid points
space.dof{1}.x_min = -pi;                % Lower bound of grid 
space.dof{1}.x_max =  pi;                % Upper bound of grid

% Hamiltonian operator 
hamilt.eigen.symmetry  = 'u';            % Odd parity eigenstates

hamilt.truncate.min    =  0.0;           % Lower truncation of energy
hamilt.truncate.max    = 0.25;           % Upper truncation of energy

hamilt.pot.handle    = @pot.pendulum;    % Intramolecular torsion
hamilt.pot.params.zeta = - 0.076757;     % Barrier height: 856*2.44 meV

% Select eigen/values/functions
psi.eigen.start        = 0;              % Lower index
psi.eigen.stop         =10;              % Upper index

% Modify settings for appearance of plots (if desired)
plots.density.type        = 'contour';

plots.expect.energies.max = 0.05;