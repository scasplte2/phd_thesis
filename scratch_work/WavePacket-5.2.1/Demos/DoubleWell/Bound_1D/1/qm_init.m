% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%               2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots psi space

util.disp ( '***********************************' )
util.disp ( 'Symmetric double well potential    ' )
util.disp ( 'Reference calculation using FFT-DVR' )
util.disp ( '***********************************' )

% Spatial discretization
space.dof{1} = grid_fft;                 % using FFT grid
space.dof{1}.mass = 1;                   % mass
space.dof{1}.n_pts = 128;                % Number of grid points
space.dof{1}.x_min = -8.2;               % Lower bound of grid 
space.dof{1}.x_max =  8.2;               % Upper bound of grid

% Hamiltonian operator 
hamilt.truncate.min    = -15.0;          % Lower truncation of energy
hamilt.truncate.max    =  85.0;          % Upper truncation of energy

hamilt.pot.handle      = @pot.taylor;    % Taylor series: Double well potential
hamilt.pot.params.v{1,1} = [0;-3;0;+1];  % Quadratic and quartic constant

% Select eigen/values/functions
hamilt.eigen.symmetry = 'g';             % Even parity only
psi.eigen.start        = 0;              % Lower index
psi.eigen.stop         = 15;             % Upper index

% Modify settings for appearance of plots (if desired)
plots.density.type = 'contour';

plots.expect.energies.max = 20;
