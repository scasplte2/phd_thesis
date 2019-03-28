% Copyright (C) 2008-2009 Ulf Lorenz
%               2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots psi space

util.disp ( '**********************************************' )
util.disp ( 'Symmetric double well potential               ' )
util.disp ( 'Reference calculation using Gauss-Hermite DVR ' )
util.disp ( '**********************************************' )

% Spatial discretization
space.dof{1}       = grid_hermite;       % Use Gauss-Hermite grid
space.dof{1}.mass  = 1;                  % mass
space.dof{1}.v_2   = 3;                  % frequency of oscillator
space.dof{1}.n_pts = 64;                 % Number of grid points

% Hamiltonian operator 
% hamilt.eigen.symmetry = 'u';             % Symmetry may be exploited

hamilt.truncate.min    = -15.0;          % Lower truncation of energy
hamilt.truncate.max    =  85.0;          % Upper truncation of energy

hamilt.pot.handle      = @pot.taylor;    % Taylor series: Double well potential
hamilt.pot.params.v{1,1} = [0;-3;0;+1];  % Quadratic and quartic constant

% Select eigen/values/functions
hamilt.eigen.symmetry = 'g';             % Even parity only
psi.eigen.start        = 0;              % Lower index
psi.eigen.stop         = 15;             % Upper index

% Modify settings for appearance of plots (if desired)
plots.density.type = 'curve';            % simple plot of wavefunction

plots.density.surface.view  = [60 75];   % View point for surface plot: [az el]

plots.expect.energies.max = 20;
