% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008 Ulf Lorenz
%               2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init (~)
global hamilt plots space wave

log.disp ( '*****************************' )
log.disp ( 'Strontium 86 Ground State' )
log.disp ( '*****************************' )

% Spatial discretization
space.dof{1} = grid.fft;                 % using fft grid
space.dof{1}.mass = 7.830152319224479e+04;            % Reduced mass
space.dof{1}.n_pts = 512;                % Number of grid points
space.dof{1}.x_min = 5.0;                % Lower bound of grid 
space.dof{1}.x_max = 3000.0;               % Upper bound of grid
space.dof{1}.periodic = false;           % Build the kinetic energy matrix
                                         % without periodic boundary conditions
% Hamiltonian operator 
hamilt.truncate.e_min  =  -1e-2;           % Lower truncation of energy
hamilt.truncate.e_max  =  1e-2;           % Upper truncation of energy

hamilt.pot{1,1}      = pot.interp;        % From Tiemann calculation

% Select eigen/values/functions
hamilt.eigen.start     = 0;             % Lower index
hamilt.eigen.stop      = 1;             % Upper index

% Plot time evolution of the density
plots.density         = vis.contour;     % Contour plot of Wigner function
plots.density.expect  = false;

% Plot time evolution of expectation values
plots.expect          = vis.expect;
plots.expect.e_max    = 0.25;            % Manually set range for energy plot

wave.sav_export = 'true';
