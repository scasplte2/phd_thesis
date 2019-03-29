% Copyright (C) 2009,2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init(symmetry)
global hamilt plots psi space

util.disp ( '*******************************************************' )
util.disp ( 'Calculation of the vibrational spectrum of the S0 state' )
util.disp ( 'of 9-(N-carbazolyl)-anthracene (C9A)                   ' )
util.disp ( 'see Z.Phys.D 34:111                                    ' )
util.disp ( '*******************************************************' )

% moment of intertia is hbar/(4pi*c * 0.035 cm^-1). Express these
% with atomic units
inertia = 1/(4 * pi * 137.12 * 0.035 * 0.529e-8);

% Spatial discretization
space.dof{1}       = grid_fft;           % Fourier grid
space.dof{1}.n_pts = 128;                % Number of grid points
space.dof{1}.x_min = 0.524;              % Lower bound of grid 30 degrees
space.dof{1}.x_max = pi-0.524;           % Upper bound of grid 150 degrees
space.dof{1}.mass  = inertia;            % Moment of inertia for the kinetic energy

% Hamiltonian operator 
hamilt.eigen.symmetry  = symmetry;       % get symmetry passed in by global variable

hamilt.truncate.min    = -1.0;           % Lower truncation of energy
hamilt.truncate.max    = 1.0;            % Upper truncation of energy

hamilt.pot.handle      = @pot.C9A;
hamilt.pot.params.state= 1;              % ground state surface of C9A

% Select eigen/values/functions
psi.eigen.start        =  0;             % Lower index
psi.eigen.stop         =  0;             % Upper index

% Modify settings for appearance of plots (if desired)
plots.density.type = 'contour';          % contour plot
