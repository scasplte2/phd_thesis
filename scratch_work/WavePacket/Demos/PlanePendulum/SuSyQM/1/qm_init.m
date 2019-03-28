% Copyright (C) 2008-2009 Ulf Lorenz
%               2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots space

util.disp ( '***************************************' )
util.disp ( 'Supersymmetry and eigensurface topology' )
util.disp ( 'of the planar quantum pendulum' )
util.disp ( 'B. Schmidt and B. Friedrich' )
util.disp ( 'Front. Phys. 2, 37' )
util.disp ( 'DOI:10.3389/fphy.2014.00037' )
util.disp ( 'Reproducing red circles in Fig. 4' )
util.disp ( 'USING HERE: FFT-DVR for THETA' )
util.disp ( '***************************************' )

% Spatial discretization
space.dof{1}       = grid_fft;           % Fourier grid
space.dof{1}.n_pts = 256;                % Number of grid points
space.dof{1}.x_min = -1*pi;              % Lower bound of grid
space.dof{1}.x_max = +1*pi;              % Upper bound of grid
space.dof{1}.mass  =  1/2;               % Mass for the kinetic energy

% cosine^2 projector
space.prj.handle = @prj.cosine;
space.prj.params.exp = 2;

% Hamiltonian operator 
hamilt.pot.handle      = @pot.pendulum;  % potential for generalized pendula
hamilt.pot.params.eta  =  5;             % Orientation: cos
hamilt.pot.params.zeta = 25;             % Alignment: cos^2

% Modify settings for appearance of plots (if desired)
plots.density.type = 'contour';          % Wigner contour plot
plots.pot.max = 100;                     % customize density plot
plots.expect.energies.max = 100;         % Set maximum for energy plot
