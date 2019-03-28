% Copyright (C) 2016 Burkhard Schmidt
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots space

util.disp ( '***************************************' )
util.disp ( 'Supersymmetry and eigensurface topology' )
util.disp ( 'of the spherical quantum pendulum' )
util.disp ( 'B. Schmidt and B. Friedrich' )
util.disp ( 'Phys. Rev. A 91, 022111' )
util.disp ( 'DOI:10.1103/PhysRev.A.91.022111' )
util.disp ( 'Reproducing red circles in Fig. 3' )
util.disp ( 'USING HERE: FFT-DVR for THETA' )
util.disp ( '***************************************' )

% Spatial discretization
space.dof{1}       = grid_fft;           % Fourier grid
space.dof{1}.n_pts = 1024;               % Number of grid points
space.dof{1}.x_min = 0;                  % Lower bound of grid
space.dof{1}.x_max = +1*pi;              % Upper bound of grid
space.dof{1}.mass  =  1/2;               % Mass for the kinetic energy

% cosine^2 projector
space.prj.handle = @prj.cosine;
space.prj.params.exp = 2;

% Hamiltonian operator 
hamilt.pot.handle      = @pot.pendulum;  % potential for generalized pendula
hamilt.pot.params.xi = 1^2 - 1/4;        % Azimuthal rotation
hamilt.pot.params.eta  =  40;            % Orientation: cos
hamilt.pot.params.zeta = 100;            % Alignment: cos^2
hamilt.pot.params.v_0 = -1/4;            % Energy shift

% Modify settings for appearance of plots (if desired)
plots.density.type = 'curve';            % polar plot
plots.pot.max = 200;                     % customize density plot
plots.expect.energies.max = 300;         % Set maximum for energy plot
