% Copyright (C) 2016 Burkhard Schmidt
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots space

util.disp ( '*********************************************' )
util.disp ( 'Supersymmetry and eigensurface topology' )
util.disp ( 'of the spherical quantum pendulum' )
util.disp ( 'B. Schmidt and B. Friedrich' )
util.disp ( 'Phys. Rev. A 91, 022111' )
util.disp ( 'DOI:10.1103/PhysRev.A.91.022111' )
util.disp ( 'Reproducing red circles in Fig. 3' )
util.disp ( 'USING HERE: Gauss-Legendre-DVR for COS(THETA)' )
util.disp ( '*********************************************' )

% Spatial discretization
space.dof{1} = grid_legendre;            % Gauss-Legendre DVR in cos(theta)
space.dof{1}.label = 'cos \Theta';
space.dof{1}.R_0 = 1;                    % constant value for R
space.dof{1}.m_0 = 1;                    % minor quantum number 
space.dof{1}.l_max = 100;                % maximum angular momentum/ number of points
space.dof{1}.mass = 0.5;                 % adjusted mass

% cosine^2 projector
space.prj.handle = @prj.cosine;
space.prj.params.exp = 2;

% Hamiltonian operator 
% hamilt.eigen.symmetry = 'u';             % Symmetry adaption for eta=0 only
hamilt.pot.handle      = @pot.taylor;    % Taylor series in cos(theta)
hamilt.pot.params.v {1,1} = [-40;-2*100];% eta and zeta parameters

% Modify settings for appearance of plots (if desired)
plots.density.type = 'curve';            % polar plot
plots.pot.max = 200;                     % customize density plot
plots.expect.energies.max = 300;         % Set maximum for energy plot
