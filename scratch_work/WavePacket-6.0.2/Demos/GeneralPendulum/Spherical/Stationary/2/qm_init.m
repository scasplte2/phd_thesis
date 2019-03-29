% Copyright (C) 2016 Burkhard Schmidt
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init (~,b,m)
global hamilt plots space

log.disp ( '*********************************************' )
log.disp ( 'Supersymmetry and eigensurface topology' )
log.disp ( 'of the spherical quantum pendulum' )
log.disp ( 'B. Schmidt and B. Friedrich' )
log.disp ( 'Phys. Rev. A 91, 022111' )
log.disp ( 'DOI:10.1103/PhysRev.A.91.022111' )
log.disp ( 'Reproducing red circles in Fig. 3' )
log.disp ( 'USING HERE: Gauss-Legendre-DVR for COS(THETA)' )
log.disp ( '*********************************************' )

% Spatial discretization
space.dof{1}       = grid.legendre;      % Gauss-Legendre DVR in cos(theta)
space.dof{1}.label = 'cos \Theta';
space.dof{1}.R_0   = 1;                  % constant value for R
space.dof{1}.m_0   = m;                  % minor quantum number 
space.dof{1}.l_max = 100;                % maximum angular momentum/ number of points
space.dof{1}.mass  = 0.5;                % adjusted mass

% Orientation
hamilt.amo{1} = amo.cosine;              % cosine projector
hamilt.amo{1}.exp = 1;                   % exponent

% Alignment
hamilt.amo{2} = amo.cosine;              % cosine^2 projector
hamilt.amo{2}.exp = 2;                   % exponent

% Pendular potential
hamilt.pot{1,1} = pot.taylor;            % Taylor series in cos(theta)
hamilt.pot{1,1}.coeffs = [-2*b*(m+1);-2*b^2];% eta and zeta parameters

% Plot time evolution of density
plots.density         = vis.curve;       % Colored curve plot
plots.density.pot_max = 200;

% Plot expectation values
plots.expect       = vis.expect;
plots.expect.p_min = -0.3;               % customize population/amo plot
plots.expect.p_max = +1.1;               % customize population/amo plot
plots.expect.e_min = -150;               % customize energy plot
plots.expect.e_max = +150;               % customize energy plot
