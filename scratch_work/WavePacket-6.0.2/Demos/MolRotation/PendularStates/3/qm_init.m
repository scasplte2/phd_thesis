% Copyright (C) 2008 Ulf Lorenz, Burkhard Schmidt
%               2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init(~,Delta)
global hamilt space time

log.disp( '***********************************************' )
log.disp( 'Time evolution of pendular states              ' )
log.disp( '\sigma = 5 pulse and varying intensity         ' )
log.disp( ' J.Chem.Phys 110, 3870 (1999)                  ' )
log.disp( '***********************************************' )


% Note that the paper uses normalised units. We choose the rotational
% constant B = \hbar^2/(2mR^2) to be one a.u., so that the time units in
% the paper are the same numbers as the time units in a.u. used in the
% WavePacket run.

% One angular degree of freedom
space.dof{1}       = grid.legendre;      % Gauss-Legendre DVR in cos theta
space.dof{1}.label = 'cos \Theta';       % axes labl
space.dof{1}.R_0   = 1;                  % constant value for R
space.dof{1}.m_0   = 0;                  % minor quantum number fixed to 0
space.dof{1}.l_max = 25;                 % maximum angular momentum
space.dof{1}.mass  = 0.5;                % adjusted mass/inertia

% Alignment: cosine^2 projector
hamilt.amo{1} = amo.cosine;
hamilt.amo{1}.exp = 2;

% Temporal discretization
time.steps.m_start  = 000;               % Index of initial time step
time.steps.m_stop   = 500;               % Index of final time step
time.steps.m_delta  = 0.1;               % Size of big time steps
time.steps.s_number = 100;               % Number of sub steps per time step

% Propagator
time.propa       = tmp.splitting;        % Split-operator method
time.propa.order = 2;                    % Strang splitting

% Electric fild with properly designed parameters
time.pulse{1}       = efi.gauss;         % Gaussian pulse
time.pulse{1}.delay = 15;                % pulse delay, measured from paper
time.pulse{1}.fwhm  = 5*sqrt(8*log(2));  % FWHM of pulse
time.pulse{1}.ampli = sqrt(Delta);       % amplitude of field

% Hamiltonian operator: polarizability interaction 
hamilt.pol{1,1}{1,1} = pol.taylor;       % Taylor series in cos Theta
hamilt.pol{1,1}{1,1}.coeffs = [0;2];     % Set the polarizability to cos^2

% Initial wave function is a spherical harmonic of degree l=1.
time.dof{1}        = init.fbr;           % use eigenstate of grid (spherical harmonic)
time.dof{1}.state  = 1;                  % use first eigenstate (for m=0, this is l=0)