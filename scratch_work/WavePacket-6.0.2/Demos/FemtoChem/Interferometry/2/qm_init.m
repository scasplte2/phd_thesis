% Copyright (C) 2017 - .... Burkhard Schmidt 
%               2009 - 2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init (~)
global atomic hamilt plots space time

log.disp ('****************************************')
log.disp ('Calculation of the pump only signal for ')
log.disp ('excitation of Na2.                      ')
log.disp ('see J.Chem.Phys. 100:5448               ')
log.disp ('****************************************')

% Number of coupled equations
hamilt.coupling.n_eqs  = 3;
hamilt.coupling.represent = 'dia';
hamilt.coupling.labels = {'X^1\Sigma_g^+', 'A^1\Sigma_u^+', '2^1\Pi_g'};

% Grid definition
space.dof{1}       = grid.fft;           % using FFT grid
space.dof{1}.n_pts = 256;                % number of points
space.dof{1}.x_min = 4;                  % lower boundary of the grid
space.dof{1}.x_max = 14;                 % upper boundary of the grid
space.dof{1}.mass  = atomic.mass.Na23/2; % reduced mass of 23Na2

% Temporal discretisation
time.steps.m_start  = 0;                 % index of first time step
time.steps.m_stop   = 210;               % index of last time step
time.steps.m_delta  = 10/atomic.t.fs;    % 10 fs per time step
time.steps.s_number = 500;               % propagation time step 20 as

% Propagator
time.propa       = tmp.splitting;        % Split operator
time.propa.order = 2;                    % Strang splitting

% Electric field 
time.pulse{1}       = efi.gauss;         % Gaussian-shaped pulse
time.pulse{1}.fwhm  = 30/atomic.t.fs;    % with 30 fs FWHM
time.pulse{1}.delay = 100/atomic.t.fs;   % delayed by 100 fs
time.pulse{1}.frequ = 0.073;             % frequency (corr. to 625 nm)
time.pulse{1}.ampli = sqrt(10/atomic.I.GW_cm2);

% Hamiltonian operator
hamilt.truncate.e_min = -0.03;           % lower truncation of energy
hamilt.truncate.e_max =  0.25;           % upper truncation of energy

for m=1:3
    hamilt.pot{m,m}          = pot.interp;   % interpolate tabulated potential
    hamilt.pot{m,m}.pos_conv = atomic.l.A;   % conversion factor for coordinates
    hamilt.pot{m,m}.pot_conv = atomic.E.eV;  % conversion factor for energies
end

hamilt.dip{1}             = cell(3);
hamilt.dip{1}{1,2}        = dip.taylor;  % transition dipole moment
hamilt.dip{1}{1,2}.vshift = 1;           % constant ONE
hamilt.dip{1}{2,3}        = dip.taylor;  % transition dipole moment
hamilt.dip{1}{2,3}.vshift = 1;           % constant ONE

% Initial wave function; Interpolate from input file.
time.corr = init.interp; 

% Plot densities
plots.density           = vis.contour;   % Draw contour lines
plots.density.expect    = false;
plots.density.cnt_min   = -0.01;         % bounds for the contour lines
plots.density.cnt_max   = +0.01;

% Plot expectation values
plots.expect       = vis.expect;
plots.expect.e_max = 0.17;




