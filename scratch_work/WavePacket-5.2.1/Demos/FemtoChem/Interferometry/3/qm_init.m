% Copyright (C) 2009,2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init(delay, states)
global hamilt plots psi space time

util.disp ('****************************************')
util.disp ('Calculation of the pump probe signal for')
util.disp ('excitation of Na2.                      ')
util.disp ('see J.Chem.Phys. 100:5448               ')
util.disp ('****************************************')

% Note that the function argument probe gives the time delay between
% pump and probe in times of 10 fs. This affects the calculation in a number of
% ways:
%
% 1. We place the 30 fs probe pulse at the given time delay.
% 2. We propagate from 100 fs before the probe to 100 fs after probe.
% 3. We take the initial wavefunction from the pump-only calculation.
%
% 2. and 3. speed up the calculation a _lot_. (Delays are up to 2 ps!)
probe_time = delay + 10;    % Necessary for some index computations. Takes care
                            % of the fact that the pump pulse is centered at
                            % 100 fs.


% Number of coupled equations
hamilt.coupling.n_eqs  = 3;
hamilt.coupling.representation = 'dia';
hamilt.coupling.labels = {'X^1\Sigma_g^+', 'A^1\Sigma_u^+', '2^1\Pi_g'};

% Grid definition
space.dof{1}       = grid_fft;  % using FFT grid
space.dof{1}.label = 'R';       % label for axis descriptions
space.dof{1}.n_pts = 256;       % number of points
space.dof{1}.x_min = 4;         % lower boundary of the grid
space.dof{1}.x_max = 14;        % upper boundary of the grid
space.dof{1}.mass  = 23 * 1836/2; % reduced mass of Na2

% Temporal discretisation
time.main.start = probe_time - 10;  % index of first time step
time.main.stop  = probe_time + 10;  % index of last time step

time.main.delta = 10 * 41.34;   % 10 fs per time step
time.sub.n      = 500;          % propagation time step 20 as

% Propagator
time.propa.handle = @ket.splitting; % Split operator
time.propa.params.order = 3;        % Strang splitting

% Electric field is a sin^2 pulse with ~25 fs FWHM
time.efield.shape   = ['gauss'; 'gauss'];   % two almost identical Gauss pulses
time.efield.fwhm    = [30 30] * 41.34;      % with 30 fs FWHM
time.efield.delay   = [100 10*probe_time] * 41.34;  % delayed by 100 fs
time.efield.polar   = [0 0];                % fixed orientation => irrelevant
time.efield.phase   = [0 -0.073*(delay*413.41105)]; % locked phase
time.efield.frequ   = [0.073 0.073];        % frequency (corr. to 625 nm)
time.efield.ampli   = [1 1] * sqrt(1e10/3.51e16); % 1e10 W/cm^2 both

% Hamiltonian operator
hamilt.truncate.min = -0.03;    % lower truncation of energy
hamilt.truncate.max =  0.25;    % upper truncation of energy

hamilt.pot.handle = @pot.interp;       % potential energy
hamilt.pot.params.pos_conv = 0.5292;   % conversion factor for coordinates
hamilt.pot.params.pot_conv = 27.212;   % conversion factor for energies
hamilt.pot.params.method   = 'spline'; % spline interpolation

hamilt.dip.handle = @FC_dip;    % dipole moments; defined further below.

% Initial wave function; Load from pump-only run.
psi.init.corr.handle = @wav.load;   % load from file
psi.init.corr.dir    = '../2/';
psi.init.corr.file   = 'pump';
psi.init.corr.index  = probe_time + 1; % index 1 is the initial wave function!
psi.init.corr.state  = states;

psi.init.representation = 'dia';
psi.init.norm        = false;       % do not normalize the initial wave function

% Tirn off plotting; they take too much time
plots.density.on = false;
plots.expect.on  = false;

end  % qm_init


% Now specify a function that gives a simple Franck-Condon transition for
% the different states. Unfortunately, the reference does not give values here,
% so we just invent some vaguely based on reality.
function FC_dip
    global hamilt space
    hamilt.d_x.grid_ND{1,1} = [];
    hamilt.d_x.grid_ND{2,2} = [];
    hamilt.d_x.grid_ND{3,3} = [];
    hamilt.d_x.grid_ND{1,2} = ones(size(space.dvr.grid_ND{1}));
    hamilt.d_x.grid_ND{2,3} = ones(size(space.dvr.grid_ND{1}));
end

