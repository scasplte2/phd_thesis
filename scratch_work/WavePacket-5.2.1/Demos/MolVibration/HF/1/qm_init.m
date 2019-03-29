% Copyright (C) 2009-2010 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots psi space time


util.disp('*********************************************')
util.disp('Control of vibrational modes of HF via       ')
util.disp('chirped infrared laser pulses.               ')
util.disp('                                             ')
util.disp('Figure 1 in the reference                    ')
util.disp('                                             ')
util.disp('see S. Chelkowski, A.D. Bandrauk, P.B. Corkum')
util.disp('Phys. Rev. Lett. 65, page 2355               ')
util.disp('*********************************************')

% Define the grid
space.dof{1}       = grid_fft;  % use an FFT grid
space.dof{1}.label = 'R';
space.dof{1}.mass  = 0.95*1823; % reduced mass of HF
space.dof{1}.n_pts = 2048;      % number of points
space.dof{1}.x_min = 0.5;       % lower boundary; not given in paper
space.dof{1}.x_max = 130.5;     % upper boundary

% Simple Morse potential
hamilt.pot.handle     = @pot.morse;
hamilt.pot.params.d_e = 6.125 / 27.21; % dissociation energy
hamilt.pot.params.r_e = 1.7329;        % equilibrium distance
hamilt.pot.params.alf = 1.1741;        % range parameter
hamilt.pot.params.t_e = 0;             % no energy shift

% linear dipole moment; simple, but not implemented, so we load it from an
% external file
hamilt.dip.handle          = @dip.interp;
hamilt.dip.params.pos_conv = 1;        % conversion factor for coordinates
hamilt.dip.params.dip_conv = 1;        % conversion factor for dipole moments
hamilt.dip.params.method   = 'spline'; % interpolation method
hamilt.dip.params.n_pts    = 131;      % number of points in data file

% truncation
hamilt.truncate.min = 0;
hamilt.truncate.max = 1;

% time step parameters
time.main.start = 0;            % index of first time step
time.main.stop  = 120;          % index of last time step
time.main.delta = 8.41 * 41.34; % main step has units of laser cycles (8.41 fs)
time.sub.n      = 1000;         % number of propagation steps per main step

% the propagator
time.propa.handle       = @ket.splitting;
time.propa.params.order = 3;    % Strang splitting

% the electric field is composed of four pulses. The turn-on pulse is
% interpolated from a file. After that, we have a plateau, and chirp the
% frequency up to time t_N. We continue for some time without chirping the
% pulse, and the final turn-off pulse is again interpolated from file.
% A major complication in the chirped part of the pulse is the fact that the
% reference expresses the chirp as a function of (t - t_1), while the WavePacket
% implementation wants to have the chirp as a function of (t - delay time), and
% also multiplies the time-dependent frequency with this shifted time,
% where the delay time is the center between t_0 (when the chirping starts) and
% t_8 (when chirping ends). Also, the phases of the electric fields are always
% defined relative to the time delay, so static phase shifts have to be added to
% avoid jumps in the electric field when switching from one field to the next.

% some values we need for the calculations; They can be extracted from the
% reference.
A  = -0.5224;
B  = 0.0419;
Q  = 7316;
t0 = 9513;
t1 = 12286;
t8 = t1 + 2 * Q * (sqrt(9) - sqrt(2));   % eq. 12 in the reference
w0 = 2 * B * hamilt.pot.params.d_e;      % harmonic frequency of the Morse potential
wi = w0 * (1 - B);                       % initial frequency of the pulse
wf = w0 - w0 * B * ( (t8-t1)^2 / (4*Q^2) + sqrt(2) * (t8-t1) / Q + 3/2);
                                         % final frequency of the pulse

T = (t0 + t8) /2;                        % delay time of the chirped pulse
quadratic = -B*w0 / (2*Q^2);             % quadratic chirp
linear    = -B*w0 / (4*Q^2) * (3*T - 2*t1) - sqrt(2)*B*w0 / Q;
                                         % linear chirp
constant  = -B*w0 / (4*Q^2) * (T^2 - t1^2 + 2*(T-t1)^2) ...
            - sqrt(2)*B*w0 / Q * (2*T - t1) + (w0 - 3/2*B*w0);
                                         % time-independent part of chirped frequency
% static phase of the second and third pulse to ensure a continuous field.
phi2 = -B * w0 / (4*Q^2) * T * (T-t1)^2 ...
        - sqrt(2) * B * w0 / Q * T * (T - t1) + (w0 - B * w0*3/2) * T;
phi3 = wf* t8 + wf*t0/10;

time.efield.shape     = ['inter'; 'recta'; 'recta'; 'inter'];   % pulse shapes
time.efield.ampli     = ones(4,1) * sqrt(1e13/3.5101e16);       % 1e13 W/cm^2
time.efield.delay     = [t0; T; t8+t0/10; t8+t0/5];             % delay times
time.efield.frequ     = [wi; constant; wf; wf];                 % central frequency
time.efield.linear    = [0; linear; 0; 0];                      % linear chirp
time.efield.quadratic = [0; quadratic; 0; 0];                   % quadratic chirp
time.efield.phase     = [wi*t0; phi2; phi3; phi3 + wf*t0/10];   % match phases
time.efield.polar     = [0; 0; 0; 0];                           % polarization
time.efield.fwhm      = [0; t8-t0; t0/5; 0];                    % FWHM (/length) of pulses
time.efield.file      = {'turn_on.dat'; ''; ''; 'turn_off.dat'};% files to load from
time.efield.method    = 'spline';                               % interpolation method


% Initial wave function is the ground state of the Morse oscillator
%  The values can be taken either from the grid or from the potential.
psi.init.dof{1}.handle = @wav.morse;
psi.init.dof{1}.m_r    = space.dof{1}.mass;     % reduced mass
psi.init.dof{1}.r_e    = hamilt.pot.params.r_e; % equilibrium distance
psi.init.dof{1}.d_e    = hamilt.pot.params.d_e; % dissociation energy
psi.init.dof{1}.alf    = hamilt.pot.params.alf; % range parameter
psi.init.dof{1}.n_q    = 0;                     % we want the groundstate

% Save the wave function to a file. We will extract data from it in the
% rundemo.m script.
psi.save.export = true;
psi.save.dir    = '.';
psi.save.file   = 'HF';


% Plots: Draw only the wave function, not the Wigner-transform (takes ages)
plots.density.type        = 'curve';

plots.density.range.on    = true;     % custom plotting ranges
plots.density.range.x_min = 0.5;
plots.density.range.x_max = 6;

plots.expect.energies.max = 0.2;
