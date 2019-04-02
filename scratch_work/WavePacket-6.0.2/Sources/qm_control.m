%--------------------------------------------------------------------------
%
% Solves the bilinear control equation using standard ODE solvers
%
% The three input variables specify from which file to read 
%  A, B, N, and C, D matrices, initial/equilib state/density, etc.
%
% Input "filename" typically specifying physical context, e.g.
% 'lvne' for Liouville-von-Neumann
% 'tdse' for time-dependent Schroedinger equation
%
% If specified, input "method" specifies 
% b: balanced
% t: truncated
% s: singular perturbation theory
% r: H2 error reduction
% If not specified, then original/unmodified data will be used.
%
% If specified, input "reduce" specifies dimensionality
% of truncated/reduced model equations 
%
% In contrast to qm_optimal, this code is based on the 
% *coarse* time discretization given in time.main.whatever
% However, internally the ODE solvers may use finer
% grids, typically chosen adaptively, whereby +oct/rhs.m
% will be calling +util/efield.m as many times as necessary.
%
%--------------------------------------------------------------------------


% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2011 Boris Schaefer-Bung, Ulf Lorenz
%               2012 Jeremy Rodriguez, Ulf Lorenz
%               2013-16 Burkhard Schmidt
%
% see the README file for license details.

function qm_control (filename, method, reduce)

% Main variables are global throughout;
global bilinear control plots time  

% Initializes general information and sets up log files.
log.init (mfilename('fullpath'));

log.disp ('***************************************************************')
log.disp ('Numerical solution of bilinear control problem    ')
log.disp ('***************************************************************')
log.disp (' ')
log.disp ('given in terms of the matrices A, B, N and C, D             ')
log.disp ('using ordinary differential equations (ODE) solvers         ')
log.disp ('https://sf.net/p/wavepacket/matlab/wiki/Reference.Programs.qm_control  ')
log.disp ('                                                            ')
log.disp ('d                       m                                   ')
log.disp ('-- |x(t)> = A x(t) + i Sum u (t) ( N |x(t)> + |b >          ')
log.disp ('dt                     k=1  k       k           k           ')
log.disp ('                                                            ')
log.disp ('y (t) = < c | x(t)+x > + <x(t)+x | D | x(t)+x >, j = 1,...,p')
log.disp ('  j        j        e           e   j        e              ')
log.disp ('                                                            ')
log.disp ('where the input is defined in terms of control field(s) u(t)') 
log.disp ('where the output is defined in terms of observable(s) y(t)  ')
log.disp ('with state vector |x(t)> ==> |x(t)> - |x_e> (e=equilibrium) ')
log.disp ('and with with |b_k> = N_k |x_e>. Note that the spectrum of  ')
log.disp ('matrix A should be in left half of the complex plane thus   ')
log.disp ('ensuring stability in the absence of control fields.        ')
log.disp (' ')

%% Initialize time and fields

% Initialize temporal discretization; use *coarse* time-stepping
init (time.steps);
control.t.steps = time.steps.m_grid;
control.t.n = length(control.t.steps);

% Initialize e-field (coarse time-stepping)
efi.init;
if isfield (time,'pulse')
    time.efield.grid = efi.eval(control.t.steps); % used for plotting only
    
    % Adding a row to control.u.forward only for non-vanishing fields
    for p=1:length(time.efield.grid)
        if any(abs(time.efield.grid{p})>0)
            control.u.forward(p,:) = time.efield.grid{p}';
        end
    end
    control.u.dim = size(control.u.forward,1);
end
    
% Load A, B, N, and C, D matrices, initial/equilib state/density, etc.
switch nargin
    case 0
        log.error ('Please specify at least one input argument: filetype, e.g. tdse or lvne')
    case 1
        myfile = filename;
    case 2
        myfile = [filename '_' method];
    case 3
        myfile = [filename '_' method int2str(reduce)];
end
load (myfile)

log.disp (' ')
log.disp ('********************************************************************')
log.disp ('Bilinear control: Using MATLAB''s (or other) ODE solvers')
log.disp ('********************************************************************')

if ~isfield(control,'solvers')
    control.solvers=[];
end
if ~isfield(control.solvers,'handle2')
    control.solvers.handle2=@ode45;
end

switch lower(func2str(control.solvers.handle2))
    case {'ode113', 'ode15s', 'ode23', 'ode23s', 'ode23t', 'ode23tb','ode45'}
        log.disp ( ['ODE solver from Matlab built-in : ' func2str(control.solvers.handle2)])
    otherwise
        log.error ( ['Chosen ODE solver ' func2str(control.solvers.handle2) ...
            ' is not available.'] )
end
log.disp ( ' ' )

options = odeset ('reltol',control.solvers.reltol)

% Initial and equilibrium state
control.x.initial = bilinear.x.initial;
control.x.equilib = bilinear.x.equilib;
control.y.initial = bilinear.y.initial;
control.y.equilib = bilinear.y.equilib;

% Initialize calculation of observables
oct.observe('initial');

% Preallocate x(t), y(t)
control.x.dim = size(bilinear.A,1);
control.y.dim = bilinear.len.CD;
control.x.forward = zeros(control.x.dim, control.t.n);
control.y.forward = zeros(control.y.dim, control.t.n);

% Initialize plotting: Equilibrium values as horizontal lines
control.title{1} = 'Propagate forward';
control.title{2} = [bilinear.title int2str(length(control.x.initial)) ' coupled ODEs'];

control.plot.mov = false;

if ~isfield (control.plot,'uxy')
    control.plot.uxy=true;
end
log.disp (['Plot evolutions u(t), x(t), y(t) : ' int2str(control.plot.uxy)])

if control.plot.uxy
    initial (plots.control);
    equilib (plots.control);
end

%% Propagate forward in time
for step = 1:control.t.n
    
    % Initial state
    if step==1
        control.x.forward(:,1) = control.x.initial;
    else
        % Propagate one step forward by chosen ODE solver
        [~,x_new] = feval ( ...
            control.solvers.handle2, ...
            @oct.rhs2, ...
            [control.t.steps(step-1) control.t.steps(step)], ...
            control.x.forward(:,step-1), ...
            options );
        control.x.forward(:,step) = x_new(end,:).';
    end
    
    % Calculating observables y(t) and plot u(t), x(t), y(t)
    oct.observe ( 'forward', step );
    
    % Plotting u(t), x(t), y(t) for every step
    if control.plot.uxy && step>1
        forward (plots.control, step-1, step );
    end
    
end
if control.plot.uxy
    clearfig (plots.control);
    closefig (plots.control);
end

% Save time dependence of u,x,y-vectors and related quantities
save ([myfile '_control'], 'control')
log.disp (['Saving simulation data to file : ' myfile '_control.mat'])

% Output clock/date/time
log.clock;

end