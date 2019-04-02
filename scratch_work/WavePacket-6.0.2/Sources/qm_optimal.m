%--------------------------------------------------------------------------
%
% Solves an optimal control problem by iterated forward-backward propagation
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
% In contrast to qm_control, this code is based on the 
% *fine* time discretization given in time.sub.whatever
% This implies that our ODE solvers are restriced to an
% equidistant time-stepping which makes forward-backward
% issues much easier to handle.
%
%--------------------------------------------------------------------------


% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2015-17 Burkhard Schmidt, FU Berlin
%
% see the README file for license details.

function qm_optimal (filename, method, reduce)

% Main variables are global throughout;
global bilinear control plots time  

% Initializes general information and sets up log files.
log.init (mfilename('fullpath'));

log.disp ('***************************************************************')
log.disp ('Numerical solution of the optimal control problem ')
log.disp ('***************************************************************')
log.disp (' ')
log.disp ('given in terms of the system matrices A, B, N and C, D      ')
log.disp ('                                                            ')
log.disp ('(1) Maximize output/observable/target at final time T       ')
log.disp ('J_1 [x] = 2 Re ( <c|x(T)> + <x_e|D|x(T)> ) + <x(T)|D|x(T)>  ')
log.disp ('with state vector x(t), called here x.forward               ')
log.disp ('                                                            ')
log.disp ('(2) Minimize cost/energy of input/control field             ')
log.disp ('J_2 [u] = \alpha \Int_0^T dt u^2 (t) / s(t)                 ')
log.disp ('with penalty factor alpha > 0                               ')
log.disp ('with shape function s(t)                                    ')
log.disp ('                                                            ')
log.disp ('(3) Minimize deviation from exact evolution                 ')
log.disp ('J_3 [u,x,z] = 2 Re \int_0^T dt <z(t)|d/dt-L|x(t)+x_e>       ')
log.disp ('with Lagrange multiplier z(t), called here x.backward       ')
log.disp ('                                                            ')
log.disp ('The combined target functional J = J_1 - J_2 - J_3 becomes  ') 
log.disp ('extremal if the following three conditions are satisfied    ')
log.disp ('(which can be derived by standard variational calculus)     ')
log.disp ('                                                            ')
log.disp ('State vector: propagated forward from x(0) = x_0 - x_e      ')
log.disp ('                                                            ')
log.disp ('d                                 m                         ')
log.disp ('--|x(t)> =  L |x(t)+x > = A + i Sum u (t) ( N |x(t)> + b >  ')
log.disp ('dt                   e          k=1  k       k          k   ')
log.disp ('                                                            ')
log.disp ('Lagrange multiplier: propagate backward from z(T) = C+Dx(T) ')
log.disp ('                                                            ')
log.disp ('d           +            +        +                         ')
log.disp ('--|z(t)> =-L |z(t)> = (-A + iu(t)N ) |z(t)> + iu(t)|b>      ')
log.disp ('dt                                                          ')
log.disp ('                                                            ')
log.disp ('with |b> = N |x_e>. Note that for the special case of       ')
log.disp ('anti-Hermition L (e.g. TDSE where L = -i/hbar H):           ')
log.disp ('Same propagation for x and z which can be tested by fb_test ')
log.disp ('                                                            ')
log.disp ('Control field                                               ')
log.disp ('                                                            ')
log.disp ('        -s(t)                                               ')
log.disp ('u (t) = ----- Im < z(t)| N | x(t) + x >                     ')
log.disp ('        alpha                        e                      ')
log.disp ('                                                            ')
log.disp ('This system of coupled equations is solved here by iterative')
log.disp ('forward-backward propagation which is known to be rapidly   ')
log.disp ('monotonically convergent (see Rabitz Ohtsuki... JCP 1998 ff.')
log.disp ('For more details see                                        ')
log.disp ('https://sf.net/p/wavepacket/matlab/wiki/Reference.Programs.qm_optimal  ')
log.disp (' ')

%% Initialize time and fields

% Initialize temporal discretization; use *fine* time-stepping
init (time.steps);
control.t.steps = time.steps.s_grid;
control.t.delta = time.steps.s_delta;
control.t.n = length(control.t.steps);

% Initialize e-field (fine time-stepping)
efi.init;
time.efield.dressed = true; % Tricky fake
time.efield.grid = efi.eval(control.t.steps); % half envelope
time.efield.dressed = false;

% Adding a row to control.u.forward etc only for non-vaishing fields 
for p=1:length(time.efield.grid)
    if any(abs(time.efield.grid{p})>0)
        control.u.shaping (p,:) = time.efield.grid{p}'*2/time.efield.max_ampli;
        control.u.forward (p,:) = time.efield.max_ampli * ones(1,control.t.n);
        control.d.forward (p,:) = zeros(1,control.t.n); % du/dt
        control.d.backward(p,:) = zeros(1,control.t.n); % du/dt
    end
end
control.u.dim = size(control.u.forward,1);

% Initialize frequency resolved optical gating (FROG)
if ~isfield(time,'frog')
    time.frog=[];
end
if ~isfield(time.frog,'choice')
    time.frog.choice='none';
end
log.disp (['Which type of FROG to be performed : ' time.frog.choice])

if strcmpi(time.frog.choice,'gauss')
    if ~isfield(time.frog,'width')
        time.frog.width=time.steps.m_delta;
    end
    log.disp (['Width of Gaussian : ' num2str(time.frog.width)])
end

if ~isfield(time.frog,'zoom')
    time.frog.zoom=1;
end
log.disp (['Zoom factor for PSD and FROG frequencies : ' num2str(time.frog.zoom)])
log.disp ('   ')


%% Check and - where necessary - set OCT variables
log.disp('********************************************************************')
log.disp('Settings for optimal control')
log.disp('********************************************************************')

if ~isfield(control,'optimal')
    control.optimal=[];
end

if ~isfield(control.optimal,'terminal')
    control.optimal.terminal=2;
end
log.disp (['Which observable to be optimized : ' int2str(control.optimal.terminal)])
if length(control.optimal.terminal)>1
    log.error('Multi-target optimization not yet(!) implemented')
end

if ~isfield(control.optimal,'tolerance')
    control.optimal.tolerance=1e-3;
end
log.disp (['Tolerance terminating iteration  : ' num2str(control.optimal.tolerance)])

if ~isfield(control.optimal,'max_iter')
    control.optimal.max_iter=0;
end
log.disp (['Maximum number of iteration step : ' int2str(control.optimal.max_iter)])

if ~isfield(control.optimal,'alpha')
    control.optimal.alpha=ones(size(control.u.forward,1),1);
end
log.disp (['Penalty factor for laser fluence : ' num2str(control.optimal.alpha)])

if ~isfield(control.optimal,'eta')
    control.optimal.eta=1;  % no admixture of fields from previous step
end
log.disp (['Eta parameter for backward field : ' num2str(control.optimal.eta)])

if ~isfield(control.optimal,'zeta')
    control.optimal.zeta=1; % no admixture of fields from previous step
end
log.disp (['Zeta parameter for forward field : ' num2str(control.optimal.zeta)])

if ~isfield(control.optimal,'order')
    control.optimal.order=1; 
end
log.disp (['Order for calculations of field  : ' int2str(control.optimal.order)])
if control.optimal.order>2
    log.error(['Wrong choice for order parameter : ' int2str(control.optimal.order)])
end

if ~isfield(control.optimal,'prefactor')
    control.optimal.prefactor='current'; 
end
log.disp (['What time to calculate prefactor : ' control.optimal.prefactor])

if ~isfield(control.optimal,'fb_test')
    control.optimal.fb_test=false;
end
log.disp (['Test only forward-backward propa : ' int2str(control.optimal.fb_test)])
if control.optimal.fb_test
    control.u.backward = control.u.forward;
    control.optimal.max_iter=1;
end

if ~isfield (control.plot,'uxy')
    control.plot.uxy=true;
end
log.disp (['Plot evolutions u(t), x(t), y(t) : ' int2str(control.plot.uxy)])

if ~isfield (control.plot,'j12')
    control.plot.j12=true;
end
log.disp (['Plot functionals j_1, j_2, j_tot : ' int2str(control.plot.j12)])

if ~isfield (control.plot,'psd')
    control.plot.psd=true;
end
log.disp (['Plot power spectral density      : ' int2str(control.plot.psd)])

if ~isfield (control.plot,'mov')
    control.plot.mov=false;
end
log.disp (['Create animation of u, x, y(t)   : ' int2str(control.plot.mov)])


if ~isfield(control,'solvers')
    control.solvers=[];
end
if ~isfield(control.solvers,'handle1')
    control.solvers.handle1=@ode.RuKu4;
end
log.disp (['ODE solver (constant step size)  : ' func2str(control.solvers.handle1)])
log.disp ( ' ' )

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


% Initial and equilibrium state
control.x.initial = bilinear.x.initial;
control.x.equilib = bilinear.x.equilib;
control.y.initial = bilinear.y.initial;
control.y.equilib = bilinear.y.equilib;

% Initialize calculation of observables 
oct.observe( 'initial' );

% Preallocate x(t) and y(t), both for forward and backward
control.x.dim = size(bilinear.A,1);
control.x.forward  = zeros(control.x.dim, control.t.n);
control.x.backward = zeros(control.x.dim, control.t.n);

control.y.dim = bilinear.len.CD;
control.y.forward  = zeros(control.y.dim, control.t.n);
control.y.backward = zeros(control.y.dim, control.t.n);

control.title{1}=[...
    bilinear.label{control.optimal.terminal} ...
    ': ord=' ...
    int2str(control.optimal.order) ...
    ', u_0=' ...
    num2str(time.pulse{1}.ampli) ...
    ' (' ...
    class(time.pulse{1}) ...
     '), \alpha=' ...
    num2str(control.optimal.alpha) ...
    ', \eta=' ...
    num2str(control.optimal.eta) ...
    ', \zeta=' ...
    num2str(control.optimal.zeta)];
if isfield(bilinear,'C') 
    if ~bilinear.Q{control.optimal.terminal}
        control.title{1} = [control.title{1} ' (C1)'];
    else
        control.title{1} = [control.title{1} ' (C2/' control.optimal.prefactor(1) ')'];
    end
elseif isfield(bilinear,'D')
    control.title{1} = [control.title{1} ' (D)'];
end
control.title{2} = [bilinear.title int2str(length(control.x.initial)) ' coupled ODEs'];
control.title{3} = 'Initialization: propagate  forward';

%% Initialize OCT scheme: Propagate FORWARD in time using "guessed" fields

% Initialize plotting: Equilibrium values as horizontal lines
if control.plot.uxy
    initial (plots.control);
    equilib (plots.control);
end
    
% Propagate FORWARD in time
for step = 1:control.t.n
    
    % Initial state
    if step==1
        control.x.forward(:,1) = control.x.initial;

    else
        % Propagate x(t) forward by chosen ODE solver: step-1 ==> step
        control.x.forward(:,step) = feval ( ...
            control.solvers.handle1, ...
            @oct.rhs1x, ...
            control.x.forward(:,step-1), ...
            control.u.forward(:,step-1), ...
            0, ...
            control.t.delta);
    end
    
    % Calculating observables y(t) for every step
	oct.observe ( 'forward', step );
    
	% Plotting u(t), x(t), y(t)  for every n-th step
    if control.plot.uxy && step>1 && mod(step-1,time.steps.s_number)==0
        forward (plots.control, step-time.steps.s_number, step );
    end
    
end
% log.disp(['norm after forward propagation = ' num2str(norm(control.x.forward(:,1)))])

% Calculate functionals and write to console and logfile
[j_ini] = oct.functionals (0);

% Preallocate
control.j.target = zeros(1,control.optimal.max_iter);
control.j.cost   = zeros(1,control.optimal.max_iter);
control.j.total  = zeros(1,control.optimal.max_iter);
control.j.addup  = zeros(1,control.optimal.max_iter);

%% Main OCT iteration loop
inc = inf;
iter = 0;
while( inc>control.optimal.tolerance && iter<control.optimal.max_iter )
    iter = iter+1;
    
    % Save results from previous forward propagation
    control.u.previous = control.u.forward;
    control.x.previous = control.x.forward;
    
    % Reset plots: Equilibrium values as horizontal lines
    if control.plot.uxy
        clearfig ( plots.control );
        equilib  ( plots.control );
    end
        
    %% Propagate BACKWARD in time
    control.title{3} = [ int2str(iter) '-th iteration: propagate backward' ];
    for step = control.t.n:-1:1
        
        % Various options to initialize backward propagation
        if step==control.t.n
            if control.optimal.fb_test % Start from last step of forward propagation
                control.x.backward(:,end) = control.x.forward(:,end);
            else
                if isfield(bilinear,'C') % linear target 
                    control.x.backward(:,end) = bilinear.C{control.optimal.terminal}';
                elseif isfield(bilinear,'D') % quadratic target
                    control.x.backward(:,end) = bilinear.D{control.optimal.terminal} * ...
                        control.x.forward(:,end);
                end
            end
            
        else
            
            % Propagate z(t) BACKWARD by chosen ODE solver: step+1 ==> step
            control.x.backward(:,step) = feval ( ...
                control.solvers.handle1, ...
                @oct.rhs1z, ...
                control.x.backward(:,step+1), ...
                control.u.backward(:,step+1), ...
                control.d.backward(:,step+1), ...
               -control.t.delta);
        end
        
        % Optimal control field (and derivative) for use in next time step
        if ~control.optimal.fb_test
            oct.u_opt ('backward',step);
            oct.u_dot ('backward',step);
        end
        
        % Calculating observables y(t) for every step
        oct.observe ('backward', step );
        
        % Plotting u(t), x(t), y(t) for every n-th step
        if control.plot.uxy && step<control.t.n && mod(step-1,time.steps.s_number)==0
            backward (plots.control, step, step+time.steps.s_number );
        end
        
    end    
    % log.disp(['norm after backward propagation = ' num2str(norm(control.x.backward(:,end)))])
    
    %% Propagate FORWARD in time
    control.title{3} = [ int2str(iter) '-th iteration: propagate  forward' ];
    
    for step = 1:control.t.n
        
        % Initial state
        if step==1
            control.x.forward(:,1) = control.x.initial;
            
        else
            % Propagate x(t) FORWARD by chosen ODE solver: step-1 ==> step
            control.x.forward(:,step) = feval ( ...
                control.solvers.handle1, ...
                @oct.rhs1x, ...
                control.x.forward(:,step-1), ...
                control.u.forward(:,step-1), ...
                control.d.forward(:,step-1), ...
                control.t.delta);
        end

        % Optimal control field (and derivative) for use in next time step
        if ~control.optimal.fb_test
            oct.u_opt ('forward',step);
            oct.u_dot ('forward',step);
        end
        
        % Calculating observables y(t) for every step 
        oct.observe ( 'forward', step );
        
        % Plotting u(t), x(t), y(t) for every n-th step
        if control.plot.uxy && step>1 && mod(step-1,time.steps.s_number)==0
            forward (plots.control, step-time.steps.s_number, step );
        end
        
    end
    % log.disp(['norm after forward propagation = ' num2str(norm(control.x.forward(:,1)))])
        
    
    % Calculate functionals and write / plot them
    [inc] = oct.increment (iter);
    [j12] = oct.functionals (iter);
    control.j.target(iter) = j12.target;
    control.j.cost  (iter) = j12.cost;
    control.j.total (iter) = j12.total;
    if iter==1
        control.j.addup (1) = j_ini.total + inc;
    else
        control.j.addup (iter) = control.j.addup (iter-1)+inc;
    end
    if control.plot.j12 && control.optimal.max_iter>1
        oct.plot_j12 (iter)
    end
    
    %  Calculate power spectral density and FROG
    oct.spectrum;
    if control.plot.psd
        oct.plot_psd (iter)
    end
    
end
log.disp('   ')
    log.disp('*********************************************************************')
if iter==control.optimal.max_iter
    log.disp('Terminate iteration after maximum number of iterations')
else
    log.disp('Terminate iteration because increment is below tolerance')
end
    log.disp('*********************************************************************')

if control.plot.uxy
    clearfig (plots.control);
    closefig (plots.control);
end
    
% Save time dependence of u,x,y-vectors and related quantities
save ([myfile '_optimal'], 'control')
log.disp (['Saving simulation data to file  : ' myfile '_optimal.mat'])

% Save time dependence of final u(t) only
data = zeros (control.t.n,2);
data(:,1) = control.t.steps;
data(:,2) = control.u.forward;
save ([myfile '_optimal_' int2str(iter) '.dat'], 'data','-ascii')
log.disp (['Saving optimized field to file  : ' myfile '_optimal_' int2str(iter) '.dat'])

% Output clock/date/time
log.clock;

end