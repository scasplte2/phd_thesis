function solOut = sr_1S0_scatter(Tinc, el, rRange, masses, varargin)
% Function to solve for the collisional wavefunction using a simple Numerov propagator. 
% The calculation is done in atomic units but some inputs are specified in easier to use 
% (or more typical AMO) units. See the INPUTS section for details.
%
% Reference on atomic units
% https://onlinelibrary.wiley.com/doi/pdf/10.1002/9781118229101.app5
% Be careful of masses, mass in atomic units (a.u.) is electron mass. Atomic mass units (a.m.u.) specifies proton mass
%
% Sample input
% sr_1S0_scatter(1e-6, 0, [5 15e3], [84 84])
% sr_1S0_scatter(1e-6, 0, [5 15e3], [84 84], 'coeffs', C, 'plot', 1)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
%--------------------------------------------------------------------------
%   Name    | Input Unit  | Description
%-----------|-------------|------------------------------------------------
%   Tinc    | Kelvin      | Incoming energy of the wavefunction - Kelvin
%   el      | (no dim)    | Partial wave of collision to consider
%   rRange  | Bohr radius | two element vector specifying start and end points of propagation, [rStart rEnd]
%   masses  | a.m.u.      | two element vector specifying masses of the isotopic combination to consider, [iso1 iso2]
%
% Optional inputs: (for use with varargin, make a cell array of strings followed by the argument. See example above)
%   coeffs  | (arb)       | Coefficients from normalizing the wavefunctions. This can save on propagation time once you've gotten the coefficients once.
%   plot    | boolean     | flag for plotting
%   figNum  | integer     | if plotting, use this figure number to plot on (useful for combining)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
%--------------------------------------------------------------------------
%   Name    | Input Unit  | Description
%-----------|-------------|-------------------------------------------------
%   x       | Bohr radius | space vector where the value of the wavefunction was determined
%   V       | Hartree     | potential energy curve evaluated at V(x)
%   R       | (arb)       | radial wavefunction, normalization is arbitrary
%   C       | (arb)       | asymptotic solution coefficients, see discussion below for details
%           |             | 

%% Function setup
% Conversions
amu2Au = 1.82289e3;  % conversion from amu to atomic units of mass (see note above)
kel2Ha = 3.16683e-6; % conversion from Kelvin to Hartree

% Handy numbers which probably don't need to be messed with but better to group them together
% than to have magic numbers down below
minGuessRange = [0 10]; % Range to expect the potential minimum in, used for setting step size
stepsPerOsc   = 12;     % Number of points per oscillation, more points = long run time. Solution may blow up if too small
scatRegion    = 150;
asymRegion    = 1e3;
asymPntSkip   = 100; % In the asymptotic region, only plot every nth point (speeds up plot manipulation)

% Specify default behavior (can be modified by varargin below)
calcNorm = 1; % Asymptotic solution will be used to normalize the numerical solution (typically must propagate out to ~15e3 - 20e3 Bohr to ensure proper matching)
flagPlot = 0; % Default to no plotting
flagFig  = 0;

% Loop through specified optional parameters (if any)
for i = 1:2:length(varargin)
    switch varargin{i}
        case 'coeffs'
            C        = varargin{i+1};
            calcNorm = 0;
        case 'plot'
            flagPlot = varargin{i+1};
        case 'figNum'
            flagFig = 1;
            figNum  = varargin{i+1};
    end
end


%% Function definitions
% Types of potentials
V_centr = @(el, m, r) el*(el + 1)./(2*m*r.^2);
funcV   = @(el, m, r) srPEC_1S0plus1S0(r)' + V_centr(el, m, r);

%Wavevector definition
k       = @(mu, E) sqrt(2*mu*E);

% Analytic aymptotic solution
posPlane = @(A, k, r) A*exp(1i*k.*r);
negPlane = @(B, k, r) B*exp(-1i*k.*r);
%fRAsym   = @(x, k, C, el) [ (-1i)^(el+1)*posPlane(C(1), k, x) (1i)^(el+1)*negPlane(C(2), k, x) ];
funAsymR  = @(r, k, C, el) [ (C(1)/2)*besselh(el, 1, k.*r) (C(2)/2)*besselh(el, 2, k.*r) ]; % Sakurai 2nd Ed. - Eq. 6.4.49

% Analytic aymptotic solution (different approach)
inTrig = @(z, k, el) k.*z - el*pi/2;
yAsym  = @(z, k, el, delL) (sin(inTrig(z, k, el)) + tan(delL)*cos(inTrig(z, k, el)));
fDelL  = @(z, k, el, numY) atan((numY(2)*sin(inTrig(z(1), k, el)) - numY(1)*sin(inTrig(z(2), k, el)))/...
                                (numY(1)*cos(inTrig(z(2), k, el)) - numY(2)*cos(inTrig(z(1), k, el))));

% General Numerov propagator (used for linear second order differential equations)
% INPUTS
%   d: grid spacing (grid must be equidistanct)
%   k: 3 element squared wavevector, [k_n-2, k_n-1, k_n] where k_i is the value of the local wavevector at position x_i
%   R: 2 element value vector,       [R_n-2, R_n-1]      where R_i is the value of the radial wavefunction at position x_i
numerov = @(d, k, R) ((2 - (10*d^2/12)*k(2))*R(2) - (1 + (d^2/12)*k(1))*R(1))/(1 + (d^2/12)*k(3));

%% Setup
% Physical specifics
mu   = prod(masses)/sum(masses)*amu2Au;
Einc = Tinc*kel2Ha;

% guess for range of minimum for determining the grid spacing
zFindMin  = linspace(minGuessRange(1), minGuessRange(2), 1e3);
lambdaMax = 2*pi/k(mu, Einc - min(funcV(el, mu, zFindMin))); % fastest oscillation is at deepest point in potential
d         = lambdaMax/stepsPerOsc; % At least 10 points per oscillation

% Grid spacing (must be equidistant), z - the interparticle spacing
r    = rRange(1):d:rRange(2);

% Get vectors of potential and local wavevector over the space of interest
V     = funcV(el, mu, r);        % potential energy curve
kLoc  = k(mu, Einc - V); % local wavevector
kAsym = k(mu, Einc);     % asymptotic wavevector
    
% Initialize solution vector
numR    = deal(nan(1,length(r)));
numR(1) = eps; %arbitrary (just needs to be really small)
numR(2) = ((1 + k(mu, V(1) - Einc)*d)*numR(1)); % assume linear first step

%% Solver
% Propagate radial function from within classically forbidden region
for i = 3:length(r)
    % This is where the magic happens
    numR(i) = numerov(d, kLoc([i-2 i-1 i]).^2, numR([i-2 i-1]));
end

if calcNorm
    % See Sakurai section 6.4 for discussion (particularly eq. 6.4.49 & .51)
    % See also http://www.physics.udel.edu/~msafrono/425-2011/Lecture%2023.pdf for discussion
    % and eq. 42-49 of http://www2.chem.umd.edu/groups/alexander/teaching/inelastic_scattering.pdf
    
    % Find the value at the asymptotic solution to get normalization by enforcing (these are matricies)
    % [ h_{el}^(1)(k*r_n)      h_{el}^(2)(k*r_n)     ] * [ A ] = [ R(r_n)/r_n        ]
    % [ h_{el}^(1)(k*r_{n-1})  h_{el}^(2)(k*r_{n-1}) ]   [ B ]   [ R(r_{n-1}/r_{n-1} ]
    C = funAsymR(r([end end-1])', kAsym, [1 1], el) \ (numR([end end-1])./r([end end-1]))';
    C = C./norm(C); % overall magnitude shouldn't matter, so convenient to work around unity
end

% Find the phase shift
S_el    = C(1)/C(2);
f_el    = (S_el - 1)./(2i*kAsym);
del_el  = real(log(S_el)/(2i));

% Get the normalized wavefunctions
asymR = sum(real(funAsymR(r', kAsym, C, el)), 2)'; 
numR  = asymR(end)/numR(end)*numR;

if flagPlot;
    scatPart = r < scatRegion;
    asymPart = r > asymRegion;
    
    if flagFig; figure(figNum); else figure; end
    axHan = subplot(2,1,1); hold on
    plot(r(scatPart), V(scatPart) , r(scatPart), numR(scatPart))
    xlim([0 scatRegion]);
    ylim([min(numR(scatPart)) max(numR(scatPart))]*1.1)
    xlabel('Interparticle distance, r [$a_0$]', 'Interpreter' , 'latex')
    ylabel('Wavefunction (arb.)', 'Interpreter' , 'latex');
    box on
    axHan.FontSize = 22;
    axHan.YTick    = [];
    axHan.TickLabelInterpreter = 'latex';
    legend({'V(r)', sprintf('$%g\\,+\\,%g, l = %g, \\frac{E_{inc}}{k_B} = %g\\times10^{%g}$ K',...
                            masses, el, Tinc/10^(floor(log10(Tinc))), floor(log10(Tinc)))},...
        'Interpreter' , 'latex' )
    
    
    axHan = subplot(2,1,2); hold on
    asymInd = length(asymPart(asymPart == 0)); % index of the position cutoff
    plot(r(asymInd:asymPntSkip:end)./1e3, V(asymInd:asymPntSkip:end) ,...
         r(asymInd:asymPntSkip:end)./1e3, numR(asymInd:asymPntSkip:end) ,...
         r(asymInd:asymPntSkip:end)./1e3, asymR(asymInd:asymPntSkip:end))
    xlim([asymRegion r(end)]./1e3)
    xlabel('Interparticle distance, r [$a_0\times10^3$]', 'Interpreter' , 'latex')
    ylabel('Asymptotic behavior (arb.)', 'Interpreter' , 'latex')
    box on
    axHan.FontSize = 22;
    axHan.YTick    = [];
    axHan.TickLabelInterpreter = 'latex';
    legend({'V(r)', 'Numerical wavefunction', 'Asymp solution:\\$\frac{R(r)}{r} = C\left[S_{l}h_{l}^{1}(kr) + h_{l}^{2}(kr)\right]$'},...
        'Interpreter' , 'latex' )
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Try another appraoch to normalizing. Useful for debugging
% This uses another form of the asymptotic solution, primarily to check that the
% calculated phase shifts are the same.
% del_el = fDelL(x([end end-1]), kAsym, el, R([end end-1]))
% A    = yAsym(z(end), kAsym, el, del_el)/R(end);
% R2   = A*R;
%      z, 0.5*Rasym,...
%      z, 0.5*sum(real(fRAsym(z', k(Einc), [1 1])), 2),...

%% Output
% Build output table row (can be used to combine with other outputs)
solOut = table(masses, Tinc, el, r, V, numR, asymR, del_el, S_el, f_el, {C});

%solOut.Properties.VariableNames        = {'masses', 'T_inc', 'el', 'r', 'PEC', 'numR', 'asymR', 'del_el', 'S_el', 'coeffs'};
solOut.Properties.VariableDescriptions = {'input masses', 'input scattering energy', 'partial wave', 'interparticle spacing',...
    '1S0 potential energy curve', 'numerical solution to R(r)', 'asymptotic solution to R(r)', 'phase shift for el partial wave',...
    'single channel S matrix', 'partial wave amplitude for el', 'asymptotic coefficients'};
solOut.Properties.VariableUnits        = {'a.m.u.', 'Kelvin', 'integer', 'a_0', 'Ha', 'arb', 'arb', 'radians', 'no dim', 'no dim', 'no dim'};
end