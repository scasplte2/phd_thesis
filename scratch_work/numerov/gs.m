function gs
%% Problem specifics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sr-Sr scattering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zStart = 0.1; 
zEnd   = 1000; %bohr
mu     = 84/2; % molecular mass

% incoming energy
%Tinc = 3000;           % Kelvin
Einc = 0.001; % Hartree

% partial wave
el = 0;

% Potentials coeffs given in Ha a0 (Hartree Bohr)
% For use with V_modLJ
% C6  = 3164; 
% C8  = 3.82e5;
% C10 = 5.05e7;
% C12 = 7.464e9;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ar-Ar scattering example from 
% http://www2.chem.umd.edu/groups/alexander/teaching/inelastic_scattering.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% zStart = 5.5; 
% zEnd   = 1000;  % bohr
% 
% mu   = 19.974; % molecular mass
% Einc = 5e-2; % incoming energy
% el   = 0;    % partial wave
% 
% %Potentials coeffs given in Ha a0 (Hartree Bohr)
% % for use with V_morse
% De   = 4.456;
% re   = 7.1053;
% beta = 0.89305;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Trial3 problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% zStart = 0.01;
% zEnd   = 100;
% 
% mu     = 0.5; % mass
% Einc   = 1;
% el     = 0;
% 
% % Potential coeffs
% % for use with V_LJ
% A      = 500;
% C6     = 1;
% C12    = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Function definitions
% Types of potentials
V_morse = @(r) De*(exp(-2*beta*(r-re)) - 2*exp(-beta*(r-re)));% Morse potential
V_LJ    = @(r) A*(C12./r.^12 - C6./r.^6) + el*(el + 1)./r.^2;
V_modLJ = @(r) C12./r.^12 - C6./r.^6 - C8./r.^8 - C10./r.^10; % modified Lennard - Jones potential

%Wavevector definition
k       = @(E) sqrt(2*mu*E);

% Analytic aymptotic solution
posPlane = @(A, k, x) A*exp(1i*k.*x);
negPlane = @(B, k, x) B*exp(-1i*k.*x);
fRAsym   = @(x, k, C) [ (-1i)^(el+1)*posPlane(C(1), k, x) (1i)^(el+1)*negPlane(C(2), k, x) ];

% Analytic aymptotic solution (different approach)
inTrig = @(z, k, el) k.*z - el*pi/2;
yAsym  = @(z, k, el, delL) (sin(inTrig(z, k, el)) + tan(delL)*cos(inTrig(z, k, el)));
fDelL  = @(z, k, el, numY) atan((numY(2)*sin(inTrig(z(1), k, el)) - numY(1)*sin(inTrig(z(2), k, el)))/...
                                (numY(1)*cos(inTrig(z(2), k, el)) - numY(2)*cos(inTrig(z(1), k, el))));

% General Numerov propagator (used for linear second order differential equations)
% INPUTS
%   d: grid spacing (grid must be equidistanct)
%   k: 3 element squared wavevector, [k_n-2, k_n-1, k_n] where k_i is the position at a specified index i, and i = n is the current index
%   z: 3 element space vector, [z_n-2, z_n-1, z_n] where z_i is the position at a specified index i
%   R: 2 element value vector, [R_n-2, R_n-1] where R(z) is the scaled radial wavefunction being solved for
numerov = @(d, k, z, R) ((2 - (10*d^2/12)*k(2))*R(2) - (1 + (d^2/12)*k(1))*R(1))/(1 + (d^2/12)*k(3));

%% Setup
% Choose the type of potential to be solved
funcV = @(r) srGSPot(r);

% guess for range of minimum for determining the grid spacing
zFindMin  = linspace(0, 10, 1e3);
lambdaMax = 2*pi/k(Einc - min(funcV(zFindMin))); % fastest oscillation is at deepest point in potential
d         = lambdaMax/25; % At least 10 points per oscillation

% Grid spacing (must be equidistant), z - the interparticle spacing
z    = zStart:d:zEnd;

% Get vectors of potential and local wavevector over the space of interest
V    = funcV(z);    % potential energy curve
kLoc = k(Einc - V); % local wavevector
    
% Initialize solution vector
R    = deal(nan(1,length(z)));
R(1) = eps; %arbitrary
R(2) = ((1 + k(V(1) - Einc)*d)*R(1)); % assume linear first step

%% Solver
% Propagate radial function from within classically forbidden region
for i = 3:length(z)
    R(i) = numerov(d, kLoc([i-2 i-1 i]).^2, z([i-2 i-1 i]), R([i-2 i-1]));
end

% Find the value at the asymptotic solution to get normalization
C = fRAsym(z([end end-1])', k(Einc), [1 1])\R([end end-1])';
C = C./norm(C);

% Find the phase shift
S_el = C(1)*k(Einc)/C(2);
delL = real(log(S_el)/(2*1i));

% Get the normalized wavefunctions
Rasym = sum(real(fRAsym(z', k(Einc), C)), 2); 
R1    = Rasym(end)/R(end)*R;

% Plot scale
ampV      = -min(V);
pEnScale  = Einc/ampV;
pAmpScale = 0.5/max(R1((end-round((length(z)/4))):end));

figure;
subplot(2,1,1)
plot(z, V./ampV ,...
     z, pAmpScale*R1)
ylim([-1.1 1.1])
xlim([0 50])
%      z, 0.5*Rasym,...
%      z, 0.5*sum(real(fRAsym(z', k(Einc), [1 1])), 2),...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Try other appraoch to normalizing
delL = fDelL(z([end end-1]), k(Einc), el, R([end end-1]));
A    = yAsym(z(end), k(Einc), el, delL)/R(end);
R2   = A*R;

% Plot scale
ampV      = -min(V);
pEnScale  = Einc/ampV;
pAmpScale = abs(0.5/max(R2((end-round((length(z)/4))):end)));

subplot(2,1,2)
plot(z, V./ampV ,...
     z, pAmpScale*R2 + pEnScale)
ylim([-1.1 max([pEnScale*1.5 1.1])])
xlim([0 50])
%      z, 0.5*Rasym,...
%      z, 0.5*sum(real(fRAsym(z', k(Einc), [1 1])), 2),...
end