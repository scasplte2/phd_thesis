function trial2
%% Problem spceifics
% Ar-Ar scattering example from 
% http://www2.chem.umd.edu/groups/alexander/teaching/inelastic_scattering.pdf
% zStart = 1; 
% zEnd   = 100;  % bohr
% d      = 0.01; % step size
% 
% mu   = 19.974; % molecular mass
% Einc = 5e-4; % incoming energy
% el   = 0;    % partial wave

% Potentials coeffs given in Ha a0 (Hartree Bohr)
% De = 4.456e-4;
% re = 7.1053;
% beta = 0.89305;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trial3 problem
zStart = 0.01;
zEnd   = 100;
d      = 0.01;

mu     = 0.5; % mass
Einc   = 1;
el     = 0;

% Potential coeffs
A      = 500;
C6     = 1;
C12    = 1;


%% Function definitions
% Types of potentials
V_morse = @(r) De*(exp(-2*beta*(r-re)) - 2*exp(-beta*(r-re)));% Morse potential
V_LJ    = @(r) A*(C12./r.^12 - C6./r.^6) + el*(el + 1)./r.^2;

%Wavevector definition
k       = @(E) sqrt(2*mu*E);

% Analytic aymptotic solution
posPlane = @(A, k, x) A*exp(1i*k.*x);
negPlane = @(B, k, x) B*exp(-1i*k.*x);
fRAsym   = @(x, k, C) [ (-1i)^(el+1)*posPlane(C(1), k, x) (1i)^(el+1)*negPlane(C(2), k, x) ];

% General Numerov propagator (used for linear second order differential equations)
% INPUTS
%   d: grid spacing (grid must be equidistanct)
%   k: 3 element squared wavevector, [k_n-2, k_n-1, k_n] where k_i is the position at a specified index i, and i = n is the current index
%   z: 3 element space vector, [z_n-2, z_n-1, z_n] where z_i is the position at a specified index i
%   R: 2 element value vector, [R_n-2, R_n-1] where R(z) is the scaled radial wavefunction being solved for
numerov = @(d, k, z, R) ((2 - (10*d^2/12)*k(2))*R(2) - (1 + (d^2/12)*k(1))*R(1))/(1 + (d^2/12)*k(3));

%% Setup
% Choose the type of potential
V    = @(r) V_LJ(r);
kloc = @(r) k(Einc - V(r)); % local wavevector

% z - the interparticle spacing
z = zStart:d:zEnd;
    
% Initialize solution vector
R    = deal(nan(1,length(z)));
R(1) = 1e-10; %arbitrary
R(2) = ((1 + k(V(z(1)) - Einc)*d)*R(1)); % assume linear first step

%% Solver
% Propagate radial function from within classically forbidden region
for i = 3:length(z)
    R(i) = numerov(d, kloc(z([i-2 i-1 i])).^2, z([i-2 i-1 i]), R([i-2 i-1]));
end

% Find the value at the asymptotic solution to get normalization
C = fRAsym(z([end end-1])', k(Einc), [1 1])\R([end end-1])';
C = C./norm(C);

% Find the phase shift
S_el = C(1)*k(Einc)/C(2);
delL = real(log(S_el)/(2*1i));

% Get the normalized wavefunctions
Rasym = sum(real(fRAsym(z', k(Einc), C)), 2); 
R     = Rasym(end)/R(end)*R;

figure;
plot(z, V(z)./(-min(V(z))),...
     z, 0.5*Rasym,...
     z, 0.5*sum(real(fRAsym(z', k(Einc), [1 1])), 2),...
     z, 0.5*R)
ylim([-1.1 1.1])
xlim([0 25])
end