function trial
%% Problem spceifics
% Ar-Ar scattering example from 
% http://www2.chem.umd.edu/groups/alexander/teaching/inelastic_scattering.pdf
zStart = 5.5; 
zEnd   = 25; %bohr
mu = 19.974; % molecular mass

% incoming energy
Einc = 5e-4;

% Potentials coeffs given in Ha a0 (Hartree Bohr)
De = 4.456e-4;
re = 7.1053;
beta = 0.89305;

% incoming energy
Einc = 5e-4;

%% Function definitions
V_morse = @(r) De*(exp(-2*beta*(r-re)) - 2*exp(-beta*(r-re)));% Morse potential

%Wavevector definition
k       = @(E) sqrt(2*mu*E);

% Spherical bessel functions
jRicBes = @(l, x) sqrt(pi*x/2).*besselj(l+(1/2),x);
yRicBes = @(l, x) -sqrt(pi*x/2).*bessely(l+(1/2),x);

% Analytic aymptotic solution
u    = @(z, k, el, coeffs) (1./k).*(coeffs(1).*jRicBes(el, k.*z) + coeffs(2).*yRicBes(el, k.*z));
uCol = @(z, k, el, coeffs) [ (1./k).*coeffs(1).*jRicBes(el, k.*z)
                             (1./k).*coeffs(2).*yRicBes(el, k.*z) ];

% General Numerov propagator (used for linear second order differential equations)
% INPUTS
%   d: grid spacing (grid must be equidistanct)
%   k: 3 element squared wavevector, [k_n-2, k_n-1, k_n] where k_i is the position at a specified index i, and i = n is the current index
%   R: 2 element value vector, [R_n-2, R_n-1] where R(z) is the scaled radial wavefunction being solved for
numerov  = @(d, k, R) ((2 - (10*d^2/12)*k(2))*R(2) - (1 + (d^2/12)*k(1))*R(1))/(1 + (d^2/12)*k(3));

%% Setup
% Choose the type of potential
V    = @(r) V_morse(r);
kLoc = @(r) k(Einc - V(r)); % local wavevector

% Grid spacing (must be equidistant)
d = 0.02;

% Partial wave
el = 0;
    
% z - the interparticle spacing
z = zStart:d:zEnd;
    
% Initialize solution vector
R    = deal(nan(1,length(z)));
R(1) = 1e-10/z(1); %arbitrary
R(2) = ((1 + k(V(z(1)) - Einc)*d)*R(1))/z(2); % assume linear first step

%% Solver
% Propagate radial function from within classically forbidden region
for i = 3:length(z)
    R(i) = numerov(d, kLoc(z([i-2 i-1 i])).^2, R([i-2 i-1]));
end

% Find the value at the asymptotic solution to get normalization
x = k(Einc).*z([end end-1]);
asymSol = [ jRicBes(el, x(1)), yRicBes(el, x(1)) 
            jRicBes(el, x(2)), yRicBes(el, x(2)) ];
            
numSol  = [ k(Einc)*R(end)
            k(Einc)*R(end-1) ];
        
% Solve for incoming and outgoing wave coefficients
C  = asymSol\numSol;

% Normalize coefficients
C    = C./norm(C);
Scat = C(1) - abs(C(2))^2/conj(C(1));

% Find the normalization for the numerical solution
A = u(z(end), kLoc(z(end)), el, C)/R(end);
uPsi = @(z) u(z, kLoc(z), el, C);

psiR = A*R.*z;
psiAsym = @(z) u(z, kLoc(z), el, coeffs).*z;

scalePsiR = psiR*1e-5;
scalePsi  = @(z) psiAsym(z)*1e-5;


figure;
plot(z, V(z), z, A*R*1e-3)
ylim([min(V(z)) max(A*R*1e-3)]*1.1)
xlim([5 40])

% figure;
% plot(z, V(z), z, scalePsiR, z(z>15), scalePsi(z(z>15)))
% ylim([min(V(z)) max(scalePsiR)]*1.1)
% %ylim([-0.5 0.5]*1e-3); xlim([10 200]);
% 
% figure;
% plot(z, V(z), z, Psi)
% ylim([-0.5 0.5]*1e-3); xlim([10 200]);
end