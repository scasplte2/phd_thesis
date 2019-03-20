function trial2
%% Problem spceifics
% Ar-Ar scattering example from 
% http://www2.chem.umd.edu/groups/alexander/teaching/inelastic_scattering.pdf
zStart = 1; 
zEnd   = 100; %bohr
mu     = 19.974; % molecular mass

% incoming energy
Einc = 5e-4;

% Potentials coeffs given in Ha a0 (Hartree Bohr)
De = 4.456e-4;
re = 7.1053;
beta = 0.89305;

%% Function definitions
% Types of potentials
V_morse = @(r) De*(exp(-2*beta*(r-re)) - 2*exp(-beta*(r-re)));% Morse potential

%Wavevector definition
k       = @(E) sqrt(2*mu*E);

% Analytic aymptotic solution
psiInf = @(z, k, C) C*exp(-1i*k.*z);

% General Numerov propagator (used for linear second order differential equations)
% INPUTS
%   d: grid spacing (grid must be equidistanct)
%   k: 3 element squared wavevector, [k_n-2, k_n-1, k_n] where k_i is the position at a specified index i, and i = n is the current index
%   z: 3 element space vector, [z_n-2, z_n-1, z_n] where z_i is the position at a specified index i
%   R: 2 element value vector, [R_n-2, R_n-1] where R(z) is the scaled radial wavefunction being solved for
numerov = @(d, k, z, R) ((2 - (10*d^2/12)*k(2))*R(2) - (1 + (d^2/12)*k(1))*R(1))/(1 + (d^2/12)*k(3));

%% Setup
% Choose the type of potential
V    = @(r) V_morse(r);
kloc = @(r) k(Einc - V(r)); % local wavevector

% Grid spacing (must be equidistant)
% rMin = 9.688; % Position of minimum defined above
% lambdaMax = 1/(2*pi*kloc(rMin)); % fastest oscillation is at deepest point in potential
% d = lambdaMax/10; % At least 10 points per oscillation

% Grid spacing (must be equidistant)
% rMin = 7.113; % Position of minimum defined above
% lambdaMax = 1/(2*pi*kloc(rMin)); % fastest oscillation is at deepest point in potential
% d = lambdaMax/50; % At least 10 points per oscillation
d = 0.01;
    
% z - the interparticle spacing
z = zEnd:-d:zStart;
    
% Initialize solution vector
psiT    = nan(1,length(z));
psiT(1) = psiInf(z(1), kloc(z(1)), 1);
psiT(2) = (1 + kloc(z(1)))*z(1);

%% Solver
% Propagate radial function from within classically forbidden region
for i = 3:length(z)
    psiT(i) = numerov(d, kloc(z([i-2 i-1 i])).^2, z([i-2 i-1 i]), psiT([i-2 i-1]));
end

% Find the value at the asymptotic solution to get normalization
x = k(Einc).*z([end end-1]);
asymSol = [ jRicBes(el, x(1)), yRicBes(el, x(1)); 
            jRicBes(el, x(2)), yRicBes(el, x(2)) ];
            
numSol  = [ k(Einc)*R(end)
            k(Einc)*R(end-1) ];
        
% Solve for incoming and outgoing wave coefficients
rcoeffs  = asymSol\numSol;

% Normalize coefficients
coeffs  = rcoeffs./norm(rcoeffs);

% Find the normalization for the numerical solution
A = u(z(end), k(Einc), el, coeffs)/R(end);

psiR = A*R./z;
psi = @(z) u(z, kloc(z), el, coeffs)./z;

scalePsi  = @(z) psi(z)*1e-3;
scalePsiR = psiR*1e-3;

plot(z, V(z), z, scalePsiR, z(z>15), scalePsi(z(z>15)))
ylim([min(V(z)) max(scalePsiR)]*1.1)

figure; plot(z, psiT)



% Last thing
end