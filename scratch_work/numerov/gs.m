function gs
%% Problem spceifics
zStart = 3; 
zEnd   = 1000; %bohr
mu = 84/2; % molecular mass

% incoming energy
Tinc = 200e-9;           % Kelvin
Einc = Tinc/3.1668e-6; % Hartree

% Potentials coeffs given in Ha a0 (Hartree Bohr)
C6  = 3164; 
C8  = 3.82e5;
C10 = 5.05e7;
C12 = 7.464e9;

%% Function definitions
% Types of potentials
V_modLJ = @(r) C12./r.^12 - C6./r.^6 - C8./r.^8 - C10./r.^10; % modified Lennard - Jones potential
V_morse = @(r) De*(exp(-2*beta*(r-re)) - 2*exp(-beta*(r-re)));% Morse potential

%Wavevector definition
k       = @(E) sqrt(2*mu*E);

% Spherical bessel functions
jRicBes = @(l, x) x.*sqrt((pi./(2.*x))).*besselj(l+(1/2),x);
yRicBes = @(l, x) x.*sqrt((pi./(2.*x))).*bessely(l+(1/2),x);

% Analytic aymptotic solution
u = @(z, k, el, coeffs) (1./k).*(coeffs(1).*jRicBes(el, k.*z) + coeffs(2).*yRicBes(el, k.*z));

% General Numerov propagator (used for linear second order differential equations)
% INPUTS
%   d: grid spacing (grid must be equidistanct)
%   V: the interaction potential between particles
%   z: 3 element space vector, [z_n-2, z_n-1, z_n] where z_i is the position at a specified index i, and i = n is the current index
%   R: 2 element value vector, [R_n-2, R_n-1] where R(z) is the scaled radial wavefunction being solved for
R_n  = @(d, k, z, R) (-1*(1 + 1/12*(d*k(z(1)))^2)*R(1) + (2 - 10/12*(d*k(z(2)))^2)*R(2))/(1 + 1/12.*(d*k(z(3))).^2);

%% Setup
% Choose the type of potential
V    = @(r) V_modLJ(r);
kloc = @(r) k(Einc - V(r)); % local wavevector

% Grid spacing (must be equidistant)
rMin = 9.688; % Position of minimum defined above
lambdaMax = 1/(2*pi*kloc(rMin)); % fastest oscillation is at deepest point in potential
d = lambdaMax/10; % At least 10 points per oscillation

% Partial wave
el = 0;
    
% z - the interparticle spacing
z = zStart:d:zEnd;
    
% Initialize solution vector
R = zeros(1,length(z));
R(1) = eps; %arbitrary
R(2) = (1 + k(V(z(1)) - Einc)*d)*R(1); % assume linear first step

%% Solver
% Propagate radial function from within classically forbidden region
for i = 3:length(z)
    R(i) = R_n(d, kloc, z([i-2 i-1 i]), R([i-2 i-1]));
end

% Find the value at the asymptotic solution to get normalization
x = kloc(z([end end-1])).*z([end end-1]);
asymSol = [ jRicBes(el, x(1)) yRicBes(el, x(1)) 
            jRicBes(el, x(2)) yRicBes(el, x(2)) ];
            
numSol  = [ kloc(z(end))*R(end)
            kloc(z(end-1))*R(end-1) ];
        
% Solve for incoming and outgoing wave coefficients
coeffs  = asymSol\numSol;

% Normalize coefficients
coeffs  = coeffs./norm(coeffs);

% Find the normalization for the numerical solution
A = u(z(end), k(Einc - V(z(end))), 0, coeffs)/R(end);

psiR = A*R./z;
psi = @(z) u(z, kloc(z), el, coeffs)./z;

scalePsi = psiR*0.1;
plot(z, V(z), z, scalePsi)
ylim([min([scalePsi V(z)]) max(scalePsi)]*1.1)
xlim([0 100])

%coVec(:,j) = coeffs;
%end

%coeffs
%R(end)/z(end)


%plot(z, u(z, k(Einc), el, coeffs)./z)
%plot(z, V(z), z, 100*R./z)

% Find minimum grid spacing
% rMin = 9.688; % Position of minimum defined above
% lambdaMax = 1/(2*pi*kloc(rMin)); % fastest oscillation is at deepest point in potential
% d = lambdaMax/10; % At least 10 points per oscillation

% numSteps = floor((rEnd-rStart)./d);
% R = linspace(rStart,rEnd,numSteps);
% % % Initial values of psi and r (m1 - minus 1, m2 - minus 2)
% psi_m2 = 1e-10;     r_m2 = rMin - 2*d;
% psi_m1 = (1 + kloc(r_m2);  r_m1 = rMin - d;
% 
% psi = zeros(1,length(R));
% for i = 1:length(R)
%     % Calculate psi at current point (based off last two points, m1 and m2)
%     psi(i) = 1/(d^2*g(R(i))) * (psi_m2*(12-d^2*g(r_m2)) - 2*psi_m1*(5*d^2*g(r_m1) + 12));
%     
%     % Shift m1 and m2 values for next iteration
%     psi_m2 = psi_m1;  r_m2 = R(i) - 2*d;
%     psi_m1 = psi(i);  r_m1 = R(i) - d;
% end
end