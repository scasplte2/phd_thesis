function trial3
%% Function definitions
V_LJ    = @(A, C6, C12, r, el) A*(C12./r.^12 - C6./r.^6) + el*(el + 1)./r.^2;

%Wavevector definition (with hBar = 1)
k       = @(m, E) sqrt(2*m*E);

% Analytic aymptotic solution
inTrig = @(z, k, el) k.*z - el*pi/2;
yAsym  = @(z, k, el, delL) (sin(inTrig(z, k, el)) + tan(delL)*cos(inTrig(z, k, el)));
fDelL  = @(z, k, el, numY) atan((numY(2)*sin(inTrig(z(1), k, el)) - numY(1)*sin(inTrig(z(2), k, el)))/...
                                (numY(1)*cos(inTrig(z(2), k, el)) - numY(2)*cos(inTrig(z(1), k, el))));

% General Numerov propagator (used for linear second order differential equations)
% INPUTS
%   d: grid spacing (grid must be equidistanct)
%   k: 3 element squared wavevector, [k_n-2, k_n-1, k_n] where k_i is the position at a specified index i, and i = n is the current index
%   R: 2 element value vector,       [R_n-2, R_n-1]      where R(z) is the scaled radial wavefunction being solved for
numerov  = @(d, k, R) ((2 - (10*d^2/12)*k(2))*R(2) - (1 + (d^2/12)*k(1))*R(1))/(1 + (d^2/12)*k(3));

%% Problem spceifics
% Region
zStart = 0.01;
zEnd   = 1000;

% mass
m      = 0.5;

%Potential coefficients
A      = 500;
C6     = 1;
C12    = 1;

%% Setup
% Partial wave
el = 0;

% Incoming energy
Einc = 1;

% Choose the type of potential
V    = @(r) V_LJ(A, C6, C12, r, el);
kLoc = @(r) k(m, Einc - V(r)); % local wavevector

d = 0.01;
    
% z - the interparticle spacing
z = zStart:d:zEnd;
    
% Initialize solution vector
yNum = nan(1,length(z));
yNum(1) = 1e-10; %arbitrary
yNum(2) = (1 + k(m, V(z(1)) - Einc)*d)*yNum(1); % assume linear first step

%% Solver
% Propagate radial function from within classically forbidden region
for i = 3:length(z)
    yNum(i) = numerov(d, kLoc(z([i-2 i-1 i])).^2, yNum([i-2 i-1]));
end

delL = fDelL(z([end end-1]), kLoc(z([end end-1])), el, yNum([end end-1]));
A    = yAsym(z(end), k(m, Einc), el, delL)/yNum(end);
yNum = A*yNum;

figure;
plot(z, V(z),...
     z, (-min(V(z))/(4*max(yAsym(z, k(m, Einc), el, 0))))*yAsym(z, k(m, Einc), el, 0) + Einc,...
     z, (-min(V(z))/(4*max(yAsym(z, k(m, Einc), el, delL))))*yAsym(z, k(m, Einc), el, delL) + Einc,...
     z, (-min(V(z))/(4*max(yNum)))*yNum + Einc)
ylim([min(V(z)) -min(V(z))/3]*1.1)
xlim([0 10])
end