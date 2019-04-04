function massScale
% Modeled after calcScatterLength 
% Function to calculate the scattering lengths by numerically estiamting the phase accumulation
% due to the interaction potential between two atoms. All that is needed is a form of the 
% potential that can be evaluated at any point r from the inner turning point to infinity.
% R_inner is defined as the point where the potentail crosses E=0.

%% Setup constants
aBohr  = 0.52917721067; % angstroms
hBar   = 1;

amu2au   = 1822.889;   % atomic mass units to atomic units
cm2Ha    = 4.55634e-6; % wavenumber to Hartree
ang2bohr = 1/aBohr;    % angstroms to bohr radii

m84 = 83.913425;
m86 = 85.9092607309;
m87 = 86.908877497;
m88 = 87.9056122571;

% Free vary model value
C6 = 1.526593e7*cm2Ha*ang2bohr^6; 
Ri = 3.963*ang2bohr; 

funcPEC = @(x) srPEC_1S0plus1S0(x, 'C6', 1.526593e7);

% find PEC zero crossing (assume near R_i)
minOpt = optimset('TolX', 1e-7); % set Tol to ensure that the integral below is negative
r0     = fminbnd(@(x) abs(funcPEC(x)), Ri*0.99, Ri*1.01, minOpt);

% van der Waals length
Rvdw = @(mu) (1/2)*(2*mu*C6/hBar^2)^(1/4);

% potential phase contribution
phiD = @(mu) integral(@(x) sqrt(-2*mu*funcPEC(x)), r0, Inf, 'ArrayValued', 1, 'AbsTol', 1e-10, 'RelTol', 1e-10);

% mean scattering length
aBar = @(mu) 4*pi*Rvdw(mu)/gamma(1/4)^2;

% s-wave scattering length
a    = @(mu) aBar(mu)*(1-tan(pi/4)*tan(phiD(mu) - pi/8));

% calculation of reduced mass
funcMu = @(m1, m2) prod([m1 m2])/sum([m1 m2]);

%% Calculate scattering lengths
% Define mass scaled vector
muVec = linspace(83, 91, 500)*amu2au/2;

aVec = nan(length(muVec));
for i = 1:length(aVec);
    aVec(i) = a(muVec(i));
end


%% Plotting
srM = []

figure;
plot(muVec/amu2au, aVec);

end