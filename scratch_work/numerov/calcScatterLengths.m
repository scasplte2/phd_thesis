function aTab = calcScatterLengths
% Function to calculate the scattering lengths by numerically estiamting the phase accumulation
% due to the interaction potential between two atoms. All that is needed is a form of the 
% potential that can be evaluated at any point r from the inner turning point to infinity.
% R_inner is defined as the point where the potentail crosses E=0.

%% Setup constants
aBohr  = 0.52917721067; % angstroms
hBar   = 1;

mConv   = 1822.889;   % atomic mass units to atomic units
enConv  = 4.55634e-6; % wavenumber to Hartree
lenConv = 1/aBohr;    % angstroms to bohr radii

funcPEC = @(x) srPEC_1S0plus1S0(x);

% Given in Tiemann specified units: cm-1 and Angstroms
% then converted to appropriate system for computation
C6 = 1.52701e7*enConv*lenConv^6; %using feely varying C6 to be consistent with potential used
%C6 = 1.525e7*enConv*lenConv^6; 
Ri = 3.963*lenConv; 

m84 = 83.913425;
m86 = 85.9092607309;
m87 = 86.908877497;
m88 = 87.9056122571;
srM = [m84 m86 m87 m88].*mConv;

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

%% Calculate scattering lengths
aMat = nan(length(srM));
for i = 1:length(srM);
    for j = i:length(srM);
        mu = prod([srM(i) srM(j)]) / sum([srM(i) srM(j)]);
        aMat(i,j) = a(mu);
    end
end

% Save into table so that read out is easier
aTab = table(aMat(:,1), aMat(:,2), aMat(:,3), aMat(:,4));
aTab.Properties.VariableNames = {'sr84' 'sr86' 'sr87' 'sr88'};
aTab.Properties.RowNames      = {'sr84' 'sr86' 'sr87' 'sr88'};
end