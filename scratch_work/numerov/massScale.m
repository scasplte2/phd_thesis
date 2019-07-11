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
muVec = linspace(83.5, 91.5, 500)*amu2au/2;

aVec = nan(1,length(muVec));
for i = 1:length(muVec);
    aVec(i) = a(muVec(i));
end

srM = [ funcMu(m84, m84) funcMu(m84, m86) funcMu(m84, m87) funcMu(m84, m88) ...
        funcMu(m86, m86) funcMu(m86, m87) funcMu(m86, m88) ...
        funcMu(m87, m87) funcMu(m87, m88) ...
        funcMu(m88, m88) ]*amu2au;
homoMix   = [1 5 8 10];
heteroMix = [2 3 4 6 7 9];

srVec = nan(1,length(srM));
for i = 1:length(srM);
    srVec(i) = a(srM(i));
end

%% Plotting
% Use this mask to avoid plotting the asymptotic line
part1 = 2*muVec/amu2au > 83.5 & 2*muVec/amu2au < 85.7;
part2 = 2*muVec/amu2au > 85.85 & 2*muVec/amu2au < 88.5;
part3 = 2*muVec/amu2au > 88.6 & 2*muVec/amu2au < 90.5;

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Want to plot the x-axis in isotopic 
plHan = plot(2*muVec(part1)/amu2au, aVec(part1),...
             2*muVec(part2)/amu2au, aVec(part2),...
             2*muVec(part3)/amu2au, aVec(part3));
plHan(1).LineWidth = 1.5;
plHan(2).LineWidth = 1.5;
plHan(3).LineWidth = 1.5;

% Create plot
plot(2*srM(homoMix)/amu2au, srVec(homoMix),...
    'MarkerFaceColor',[0.929411768913269 0.694117665290833 0.125490203499794],...
    'MarkerSize',8,...
    'Marker','o',...
    'LineWidth',1,...
    'LineStyle','none',...
    'Color',[0 0 0]);

% Create plot
plot(2*srM(heteroMix)/amu2au, srVec(heteroMix),...
    'MarkerFaceColor',[0.850980401039124 0.325490206480026 0.0980392172932625],...
    'MarkerSize',8,...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 0]);

% Create xlabel
xlabel('Mass [amu]');

% Create ylabel
ylabel('s-wave scattering length [a_0]');

xlim(axes1,[83.5 90.5]);
ylim(axes1,[-250 2000]);
box(axes1,'on');

% Set the remaining axes properties
set(axes1,'FontSize',20,'LineStyleOrderIndex',72);
end