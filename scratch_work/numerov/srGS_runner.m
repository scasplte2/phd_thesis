% Use atomic units throughout 

%% Physical definitions
amu2Au = 1822.889;     % conversion from amu to atomic units of mass (see note above)
kel2Ha = 3.16683e-6;   % conversion from Kelvin to Hartree
cm2Ha  = 4.55634e-6;
aBohr  = 0.52917721067;

% Strontium masses
% should use these for accurate measurements, even these slight differences can
% result in a factor of 2 difference in estimated scattering length
m84 = 83.913425;
m86 = 85.9092607309;
m87 = 86.908877497;
m88 = 87.9056122571;

%% Function definitions
k        = @(mu, E) sqrt(2*mu*E);
dBLambda = @(mu, E) 2*pi/k(mu, E);

% Useful potential forms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%C = [ 3.164315349415183e+03 3.753373129228250e+05 5.226574577253875e+07 7.435*1e9  ]; % Large negative (84+84)
%C = [ 3.164315349415183e+03 3.753373129228250e+05 5.226574577253875e+07 7.4345*1e9 ]; % Large positive (84+84)
%mod_LJ = @(r) C(4)./r.^12 - C(1)./r.^6 - C(2)./r.^8 - C(3)./r.^10;

C = [948*cm2Ha 1/aBohr];
LJ = @(r) C(1)*( (C(2)./r).^12 - 2*(C(2)./r).^6 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Problem specifics
T    = 10e-9; 
mass = [m86 m86];

maxRad = 2*dBLambda(prod(mass)/sum(mass)*amu2Au, T*kel2Ha); 
range1 = maxRad/1000; % 1500 for 86 real potential, default to 150 for 84
range2 = maxRad/10;

%% Evaluate wavefunction
tic; 
%t = sr_1S0_scatter(T, 0, [5 maxRad], mass, 'plot', 1, 'region', [range1 range2]); 
%t = sr_1S0_scatter(T, 0, [5 maxRad], mass, 'plot', 1, 'region', [range1 range2], 'funcPEC', mod_LJ);
t = sr_1S0_scatter(T, 0, [0.1 maxRad], mass, 'plot', 1, 'region', [range1 range2], 'funcPEC', LJ); 
toc

%% Estimate scattering length
ov      = t.outVecs;
rLim    = length(nonzeros(ov.r < range1));
fitInds = [rLim rLim+1];

P = polyfit(ov.r(fitInds), ov.numR(fitInds), 1);
subplot(2,1,1); hold on
plot(ov.r(1:rLim), polyval(P, ov.r(1:rLim)))

tmp = polyfit(ov.numR(fitInds), ov.r(fitInds), 1);
sprintf('Estimated scattering length: %g', tmp(2))