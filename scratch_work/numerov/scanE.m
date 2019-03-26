function solTable = scanE

k        = @(mu, E) sqrt(2*mu*E);
dBLambda = @(mu, E) 2*pi/k(mu, E);

kel2Ha = 3.16683e-6; % conversion from Kelvin to Hartree
amu2Au = 1.82289e3;  % conversion from amu to atomic units of mass

% Input values
% T_inc  = 10.^(-1*(3:9)); % K
% masses = [84 84];
% el     = 4;

T_inc  = 100e-9; % K
masses = [86 86];
el     = 0;

% starting radius (haven't needed to tweak)
startRad = 5;

% plot output on same figure
%figNum = figure;

for i = 1:length(T_inc)
    maxRad = 2*dBLambda((prod(masses)/sum(masses))*amu2Au, T_inc(i)*kel2Ha);
    
    tic;
    out = sr_1S0_scatter(T_inc(i), el, [startRad maxRad], masses, 'region', [200 maxRad*0.1]);
    toc
    
    % I don't know how to preallocate tables in this version of Matlab.
    % Later version seem to do it, use varfun(@class,t,'OutputFormat','cell') to get varTypes
    if i == 1; solTable = out; else solTable = [solTable; out]; end
end
end