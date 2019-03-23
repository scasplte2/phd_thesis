function aMat = calcScatterLengths
%% Specify units
% default to atomic units but can specify 'cgs' to use their values
units = 'au';

aBohr = 0.52917721067; % angstroms

%% Setup constants
switch units
    case 'cgs'
        h       = 6.626070150e-11; % gram*ang^2/s
        c       = 299792458e10;    % ang/s
        hBar    = h/(2*pi);
        
        mConv   = 1.66054e-24; % atomic mass unit to grams
        enConv  = h*c*1e-8;    % wavenumber to erg (cgs energy unit)
        lenConv = 1;           % keep distances in angstroms
        
        % Call PEC and specify cgs unit
        funcPEC = @(x) srPEC_1S0plus1S0(x, 'cgs');
    otherwise
        hBar   = 1;
        
        mConv   = 1.82289e3;  % atomic mass units to atomic units
        enConv  = 4.55634e-6; % wavenumber to Hartree
        lenConv = 1/aBohr;    % angstroms to bohr radii
        
        funcPEC = @(x) srPEC_1S0plus1S0(x);
end

% Given in Tiemann specified units: cm-1 and Angstroms
% then converted to appropriate system for computation
C6 = 1.525e7*enConv*lenConv^6;
Ri = 3.963*lenConv; 

m84 = 83.913425;
m86 = 85.9092607309;
m87 = 86.908877497;
m88 = 87.9056122571;
srM = [m84 m86 m87 m88].*mConv;

% find PEC zero crossing (assume near R_i)
minOpt = optimset('TolX', 1e-6); % set Tol to ensure that the integral below is negative
r0    = fminbnd(@(x) abs(funcPEC(x)), Ri*0.99, Ri*1.01, minOpt) + 1e-6;

% van der Waals length
Rvdw = @(mu) (1/2)*(2*mu*C6/hBar^2)^(1/4);

% potential phase contribution
phiD = @(mu) integral(@(x) sqrt(-2*mu*funcPEC(x)), r0, Inf, 'ArrayValued', 1, 'AbsTol', 1e-8, 'RelTol', 1e-10);

% mean scattering length
aBar = @(mu) 4*pi*Rvdw(mu)/gamma(1/4)^2;

% s-wave scattering length
a    = @(mu) aBar(mu)*(1-tan(pi/4)*tan(phiD(mu) - pi/8));

%% Calculate scattering lengths
aMat = zeros(length(srM));
for i = 1:length(srM);
    for j = i:length(srM);
        mu = prod([srM(i) srM(j)]) / sum([srM(i) srM(j)]);
        aMat(i,j) = a(mu);
    end
end

%always give output in Bohr (just to confuse you more)
%   ^^^ this is how Tiemann 2010 does it though so blame them
if strcmpi(units, 'cgs')
    aMat = aMat./aBohr;
end
end