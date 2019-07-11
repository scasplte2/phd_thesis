function eBind86_kHz = halo86(model, varargin)
% Given a particular model of the potential, claculate some relevant quantities for the halo state

%% Determine model to use
% Aman '18 value from halo binding energy, Tiemann 2010 value from fourier spec.
if ~exist('model','var'); model = ''; end
switch model
    case {'Tiemann - Free', 'TF'}
        c6Model = 1.52701; 
    case {'Tiemann - Rec', 'TR'}
        c6Model = 1.525; % 3164 au
    case {'Aman - Free', 'AF'}
        c6Model = 1.526593; % 3167.61 au
    case 'varyC6'
        c6Model = varargin{2};
    otherwise %default to 'Aman - Rec'
        fprintf('\n Using defaults - Model: Recommended, C6: Aman\n')
        c6Model = 1.523967; % 3162.16 au
        model = 'AR';
end

%% Setup constants
ang2Bohr = 1/0.52917721067;
cm2H     = 4.556335270832134e-06;
amu2Au   = 1.822889e+03;

m86  = 85.9092607309*amu2Au;
mu86 = prod([m86 m86])/sum([m86 m86]);

C6 = c6Model * 1e7 * (cm2H*ang2Bohr^6)  % Ha * (aBohr)^6

Rvdw  = @(mu) (1/2)*(2*mu*C6)^(1/4);
aBar  = @(mu, Rvdw) 4*pi*Rvdw/gamma(1/4)^2;
rClas = @(a, Rvdw) (a*(2*Rvdw)^2)^(1/3);

%% Calculate quantities
g1 = gamma(1/4)^4/(6*pi^2)-2;
g2 = (5/4)*g1^2-2; 
% From Chin 2010 eq. 33
funcEhalo = @(mu, a, aBar) 1./(2*mu*(a-aBar).^2).*(1 + g1*aBar/(a-aBar) + g2*aBar^2/(a-aBar).^2);

Rvdw86 = Rvdw(mu86);
aBar86 = aBar(mu86, Rvdw(mu86));

srScat = calcScatterLengths(model, varargin{:}, 'noPrint');
a86 = srScat.sr86(2)

rClas86 = rClas(a86, Rvdw86)

eBind86_kHz = funcEhalo(mu86, a86, aBar86)*6.57968e12;
end