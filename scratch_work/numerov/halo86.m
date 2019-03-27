ang2Bohr = 1.889726877774355;
cm2H   = 4.556335270832134e-06;
amu2Au = 1.822889e+03;

m86 = 85.9092607309*amu2Au;
mu86 = prod([m86 m86])/sum([m86 m86]);

C6  = 1.52701e7 * (cm2H*ang2Bohr^6);  % Ha * (aBohr)^6

Rvdw = @(mu) (1/2)*(2*mu*C6)^(1/4);
aBar = @(mu, Rvdw) 4*pi*Rvdw/gamma(1/4)^2;

g1 = gamma(1/4)^4/(6*pi^2)-2;
g2 = (5/4)*g1^2-2; 
% From Chin 2010 eq. 33
funcEhalo = @(mu, a, aBar) 1./(2*mu*(a-aBar).^2).*(1 + g1*aBar/(a-aBar) + g2*aBar^2/(a-aBar).^2);

Rvdw86 = Rvdw(mu86);
aBar86 = aBar(mu86, Rvdw(mu86));

srScat = calcScatterLengths;
a86 = srScat.sr86(2);

eBind86_kHz = funcEhalo(mu86, a86, aBar86)*6.57968e12;