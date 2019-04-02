function param = toyLJ
% Complimentary function to the jupyter notebook

%% Physical constants
ev2Ha = 3.67493e-2;
cm2Ha = 4.55634e-6;
cm2MHz = 2.99792e4;
ev2cm = 8065.541;
ev2MHz = 2.41799e8;
aBohr = 0.52917721067;
kg2au = 1/9.10939e-31;
Ha2MHz = 6.57968e9;
amu2Au   = 1.822889e+03;

%% Problem specifics
rm = 1/aBohr;
%De = 1005*cm2Ha;

m86 = 85.9092607309*amu2Au;
mu = prod([m86 m86])/sum([m86 m86]);

%% Function definitions
% Lennard-Jones potential
funcV = @(R, De) De*( (rm./R).^12 - 2*(rm./R).^6 );

% mean scattering length
aBar = @(Rvdw) 4*pi*Rvdw/gamma(1/4)^2;

% s-wave scattering length
a    = @(aBar, phiD) aBar*(1-tan(pi/4)*tan(phiD - pi/8));

% van der Waals length
Rvdw = @(mu, C6) (1/2)*(2*mu*C6)^(1/4);

% van der Waals Energy
Evdw = @(mu, Rvdw) 1/(2*mu*Rvdw^2);

% classical turning point
rClas = @(a, Rvdw) (a*(2*Rvdw)^2)^(1/3);

% potential phase contribution
phiD = @(mu, r0, De) integral(@(x) sqrt(-2*mu*funcV(x, De)), r0, Inf, 'ArrayValued', 1, 'AbsTol', 1e-10, 'RelTol', 1e-10);

% From Chin 2010 eq. 33
g1 = gamma(1/4)^4/(6*pi^2)-2;
g2 = (5/4)*g1^2-2; 
funcEhalo = @(mu, a, aBar) 1./(2*mu*(a-aBar).^2).*(1 + g1*aBar/(a-aBar) + g2*aBar^2/(a-aBar).^2);

%% Main loop
%in = linspace(750, 1150, 2e3);
%in = linspace(900, 975, 1e3);
in = 948;
[ param.a, param.aBar, param.eBind_MHz, param.rClas,...
  param.Rvdw, param.Evdw, param.phiD ] = deal(zeros(1,length(in)));
for i = 1:length(in)

    De = in(i)*cm2Ha;
    C6 = 2*De*rm^6;

    % find PEC zero crossing (assume near R_i)
    minOpt = optimset('TolX', 1e-7); % set Tol to ensure that the integral below is negative
    r0     = fminbnd(@(x) abs(funcV(x, De)), rm*0.5, rm, minOpt);


    param.Rvdw(i)      = Rvdw(mu, C6);
    param.Evdw(i)      = Evdw(mu, param.Rvdw(i));
    param.aBar(i)      = aBar(param.Rvdw(i));
    param.phiD(i)      = phiD( mu, r0, De );
    
    param.a(i)         = a(param.aBar(i), param.phiD(i));
    param.eBind_MHz(i) = funcEhalo(mu, param.a(i), param.aBar(i))*Ha2MHz;
    param.rClas(i)     = rClas( param.a(i), param.Rvdw(i) );
    
end


figure;
subplot(2, 1, 1)
plot(in, param.a);
axis tight
ylim([-500 500])
ylabel('Scattering Length [a_0]')

subplot(2, 1, 2)
plot(in, param.eBind_MHz);
axis tight
ylim([0 max(param.eBind_MHz)])
ylabel('Halo Binding energy [MHz]')
end