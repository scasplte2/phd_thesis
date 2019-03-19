#!/bin/csh -f

time matlab -nojvm -nosplash << EOF

% Secrest and Johnson test problem for numerov
format long
mu    = 2/3;

npoint=200;
Rgrid=logspace(log10(20),log10(100),npoint);
Rgrid=20:80/(npoint-1):100;
npoint=length(Rgrid);
fprintf('Grid ranges from %f to %f and contains %d points.\n',Rgrid(1),Rgrid(npoint),npoint);

[Vymat, alpha, lambda] = Vy();
pars.mu=mu;
pars.alpha=alpha;
pars.Vymat=Vymat;
pars.lambda=lambda;

E     = 3;
tic();
Qnnum = numerov(Rgrid, E, mu, @wmat_sj,pars);
toc()
[nchan,~] = size(Qnnum);
Ls    = zeros(nchan, 1);
tic();
Qnsin{1}=zeros(nchan,nchan);
Qnsin = Qsin(Rgrid,E,mu,Qnsin,@wmat_sj,pars);
toc()
tic();
Qnair{1}=zeros(nchan,nchan);
Qnair = Qairy(Rgrid,E,mu,Qnair,@wmat_sj,pars);
toc()

Qnnum
Qnsin{1}
Qnair{1}

%match numerov
[S1,Kopen] = matchBoundaryConditions(Rgrid,Ls,mu,Qnnum,E,lambda);
I     = eye(size(Kopen));
S2     = (I + sqrt(-1)*Kopen)\(I-sqrt(-1)*Kopen);
Pnum     = abs(S1).^2

%match sine
[S1,Kopen] = matchBoundaryConditions(Rgrid,Ls,mu,Qnsin{1},E,lambda);
I     = eye(size(Kopen));
S2     = (I + sqrt(-1)*Kopen)\(I-sqrt(-1)*Kopen);
Psin     = abs(S1).^2

%match airy
[S1,Kopen] = matchBoundaryConditions(Rgrid,Ls,mu,Qnair{1},E,lambda);
I     = eye(size(Kopen));
S2     = (I + sqrt(-1)*Kopen)\(I-sqrt(-1)*Kopen);
Pair     = abs(S1).^2


%accurate value from 1000 points Numerov
Paccurate = [0.977885643570930   0.022109316961861   0.000005039475275
   0.022109316945734   0.976992651846316   0.000898031200049
   0.000005039475270   0.000898031199724   0.999096929324841];
accu_numerov=norm(Pnum-Paccurate)
accu_sin=norm(Psin-Paccurate)
accu_air=norm(Pair-Paccurate)

quit
EOF
