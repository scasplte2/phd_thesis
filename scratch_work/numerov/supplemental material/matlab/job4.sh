#!/bin/csh -f

time matlab -nojvm -nosplash << EOF

% Secrest and Johnson test problem for numerov
format long
mu    = 2/3;

npoint=1000;
Rgrid=logspace(log10(20),log10(100),npoint);
Rgrid=20:80/(npoint-1):100;
npoint=length(Rgrid);
fprintf('Grid ranges from %f to %f and contains %d points.\n',Rgrid(1),Rgrid(npoint),npoint);

[Vymat, alpha, lambda] = Vy();
pars.mu=mu;
pars.alpha=alpha;
pars.Vymat=Vymat;
pars.lambda=lambda;

E     = [2 3 5];
tic();
for(iE=1:length(E))
	Qntmp = numerov(Rgrid, E(iE), mu, @wmat_sj,pars);
	Qnnum{iE}=Qntmp;

	[nchan,~] = size(Qnnum{iE});
	Ls    = zeros(nchan, 1);
	[S1,Kopen] = matchBoundaryConditions(Rgrid,Ls,mu,Qnnum{iE},E(iE),lambda);
	I     = eye(size(Kopen));
	S2{iE}     = (I + sqrt(-1)*Kopen)\(I-sqrt(-1)*Kopen);
	Pnum{iE}   = abs(S1).^2;
end
toc()

RgridSHORT=20:80/(npoint-1):60;
npoint_short=length(RgridSHORT);
RgridLONG=[Rgrid(npoint_short-1) logspace(log10(Rgrid(npoint_short)),log10(100),npoint_short/2)];
npoint_long=length(RgridLONG);
fprintf('Variable grid ranges from %f to %f and contains %d points.\n',RgridSHORT(1),RgridLONG(npoint_long),npoint_short+npoint_long);

for(iE=1:length(E))
        Qntmp = numerov(RgridSHORT, E(iE), mu, @wmat_sj,pars);
        Qnnum{iE}=Qntmp;
end
tic();
Qnsin = Qsin(RgridLONG,E,mu,Qnnum,@wmat_sj,pars);
for(iE=1:length(E))
        [nchan,~] = size(Qnsin{iE});
        Ls    = zeros(nchan, 1);
        [S1,Kopen] = matchBoundaryConditions(RgridLONG,Ls,mu,Qnsin{iE},E(iE),lambda);
        I     = eye(size(Kopen));
        S2{iE}     = (I + sqrt(-1)*Kopen)\(I-sqrt(-1)*Kopen);
        Psin{iE}   = abs(S1).^2;
end
toc();
tic();
Qnair = Qairy(RgridLONG,E,mu,Qnnum,@wmat_sj,pars);
for(iE=1:length(E))
        [nchan,~] = size(Qnair{iE});
        Ls    = zeros(nchan, 1);
        [S1,Kopen] = matchBoundaryConditions(RgridLONG,Ls,mu,Qnair{iE},E(iE),lambda);
        I     = eye(size(Kopen));
        S2{iE}     = (I + sqrt(-1)*Kopen)\(I-sqrt(-1)*Kopen);
        Pair{iE}   = abs(S1).^2;
end
toc();

%accurate value from 1000 points Numerov
Paccurate{1} = [0.999547223917325   0.000452776082757
   0.000452776082593   0.999547223917324];
Paccurate{2} = [0.977885643570930   0.022109316961861   0.000005039475275
   0.022109316945734   0.976992651846316   0.000898031200049
   0.000005039475270   0.000898031199724   0.999096929324841];
Paccurate{3} = [0.732562099675081   0.252051177754239   0.015269998373490 0.000116708189670   0.000000016210761
   0.252051177387060   0.562235710636136   0.182403644363880 0.003308724237944   0.000000743293936
   0.015269998334567   0.182403644164815   0.743026369721991 0.059271964293274   0.000028023387949
   0.000116708189292   0.003308724231918   0.059271964250125 0.935538908548919   0.001763694755294
   0.000000016210761   0.000000743293934   0.000028023387918 0.001763694754646   0.998207522352400];

for(iE=1:length(E))
	accu_numerov(iE)=norm(Pnum{iE}-Paccurate{iE});
        accu_sin(iE)=norm(Psin{iE}-Paccurate{iE});
        accu_airy(iE)=norm(Pair{iE}-Paccurate{iE});
end

accu_numerov
accu_sin
accu_airy

quit
EOF
