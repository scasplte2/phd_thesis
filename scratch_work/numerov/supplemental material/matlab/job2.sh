#!/bin/csh -f

time matlab -nojvm -nosplash << EOF

% Linear potential test. Airy should be exact, Qsin and Numerov should not. (except for if you take W1=0, in which case this should correspond to job1.sh)
format long
mu    = 0.5;

npoint=5;
Rgrid=0:2*pi/(npoint-1):2*pi;
npoint=length(Rgrid);
fprintf('Grid ranges from %f to %f and contains %d points.\n',Rgrid(1),Rgrid(npoint),npoint);

pars.W0=-5;
pars.W1=1;

tic();
Qnnum = numerov(Rgrid, 0, mu, @wmat_lin,pars);
toc()
[nchan,~] = size(Qnnum);
Ls    = zeros(nchan, 1);
tic();
Qnsin{1}=zeros(nchan,nchan);
Qnsin = Qsin(Rgrid,0,mu,Qnsin,@wmat_lin,pars);
toc()
tic();
Qnair{1}=zeros(nchan,nchan);
Qnair = Qairy(Rgrid,0,mu,Qnair,@wmat_lin,pars);
toc()

%match numerov
[S1,Kopen] = matchBoundaryConditions(Rgrid,Ls,mu,Qnnum,0,0);
I     = eye(size(Kopen));
Snum     = (I + sqrt(-1)*Kopen)\(I-sqrt(-1)*Kopen);

%match sine
[S1,Kopen] = matchBoundaryConditions(Rgrid,Ls,mu,Qnsin{1},0,0);
I     = eye(size(Kopen));
Ssin     = (I + sqrt(-1)*Kopen)\(I-sqrt(-1)*Kopen);

%match airy
[S1,Kopen] = matchBoundaryConditions(Rgrid,Ls,mu,Qnair{1},0,0);
I     = eye(size(Kopen));
Sair     = (I + sqrt(-1)*Kopen)\(I-sqrt(-1)*Kopen);

Qnnum
Qnsin=Qnsin{1}
Qnair=Qnair{1}
b=-airy(0,nthroot(pars.W1,3)*(pars.W0)/pars.W1)/airy(2,nthroot(pars.W1,3)*(pars.W0)/pars.W1);
Qexact=(airy(0,nthroot(pars.W1,3)*((pars.W0)/pars.W1+Rgrid(npoint-1)-Rgrid(1)))+b*airy(2,nthroot(pars.W1,3)*((pars.W0)/pars.W1+Rgrid(npoint-1)-Rgrid(1))))/(airy(0,nthroot(pars.W1,3)*((pars.W0)/pars.W1+Rgrid(npoint)-Rgrid(1)))+b*airy(2,nthroot(pars.W1,3)*((pars.W0)/pars.W1+Rgrid(npoint)-Rgrid(1))))

quit
EOF
