#!/bin/csh -f

time scibat5 << EOF

// Secrest and Johnson test problem for numerov

exec('link.sci');
format(16,'v')

function [Vy, alpha, lambda]=Vy()
  delta = 0.3
  y     = (-6:delta:6)';
  V     = 0.5*y.^2;
  muy   = 1
  [lambda, U] = sinc_dvr(muy, delta, V);
  nchan = 6
  U     = U(:,1:nchan);
  lambda= lambda(1:nchan);
  
  A     = 41000;
  alpha = 0.3;

  vvy   = A*exp(alpha*y);
  Vy    = U'*diag(vvy)*U;
endfunction

function W=w_mat(R, pars)
	W=2*pars.mu*(exp(-pars.alpha*R)*pars.Vymat+diag(pars.lambda));
endfunction

mu    = 2/3

Ra    = 20
Rb    = 100
npoint= 200;
grid=Ra:(Rb-Ra)/(npoint-1):Rb;
dr=grid(2)-grid(1);
npoint=length(grid);

disp(npoint)
disp(dr)

[Vymat, alpha, lambda] = Vy();
Wpars.mu=mu;
Wpars.alpha=alpha;
Wpars.Vymat=Vymat;
Wpars.lambda=lambda;

E     = 3
tic()
Qnnum = numerov(grid, E, mu, w_mat,Wpars)
disp(toc())
tic()
Qnsin = list();
Qnsin(1) = zeros(Qnnum);
Qnsin = Qsin(grid,E,mu,Qnsin,w_mat,Wpars)
disp(toc())
tic()
Qnair = list();
Qnair(1) = zeros(Qnnum);
Qnair = Qairy(grid,E,mu,Qnair,w_mat,Wpars)
disp(toc())

disp(Qnnum)
disp(Qnsin(1))
disp(Qnair(1))

// match numerov
nchan = size(Qnnum, 1)
Ls    = zeros(nchan, 1)
Kopen = get_psiK(E, mu, lambda, Ls, Qnnum, grid(npoint), grid(npoint-1))
I     = eye(Kopen);
S     = (I + %i*Kopen)\(I-%i*Kopen)
Pnum  = abs(S).^2
// match sin
nchan = size(Qnsin(1), 1)
Ls    = zeros(nchan, 1)
Kopen = get_psiK(E, mu, lambda, Ls, Qnsin(1), grid(npoint), grid(npoint-1))
I     = eye(Kopen);
S     = (I + %i*Kopen)\(I-%i*Kopen)
Psin  = abs(S).^2
// match airy
nchan = size(Qnair(1), 1)
Ls    = zeros(nchan, 1)
Kopen = get_psiK(E, mu, lambda, Ls, Qnair(1), grid(npoint), grid(npoint-1))
I     = eye(Kopen);
S     = (I + %i*Kopen)\(I-%i*Kopen)
Pair  = abs(S).^2

disp(Pnum)
disp(Psin)
disp(Pair)

Paccurate = [0.977885643570930   0.022109316961861   0.000005039475275
   0.022109316945734   0.976992651846316   0.000898031200049
   0.000005039475270   0.000898031199724   0.999096929324841];
accu_numerov=norm(Pnum-Paccurate);
accu_sin=norm(Psin-Paccurate);
accu_air=norm(Pair-Paccurate);

disp(accu_numerov)
disp(accu_sin)
disp(accu_air)

quit
EOF
