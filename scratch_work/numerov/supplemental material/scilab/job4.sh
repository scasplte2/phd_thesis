#!/bin/csh -f

time scibat5 << EOF

// Secrest and Johnson test problem for numerov

exec('link.sci');
format(16,'v');

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
Rb    = 60
npoint1= 500;
grid1=Ra:(Rb-Ra)/(npoint1-1):Rb;
dr=grid1(2)-grid1(1);
npoint1=length(grid1);
grid2=[grid1(npoint1-1) logspace(log10(60),log10(100),100)];
npoint2=length(grid2);

disp(npoint1)
disp(npoint2)
disp(dr)

[Vymat, alpha, lambda] = Vy();
Wpars.mu=mu;
Wpars.alpha=alpha;
Wpars.Vymat=Vymat;
Wpars.lambda=lambda;

E     = [2 3 5];
tic()
Qnnum=list();
Qtmp = numerov(grid1, E(1), mu, w_mat,Wpars)
Qnnum(1)=Qtmp
Qtmp = numerov(grid1, E(2), mu, w_mat,Wpars)
Qnnum(2)=Qtmp
Qtmp = numerov(grid1, E(3), mu, w_mat,Wpars)
Qnnum(3)=Qtmp
disp(toc())
tic()
Qnsin=Qnnum;
Qnsin = Qsin(grid2,E,mu,Qnsin,w_mat,Wpars)
disp(toc())
tic()
Qnair=Qnnum;
Qnair = Qairy(grid2,E,mu,Qnair,w_mat,Wpars)
disp(toc())

// match sin
Psin=list();
for(iE=1:length(E))
	nchan    = size(Qnsin(iE), 1)
	Ls       = zeros(nchan, 1)
	Kopen    = get_psiK(E(iE), mu, lambda, Ls, Qnsin(iE), grid2(npoint2), grid2(npoint2-1))
	I        = eye(Kopen);
	S        = (I + %i*Kopen)\(I-%i*Kopen)
	Psin(iE) = abs(S).^2
end
// match airy
Pair=list();
for(iE=1:length(E))
	nchan    = size(Qnair(iE), 1)
	Ls       = zeros(nchan, 1)
	Kopen    = get_psiK(E(iE), mu, lambda, Ls, Qnair(iE), grid2(npoint2), grid2(npoint2-1))
	I        = eye(Kopen);
	S        = (I + %i*Kopen)\(I-%i*Kopen)
	Pair(iE) = abs(S).^2
end

Paccurate=list();
Paccurate(1) = [0.999547223917325   0.000452776082757
   0.000452776082593   0.999547223917324];
Paccurate(2) = [0.977885643570930   0.022109316961861   0.000005039475275
   0.022109316945734   0.976992651846316   0.000898031200049
   0.000005039475270   0.000898031199724   0.999096929324841];
Paccurate(3) = [0.732562099675081   0.252051177754239   0.015269998373490 0.000116708189670   0.000000016210761
   0.252051177387060   0.562235710636136   0.182403644363880 0.003308724237944   0.000000743293936
   0.015269998334567   0.182403644164815   0.743026369721991 0.059271964293274   0.000028023387949
   0.000116708189292   0.003308724231918   0.059271964250125 0.935538908548919   0.001763694755294
   0.000000016210761   0.000000743293934   0.000028023387918 0.001763694754646   0.998207522352400];

for(iE=1:length(E))
        accu_sin(iE)=norm(Psin(iE)-Paccurate(iE));
        accu_airy(iE)=norm(Pair(iE)-Paccurate(iE));
end

disp(accu_sin)
disp(accu_airy)

quit
EOF
