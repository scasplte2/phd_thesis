#!/bin/csh -f

time scibat5 << EOF

exec('link.sci');

function W=w_mat(R, pars)
	W=pars.W0+pars.W1*R;
endfunction

mu    = 1/2;

Ra    = 0;
Rb    = 2*%pi;
npoint= 5;
grid=Ra:(Rb-Ra)/(npoint-1):Rb;
dr=grid(2)-grid(1);
npoint=length(grid);

disp(npoint)
disp(dr)

Wpars.mu=mu;
Wpars.W0=-5;
Wpars.W1=1;

E     = 0;
Qnnum = numerov(grid, E, mu, w_mat,Wpars)
Qnsin = list();
Qnsin(1) = zeros(Qnnum);
Qnsin = Qsin(grid,E,mu,Qnsin,w_mat,Wpars)
Qnair = list();
Qnair(1) = zeros(Qnnum);
Qnair = Qairy(grid,E,mu,Qnair,w_mat,Wpars)

disp(Qnnum)
disp(Qnsin(1))
disp(Qnair(1))

b=-airyai(nthroot(Wpars.W1,3)*(Wpars.W0/Wpars.W1+grid(1)))/airybi(nthroot(Wpars.W1,3)*(Wpars.W0/Wpars.W1+grid(1)));
Qexact=(airyai(nthroot(Wpars.W1,3)*(Wpars.W0/Wpars.W1+grid(npoint-1)))+b*airybi(nthroot(Wpars.W1,3)*(Wpars.W0/Wpars.W1+grid(npoint-1))))/(airyai(nthroot(Wpars.W1,3)*(Wpars.W0/Wpars.W1+grid(npoint)))+b*airybi(nthroot(Wpars.W1,3)*(Wpars.W0/Wpars.W1+grid(npoint))))
disp(Qexact)

quit
EOF
