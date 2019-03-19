#!/bin/csh -f

time scibat5 << EOF

exec('link.sci');

function W=w_mat(R, pars)
	W=pars.kk;
endfunction

mu    = 1/2;

Ra    = 0;
Rb    = 6;
npoint= 5;
grid=Ra:(Rb-Ra)/(npoint-1):Rb;
dr=grid(2)-grid(1);
npoint=length(grid);

disp(npoint)
disp(dr)

Wpars.mu=mu;
Wpars.kk=-5;

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

if ( Wpars.kk < 0 )
	Qexact=sin(sqrt(abs(Wpars.kk))*grid(npoint-1))/sin(sqrt(abs(Wpars.kk))*grid(npoint));
else
        Qexact=sinh(sqrt(abs(Wpars.kk))*grid(npoint-1))/sinh(sqrt(abs(Wpars.kk))*grid(npoint));
end
disp(Qexact)

quit
EOF
