function Q=numerov(Rgrid,E,mu,wmat_fun,wmat_vars)
%tkarman 29 feb 2012
%fixed step numerov propagator.
%
%npoints, delta the stepsize, E scattering energy. mu the reduced mass and wmat, the W-matrix at zero scattering energy.
%inputs all wmats (costs some memory, but i might be able to save them more cheaply as a cell array of sparse matrices)
%wmat a cell array of matrices, as function of R.
%but does not keep track of all the F and Q mats.

npoints=length(Rgrid);
delta=Rgrid(2)-Rgrid(1);

wmat=feval(wmat_fun,Rgrid(1),wmat_vars);
siz=size(wmat);
Q=zeros(siz);
Fmineen=eye(siz)-delta^2/12*(wmat-2*mu*E*eye(siz));
wmat=feval(wmat_fun,Rgrid(2),wmat_vars);
F=eye(siz)-delta^2/12*(wmat-2*mu*E*eye(siz));

for(i=3:npoints)
	Fmintwee=Fmineen;
	Fmineen=F;
	wmat=feval(wmat_fun,Rgrid(i),wmat_vars);
	F=eye(siz)-delta^2/12*(wmat-2*mu*E*eye(siz));
	Q=(12*eye(siz)-10*Fmineen-Fmintwee*Q)\F;
end

end
