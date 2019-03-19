function [Q]=numerov6(Rgrid,E,mu,wmat_fun,wmat_vars)
%tkarman 23 juli 2012
%fixed step numerov propagator.

delta=Rgrid(2)-Rgrid(1);
npoints=length(Rgrid);
wmat1=feval(wmat_fun,Rgrid(1),wmat_vars);
siz=size(wmat1);
Q=zeros(siz);
Fmineen=eye(siz)-delta^2/12*(wmat1-2*mu*E*eye(siz));
F=eye(siz)-delta^2/12*(feval(wmat_fun,Rgrid(2),wmat_vars)-2*mu*E*eye(siz));
clear wmat1;

for(i=3:npoints)
Fmintwee=Fmineen;
Fmineen=F;
F=eye(siz)-delta^2/12*(feval(wmat_fun,Rgrid(i),wmat_vars)-2*mu*E*eye(siz));
Q=(12*eye(siz)-10*Fmineen-Fmintwee*Q)\F;
end

end
