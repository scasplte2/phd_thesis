function [S,Kopen] = matchBoundaryConditions(Rgrid,ls,mu,Qn,E,Echannel)
%function [S] = matchBoundaryConditions(Rgrid,ls,mu,Qn,E,Echannel)
% improved scaling of columns and bessel functions.

%read channels
k=(2*mu*(E-Echannel)).^(1/2);
closed=find(E<Echannel);
open=find(E>=Echannel);
nclosed=length(closed);
nopen=length(open);
npoints=length(Rgrid);
nstates=length(Echannel);

%init
Fn=zeros(nstates,nopen);
Gn=zeros(nstates,nstates);
Fn1=Fn;
Gn1=Gn;

%build bessels. 
%improve by some scaling option exp(kr) or something. besselx(l,k*R,1)
for(i1=1:nopen)
       Fn(open(i1),open(i1))=sqrt(mu*pi*Rgrid(npoints)/2)*besselj(ls(open(i1))+1/2,k(open(i1))*Rgrid(npoints));
       Fn1(open(i1),open(i1))=sqrt(mu*pi*Rgrid(npoints-1)/2)*besselj(ls(open(i1))+1/2,k(open(i1))*Rgrid(npoints-1));
       Gn(open(i1),open(i1))=sqrt(mu*pi*Rgrid(npoints)/2)*bessely(ls(open(i1))+1/2,k(open(i1))*Rgrid(npoints));
       Gn1(open(i1),open(i1))=sqrt(mu*pi*Rgrid(npoints-1)/2)*bessely(ls(open(i1))+1/2,k(open(i1))*Rgrid(npoints-1));
end
for(i1=1:nclosed)
       Gn(closed(i1),closed(i1))=sqrt(mu*pi*Rgrid(npoints)/2)*besselk(ls(closed(i1))+1/2,k(closed(i1))*Rgrid(npoints),1);
       Gn1(closed(i1),closed(i1))=sqrt(mu*pi*Rgrid(npoints-1)/2)*besselk(ls(closed(i1))+1/2,k(closed(i1))*Rgrid(npoints-1),1)*exp(k(closed(i1))*(Rgrid(npoints)-Rgrid(npoints-1)));
end

% K matrix 
%K=(Gnmineen-Qn*Gn)\(Fnmineen-Qn*Fn);
B=Gn1-Qn*Gn;
nrm=zeros(nstates,1);
for(i1=1:nstates)
	nrm(i1)=norm(B(:,i1));
	B(:,i1)=B(:,i1)/nrm(i1);
end
K=B\(Fn1-Qn*Fn);
for(i1=1:nstates)
	K(i1,:)=K(i1,:)/nrm(i1);
end

% get scattering matrix.
Kopen=K(open,open);
S=(eye(size(Kopen))+sqrt(-1)*Kopen)\(eye(size(Kopen))-sqrt(-1)*Kopen);


end
