function [V]=magndip(qnums)
%function [V]=magndip(qnums)
%tijs karman 13 feb 2012
%get matrix elements of the magnetic dipolar interaction in the scattering basis (coupled)
%this excludes a factor 1/R^3

%constants
gs=2.0023193043718;
mub=1/2;
alpha=1/(137.035999074);%CODATA 2010

%init.
[n,m]=size(qnums);
V=zeros(n,n);

%loop over the basis
for(i1=1:n)
for(i2=1:n)

%qnumbers associated with matrix element
J=qnums(i1,1);
Jp=qnums(i2,1);
j=qnums(i1,3);
jp=qnums(i2,3);
l=qnums(i1,4);
lp=qnums(i2,4);
ja=qnums(i1,5);
jap=qnums(i2,5);
jb=qnums(i1,6);
jbp=qnums(i2,6);

%actual matrix element.
if(J==Jp)
fac1=(-1)^(3/2+jap)*sqrt(30)*ff_6j([ja jap 1; 2 2 1/2])+gs*(-1)^(3/2+ja)*sqrt(3/2)*ff_6j([ja jap 1; 1/2 1/2 2]);
fac2=(-1)^(3/2+jbp)*sqrt(30)*ff_6j([jb jbp 1; 2 2 1/2])+gs*(-1)^(3/2+jb)*sqrt(3/2)*ff_6j([jb jbp 1; 1/2 1/2 2]);

V(i1,i2)=(-1)^(J+jp)*sqrt((2*ja+1)*(2*jap+1)*(2*jb+1)*(2*jbp+1)*(2*j+1)*(2*jp+1)*(2*l+1)*(2*lp+1))*ff_3jm([l 2 lp; 0 0 0])*ff_6j([j jp 2; lp l J])*ff_9j([ja jap 1; jb jbp 1; j jp 2])*fac1*fac2;
end

end
end

V=-V*(mub*alpha)^2*sqrt(30);

end
