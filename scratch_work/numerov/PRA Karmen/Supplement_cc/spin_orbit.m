function [V] = spin_orbit(qnums)
%tkarman 29 feb 2012
%spin orbit in scattering basis (coupled)

%constants.
A=168.34/219474.63137098*(2/5);% experimental spin-orbit splitting for scandium in a.u.

%init.
[n,m]=size(qnums);
V=zeros(n,n);

%loop over matrix elements (all non-zero ones are on the diagonal)
for(i=1:n)
ja=qnums(i,5);
jb=qnums(i,6);

V(i,i)=A/2*(ja*(ja+1)+jb*(jb+1)-27/2);
end

end

