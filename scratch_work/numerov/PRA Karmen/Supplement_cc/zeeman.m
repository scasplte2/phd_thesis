function [V]=zeeman(qnums)
%function [V]=zeeman(qnums)
%
% calculate matrix representation of the zeeman interaction (excludes scaling with the field-strength)

[n,~]=size(qnums);
V=zeros(n,n);

gs=2.0023193043718;

for(i1=1:n)
for(i2=1:n)

J=qnums(i1,1);
M=qnums(i1,2);
j=qnums(i1,3);
l=qnums(i1,4);
ja=qnums(i1,5);
jb=qnums(i1,6);
Jp=qnums(i2,1);
Mp=qnums(i2,2);
jp=qnums(i2,3);
lp=qnums(i2,4);
jap=qnums(i2,5);
jbp=qnums(i2,6);

if(l==lp)
if(jb==jbp)
	V(i1,i2)=V(i1,i2)+(-1)^(ja+jb+jp+1/2+jap)*sqrt((2*ja+1)*(2*jap+1)*30)*ff_6j([j jp 1; jap ja jb])*ff_6j([ja jap 1; 2 2 1/2])+gs*(-1)^(ja+jb+jp+1/2+ja)*sqrt((2*ja+1)*(2*jap+1)*3/2)*ff_6j([j jp 1; jap ja jb])*ff_6j([ja jap 1; 1/2 1/2 2]);
end
if(ja==jap)
	V(i1,i2)=V(i1,i2)+(-1)^(ja+jbp+j+1/2+jbp)*sqrt((2*jb+1)*(2*jbp+1)*30)*ff_6j([j jp 1; jbp jb ja])*ff_6j([jb jbp 1; 2 2 1/2])+gs*(-1)^(ja+jbp+j+1/2+jb)*sqrt((2*jb+1)*(2*jbp+1)*3/2)*ff_6j([j jp 1; jbp jb ja])*ff_6j([jb jbp 1; 1/2 1/2 2]);
end
V(i1,i2)=V(i1,i2)/2*(-1)^(J+Jp-M+j+l+1)*sqrt((2*J+1)*(2*Jp+1)*(2*j+1)*(2*jp+1))*ff_3jm([J 1 Jp; -M 0 Mp])*ff_6j([J Jp 1; jp j l]);
end

end
end


end
