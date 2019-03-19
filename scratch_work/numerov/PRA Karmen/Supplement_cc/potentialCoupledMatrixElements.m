function [V,ktabs]=potentialCoupledMatrixElements(qnums)
%29 feb 2012
%matrix elements of the interatomic interaction potential in the coupled scattering basis.
%jpc 108 41 8945
%returns V with V(:,:,i) the potential except for the R depedent expansion coefficient.
%i labels S, k1, k2 and k, defined in ktabs.

%init
[n,~]=size(qnums);
load eightteenj.mat

%maak ktabs
ktabs=zeros(0);
for(S=0:1)
for(k1=0:4)
for(k2=0:4)
for(k=abs(k1-k2):(k1+k2))

	if(bitget(k,1)==0)
		ktabs=[ktabs; S k1 k2 k];
	end

end
end
end
end
[m,~]=size(ktabs);

%init V
V=zeros(n,n,m);

%loop ktabs
for(ik=1:m)
S=ktabs(ik,1);
k1=ktabs(ik,2);
k2=ktabs(ik,3);
k=ktabs(ik,4);

%loop over matrix
for(i1=1:n)
for(i2=1:n)

%read qnums
if(qnums(i1,1)==qnums(i2,1) && qnums(i1,2)==qnums(i2,2)) %delta_JJ' delta_MM'
	J=qnums(i1,1);
	j=qnums(i1,3);
	l=qnums(i1,4);
	ja=qnums(i1,5);
	jb=qnums(i1,6);
	jp=qnums(i2,3);
	lp=qnums(i2,4);
	jap=qnums(i2,5);
	jbp=qnums(i2,6);

	if(k>=abs(j-jp) && k>=abs(l-lp) && k <= j+jp && k <= l+lp)

               V(i1,i2,ik)=(-1)^(ja+jap+3-jp+k1-k2-k-J)*sqrt((2*ja+1)*(2*jap+1)*(2*jb+1)*(2*jbp+1)*(2*j+1)*(2*jp+1)*(2*l+1)*(2*lp+1)*(2*k1+1)*(2*k2+1)*(2*k+1))*(2*S+1)*ff_3jm([l k lp; 0 0 0])*ff_6j([j jp k; lp l J])*eightteenj(S+1,k1+1,k2+1,k+1,ja-1/2,jb-1/2,jap-1/2,jbp-1/2,j+1,jp+1);

	end
end

end
end

end

end
