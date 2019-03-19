function U=transformAdiabats2SFExpansion(ktabs,ctabs)
%function U=transformAdiabats2SFExpansion(ktabs,ctabs)
%
%determines the transformation between the adiabats,
%defined by c and lambda, and the spacefixed expansion
%defined by k_1, k_2, k
%
%ktabs contains S, k_1, k_2 and k per row.
%ctabs contains S, c, Lambda per row.

[nk,~]=size(ktabs);
[nc,~]=size(ctabs);
U=zeros(nk,nc);

for(ic=1:nc)
S=ctabs(ic,1);
c=ctabs(ic,2);
lambda=ctabs(ic,3);

Uqq=qqcoupled(lambda);

for(ik=1:nk)
if(S==ktabs(ik,1))
k1=ktabs(ik,2);
k2=ktabs(ik,3);
k=ktabs(ik,4);

for(L=abs(lambda):4)
for(Lp=abs(lambda):4)
U(ik,ic)=U(ik,ic)+conj(Uqq(L-abs(lambda)+1,c))*Uqq(Lp-abs(lambda)+1,c)*(-1)^(L-lambda)*ff_3jm([L k Lp; -lambda 0 lambda])*ff_9j([2 2 k1; 2 2 k2; L Lp k])*sqrt((2*k1+1)*(2*k2+1)*(2*k+1)*(2*L+1)*(2*Lp+1));
end
end

end

end
end

end
