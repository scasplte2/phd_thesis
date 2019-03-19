function [Udiag] = numDiagHas(qnumsunc,Uatom,permsym)
%3 april tkarman
%
%function [Udiag] = numDiagHas(qnumsunc,Uatom,permsym)
% returns the transformation U that transforms ja->na jb->nb
% the new quantum numbers have the same values as the old ones
% stored in qnumsunc. In this sense, na=3/2 means the eigenstate
% which corresponds almost exclusively to ja=3/2.
%
% Uatom is a struct which contains the definitions of na in terms of
% ja, at the current field strength, at a certain value of mja
% which labels the struct Uatom{mja+7/2}.

[n,~]=size(qnumsunc);

Udiag=zeros(n);
for(i1=1:n)
for(i2=1:n)

if( qnumsunc(i1,5)==qnumsunc(i2,5) && qnumsunc(i1,6)==qnumsunc(i2,6) )
	%read qnums.
	ja=qnumsunc(i1,1);
	mja=qnumsunc(i1,2);
	jb=qnumsunc(i1,3);
	mjb=qnumsunc(i1,4);
	na=qnumsunc(i2,1);
	mna=qnumsunc(i2,2);
	nb=qnumsunc(i2,3);
	mnb=qnumsunc(i2,4);
	l=qnumsunc(i1,5);
	%normalization
	if(ja==jb && mja==mjb)
		if(na==nb&&mja==mjb)
			fac=1/2;
		else
			fac=sqrt(2)/2;
		end
	else
		if(na==nb&&mja==mjb)
                        fac=sqrt(2)/2;
                else
                        fac=1;
                end
	end
	
	%matrix element.
	if(mja==mna && mjb==mnb)
		Udiag(i1,i2)=Uatom{mja+7/2}(ja-1/2,na-1/2)*Uatom{mjb+7/2}(jb-1/2,nb-1/2);
	end
	if(mja==mnb && mjb==mna)
		Udiag(i1,i2)=Udiag(i1,i2)+permsym*(-1)^l*Uatom{mja+7/2}(ja-1/2,nb-1/2)*Uatom{mjb+7/2}(jb-1/2,na-1/2);
	end
	Udiag(i1,i2)=Udiag(i1,i2)*fac;
end

end
end

end
