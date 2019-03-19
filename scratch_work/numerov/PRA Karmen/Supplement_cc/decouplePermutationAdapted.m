function [U,qnumsuncoupled]=decouplePermutationAdapted(qnumscoupled,permsym)
%tkarman 15 maart 2012
%function [U,qnumsuncoupled]=decouplePermutationAdapted(qnumscoupled,permsym)
%
% same as decouple but for a basis adapted to permutation symmetry. (ja >= jb)

[n,~]=size(qnumscoupled);
M=qnumscoupled(1,2); %interaction doesnt couple M, so there should be only one M, equal to ml+mja+mjb in the uncoupled basis.

qnumsuncoupled=zeros(0);
for(l=unique(qnumscoupled(:,4))')
	for(ja=[3/2 5/2])
		for(jb=3/2:ja)
			for(ml=-l:l)
			for(mja=-ja:ja)
			for(mjb=-jb:jb)
		
				if(mja+mjb+ml==M)

					if(ja==jb)
						if(mja==mjb)
							if((-1)^l==permsym)
								qnumsuncoupled=[qnumsuncoupled; ja mja jb mjb l ml];
							end
						elseif(mja>mjb)
							qnumsuncoupled=[qnumsuncoupled; ja mja jb mjb l ml];	
						end
					else
						qnumsuncoupled=[qnumsuncoupled; ja mja jb mjb l ml];
					end
				end

			end
			end
			end
		end
	end
end
[m,~]=size(qnumsuncoupled);

%U=zeros(n,m);
%for(i1=1:n)
%for(i2=1:m)

%	l=qnumscoupled(i1,4);
%	ja=qnumscoupled(i1,5);
%	jb=qnumscoupled(i1,6);
%	if(qnumsuncoupled(i2,1) == ja && qnumsuncoupled(i2,3) == jb && qnumsuncoupled(i2,5) == l)
 %       	J=qnumscoupled(i1,1);
%	        j=qnumscoupled(i1,3);
%	        mja=qnumsuncoupled(i2,2);
%	        mjb=qnumsuncoupled(i2,4);
%	        ml=qnumsuncoupled(i2,6);

        	
%		U(i1,i2)=c_gm([ja jb j; mja mjb mja+mjb])*c_gm([j l J; mja+mjb ml M]);
%	end

%	if(qnumsuncoupled(i2,1) == jb && qnumsuncoupled(i2,3) == ja && qnumsuncoupled(i2,5) == l)
%		J=qnumscoupled(i1,1);
%               j=qnumscoupled(i1,3);
%		mja=qnumsuncoupled(i2,2);
%               mjb=qnumsuncoupled(i2,4);
%               ml=qnumsuncoupled(i2,6);

%		U(i1,i2)=U(i1,i2)+permsym*(-1)^(l+ja+jb-j)*c_gm([ja jb j; mja mjb mja+mjb])*c_gm([j l J; mja+mjb ml M]);
%	end

%	if(ja==jb)
%		U(i1,i2)=U(i1,i2)/sqrt(2);
%	end
%	if(qnumsuncoupled(i2,1)==qnumsuncoupled(i2,3) && qnumsuncoupled(i2,2) == qnumsuncoupled(i2,4))
%		U(i1,i2)=U(i1,i2)/sqrt(2);
%	end
%end
%end
%V=U(:,sum(U~=0,1)~=0);
%Tuncoupled=V'*Tcoupled*V;

%retry:
U=zeros(n,m);
for(i1=1:n)
for(i2=1:m)
	J=qnumscoupled(i1,1);
	M=qnumscoupled(i1,2);
	jab=qnumscoupled(i1,3);
        l=qnumscoupled(i1,4);
        ja=qnumscoupled(i1,5);
        jb=qnumscoupled(i1,6);

	jap=qnumsuncoupled(i2,1);
        map=qnumsuncoupled(i2,2);
        jbp=qnumsuncoupled(i2,3);
        mbp=qnumsuncoupled(i2,4);
        lp=qnumsuncoupled(i2,5);
        mlp=qnumsuncoupled(i2,6);

	if(l==lp)
		U(i1,i2)=c_gm([ja jb jab; map mbp map+mbp])*c_gm([jab l J; map+mbp mlp M]);
		fac=0;
		if(ja==jap && jb == jbp)
			fac = 1;
		end
		if(ja==jbp && jb == jap)
			fac = fac+ permsym*(-1)^l;
		end
		if(ja==jb)
			fac=fac/sqrt(2);
		end
		if(jap==jbp && map == mbp)
			fac=fac/sqrt(2);
		end
		U(i1,i2)=U(i1,i2)*fac;
	end
end
end


end
