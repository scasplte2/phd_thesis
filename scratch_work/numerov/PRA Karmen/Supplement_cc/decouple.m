function [U,qnumsuncoupled]=decouple(qnumscoupled)
%tkarman 15 maart 2012

[n,~]=size(qnumscoupled);
M=qnumscoupled(1,2); %interaction doesnt couple M, so there should be only one M, equal to ml+mja+mjb in the uncoupled basis.

qnumsuncoupled=zeros(0);
for(l=unique(qnumscoupled(:,4))')
	for(ja=[3/2 5/2])
		for(jb=[3/2 5/2])
			for(ml=-l:l)
			for(mja=-ja:ja)
			for(mjb=-jb:jb)
		
				if(mja+mjb+ml==M)
					qnumsuncoupled=[qnumsuncoupled; ja mja jb mjb l ml];
				end

			end
			end
			end
		end
	end
end
[m,~]=size(qnumsuncoupled);

U=zeros(n,m);
for(i1=1:n)
for(i2=1:m)

	l=qnumscoupled(i1,4);
	ja=qnumscoupled(i1,5);
	jb=qnumscoupled(i1,6);
	if(qnumsuncoupled(i2,1) == ja && qnumsuncoupled(i2,3) == jb && qnumsuncoupled(i2,5) == l)
        	J=qnumscoupled(i1,1);
	        j=qnumscoupled(i1,3);
	        mja=qnumsuncoupled(i2,2);
	        mjb=qnumsuncoupled(i2,4);
	        ml=qnumsuncoupled(i2,6);

        	
		U(i1,i2)=c_gm([ja jb j; mja mjb mja+mjb])*c_gm([j l J; mja+mjb ml M]);
	end

end
end

end
