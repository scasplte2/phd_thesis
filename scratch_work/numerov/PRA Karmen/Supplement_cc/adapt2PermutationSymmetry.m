function [U]=adapt2PermutationSymmetry(qnums,qnumssym,permsym)
%function [U]=adapt2PermutationSymmetry(qnums,qnumssym,permsym)
%
% generates the transformation between a primitive basis (defined in qnums)
% and a permuation symmetry adapted basis (qnumssym).
% the total (not nuclear, difference is a factor (-1)^n_e=-1 for Sc) permuation symmetry is provided in permsym.

[n,~]=size(qnums);
[m,~]=size(qnumssym);

U=zeros(n,m);
for(i1=1:n)
for(i2=1:m)
	ja1=qnums(i1,5);
	jb1=qnums(i1,6);
	ja2=qnumssym(i2,5);
	jb2=qnumssym(i2,6);
	if(qnums(i1,1)==qnumssym(i2,1) && qnums(i1,2)==qnumssym(i2,2) && qnums(i1,3)==qnumssym(i2,3) && qnums(i1,4)==qnumssym(i2,4))
		if(ja1==ja2 && jb1==jb2)
			if(ja1==jb1)
				U(i1,i2)=1/2;
			else
				U(i1,i2)=sqrt(2)/2;
			end
		end
		if(ja1==jb2 && jb1==ja2)
			if(ja1==jb1)
				U(i1,i2)=U(i1,i2)+permsym*(-1)^(ja1+jb1-qnums(i1,3)+qnums(i1,4))*1/2;
			else
				U(i1,i2)=U(i1,i2)+permsym*(-1)^(ja1+jb1-qnums(i1,3)+qnums(i1,4))*sqrt(2)/2;
			end
		end
	end
end	
end

end
