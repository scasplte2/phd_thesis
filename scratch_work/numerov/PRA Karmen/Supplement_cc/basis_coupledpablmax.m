function qnums = basis_coupledpablmax(lmax,M,parity,permsym)
%function qnums = basis_coupledpablmax(lmax,M,parity,permsym)
%
% Generates permutation symmetry adapted basis, containing multiple Jtot

qnums=[];
for(J=abs(M):(lmax+5))
for(l=0:lmax)
for(j=abs(J-l):(J+l))
for(ja=3/2:5/2)
for(jb=3/2:ja)

if(j <= (ja+jb) && j >= abs(ja-jb) && (-1)^l==parity)
if(ja~=jb || (-1)^(ja+jb-j+l) == permsym)
        qnums=[qnums; J M j l ja jb];
end
end

end
end
end
end


end

