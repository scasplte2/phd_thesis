function qnums = basis_coupledpablmax_J(Jslct,lmax,parity,permsym)
%function qnums = basis_coupledpablmax_J(Jslct,lmax,parity,permsym)
%
% generates permuation adapted basis for single Jtot, Jslct

qnums=[];
J=Jslct;
for(l=0:lmax)
for(j=abs(J-l):(J+l))
for(ja=3/2:5/2)
for(jb=3/2:ja)

if(j <= (ja+jb) && j >= abs(ja-jb) && (-1)^l==parity)
if(ja~=jb || (-1)^(ja+jb-j+l) == permsym)
        qnums=[qnums; J 0 j l ja jb];
end
end

end
end
end
end


end

