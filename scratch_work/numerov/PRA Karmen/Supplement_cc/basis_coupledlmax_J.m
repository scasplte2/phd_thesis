function qnums = basis_coupledlmax_J(Jslct,lmax,parity)
%function qnums = basis_coupledlmax_J(Jslct,lmax,parity)
%generates basis for a single selected Jtot, Jslct
qnums=[];

J=Jslct; 
for(l=0:lmax)
for(j=abs(J-l):(J+l))
for(ja=3/2:5/2)
for(jb=3/2:5/2)

if(j <= (ja+jb) && j >= abs(ja-jb) && (-1)^l==parity)
        qnums=[qnums; J 0 j l ja jb];
end

end
end
end
end

end

