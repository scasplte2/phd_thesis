function qnums = basis_coupledlmax(lmax,M,parity)
%function qnums = basis_coupledlmax(lmax,M,parity)
%
% generates coupled basis containing multiple Jtot

qnums=[];

for(J=abs(M):lmax+5)
for(l=0:lmax)
for(j=abs(J-l):(J+l))
for(ja=3/2:5/2)
for(jb=3/2:5/2)

if(j <= (ja+jb) && j >= abs(ja-jb) && (-1)^l==parity)
        qnums=[qnums; J M j l ja jb];
end

end
end
end
end
end

end

