function e = eenxkm(x,k)

i=1;
c=k;
e(i)=c*x;

while(abs(e(i)/e(1))>2*eps)
        c=c*(k-i)/(i+1);
	i=i+1;

	e(i)=c*x^i;
end

e=sum(e);

end
