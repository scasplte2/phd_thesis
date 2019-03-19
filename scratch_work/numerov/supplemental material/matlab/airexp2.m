function expans = airexp2(x)

k=0;
expans(1)=1;

while(abs(expans(k+1)/expans(1))>2*eps)
	k=k+1;
	c=gamma(3*k+1/2)/(factorial(k)*54^k*gamma(k+1/2));
	d=-c*(6*k+1)/(6*k-1);

	expans(k+1)=d/(x)^k;
end

expans=sum(expans);

end
