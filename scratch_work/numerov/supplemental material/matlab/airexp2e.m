function expans = airexp2e(x)

k=0;
expans(1)=1;

while(abs(expans(k+1)/expans(1))>2*eps)
	k=k+1;
	c=gamma(3*(2*k)+1/2)/(factorial(2*k)*54^(2*k)*gamma((2*k)+1/2));
	d=-c*(6*2*k+1)/(6*2*k-1);

	expans(k+1)=(-1)^k*d/(x)^(2*k);
end

expans=sum(expans);

end
