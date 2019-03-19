function expans = airexp1e(x)

k=0;
expans(1)=1;

while(abs(expans(k+1)/expans(1))>2*eps)
	k=k+1;
	c=gamma(3*(2*k)+1/2)/(factorial(2*k)*54^(2*k)*gamma((2*k)+1/2));

	expans(k+1)=(-1)^k*c/(x^(2*k));
end

expans=sum(expans);

end
