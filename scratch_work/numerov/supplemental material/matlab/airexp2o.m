function expans = airexp2o(x)

k=0;
%c=gamma(3*(2*k+1)+1/2)/(factorial(2*k+1)*54^(2*k+1)*gamma((2*k+1)+1/2));
c=3.75/54;
d=-c*7/5;
expans(1)=d/x;

while(abs(expans(k+1)/expans(1))>2*eps)
	k=k+1;
	c=gamma(3*(2*k+1)+1/2)/(factorial(2*k+1)*54^(2*k+1)*gamma(2*k+1+1/2));
	d=-c*(6*(2*k+1)+1)/(6*(2*k+1)-1);

	expans(k+1)=(-1)^k*d/(x)^(2*k+1);
end

expans=sum(expans);

end
