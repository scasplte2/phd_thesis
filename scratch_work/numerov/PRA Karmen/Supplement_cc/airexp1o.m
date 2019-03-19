function expans = airexp1o(x)

k=0;
%c=gamma(3*(2*k+1)+1/2)/(factorial(2*k+1)*54^(2*k+1)*gamma((2*k+1)+1/2));
c=3.75/54;
expans(1)=c/x;

while(abs(expans(k+1)/expans(1))>2*eps)
	k=k+1;
	c=gamma(3*(2*k+1)+1/2)/(factorial(2*k+1)*54^(2*k+1)*gamma(2*k+1+1/2));

	expans(k+1)=(-1)^k*c/(x^(2*k+1));
end

expans=sum(expans);

end
