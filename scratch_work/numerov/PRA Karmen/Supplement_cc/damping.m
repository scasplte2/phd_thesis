function [y]=damping(x,n)
%y=damping(x,n)
%generalized tang-toennies damping for arbitrary n.
y=ones(size(x));%(ones takes care of the k=0 term)
for(k=1:n)
	y=y+x.^k/(factorial(k));
end
y=ones(size(x))-exp(-x).*y;

end
