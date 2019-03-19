function X = Xairy(W0,W1,Delta)

x=W0/(W1^2);
Delta=W1*Delta;

if(abs(Delta)>=abs(x)/2 || abs(x) <= 10) % use airy routines for small arguments or if delta not small compared to x.
	X=real(airy(0,x)*airy(2,x+Delta)-airy(0,x+Delta)*airy(2,x));
elseif x > 0	%else, use asymptotic expansion, which looks like this for possitive x
	z=x;

	sqrtz=sqrt(z);
	zeta=2/3*sqrtz^3;
	dz=Delta/x;

	exparg=eenxkm(dz,3/2)*zeta;
	X=exp(exparg)*airexp1(-zeta)*airexp1(exparg+zeta)-exp(-exparg)*airexp1(zeta)*airexp1(-exparg-zeta);
	
	X=X/(2*pi*sqrtz)*eenxk(dz,-1/4);

else		%else, use asymptotic expansion, which looks like this for negative x
	z=abs(x);

	sqrtz=sqrt(z);
        zeta=2/3*sqrtz^3;
        dz=Delta/x;

        trigarg=eenxkm(dz,3/2)*zeta;
	X=-sin(trigarg)*(airexp1o(trigarg+zeta)*airexp1o(zeta)+airexp1e(zeta)*airexp1e(zeta+trigarg))+cos(trigarg)*(airexp1e(zeta)*airexp1o(zeta+trigarg)-airexp1e(trigarg+zeta)*airexp1o(zeta));

	X=X/(pi*sqrtz)*eenxk(dz,-1/4);
end

end
