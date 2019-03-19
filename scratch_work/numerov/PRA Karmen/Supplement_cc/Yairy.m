function Y = Yairy(W0,W1,Delta)

x=W0/(W1^2);
Delta=W1*Delta;

if(abs(Delta)>=abs(x)/2 || abs(x) <= 10) % use airy routines for small arguments or if delta not small compared to x.
	Y=real(airy(1,x)*airy(2,x+Delta)-airy(0,x+Delta)*airy(3,x));
elseif x > 0	%else, use asymptotic expansion, which looks like this for possitive x
	z=x;

	sqrtz=sqrt(z);
	zeta=2/3*sqrtz^3;
	dz=Delta/x;

	exparg=eenxkm(dz,3/2)*zeta;
	Y=-exp(exparg)*airexp2(-zeta)*airexp1(exparg+zeta)-exp(-exparg)*airexp2(zeta)*airexp1(-exparg-zeta);

	Y=Y*eenxk(dz,-1/4)/(2*pi);
else		%else, use asymptotic expansion, which looks like this for negative x
	z=abs(x);

        sqrtz=sqrt(z);
        zeta=2/3*sqrtz^3;
        dz=Delta/x;

        trigarg=eenxkm(dz,3/2)*zeta;
        Y=cos(trigarg)*(airexp2e(zeta)*airexp1e(trigarg+zeta)+airexp1o(zeta+trigarg)*airexp2o(zeta))+sin(trigarg)*(airexp1o(zeta+trigarg)*airexp2e(zeta)-airexp1e(trigarg+zeta)*airexp2o(zeta));

        Y=-Y/pi*eenxk(dz,-1/4);
end

Y=W1*Y;

end
