// calculates the X_i^\pm from the appendix
function X = Xairy(W0,W1,Delta)

	x=W0/(W1^2);
	Delta=W1*Delta;

	if(abs(Delta)>=abs(x)/2 | abs(x) <= 10) // use airy routines for small arguments or if delta not small compared to x.
	        X=real(airyai(x)*airybi(x+Delta)-airyai(x+Delta)*airybi(x));
	elseif x > 0    //else, use asymptotic expansion, which looks like this for possitive x
        	z=x;

	        sqrtz=sqrt(z);
        	zeta=2/3*sqrtz^3;
	        dz=Delta/x;

        	exparg=eenxkm(dz,3/2)*zeta;
	        X=exp(exparg)*airexp1(-zeta)*airexp1(exparg+zeta)-exp(-exparg)*airexp1(zeta)*airexp1(-exparg-zeta);

	        X=X/(2*%pi*sqrtz)*eenxk(dz,-1/4);

	else            //else, use asymptotic expansion, which looks like this for negative x
        	z=abs(x);

	        sqrtz=sqrt(z);
        	zeta=2/3*sqrtz^3;
	        dz=Delta/x;

        	trigarg=eenxkm(dz,3/2)*zeta;
	        X=-sin(trigarg)*(airexp1o(trigarg+zeta)*airexp1o(zeta)+airexp1e(zeta)*airexp1e(zeta+trigarg))+cos(trigarg)*(airexp1e(zeta)*airexp1o(zeta+trigarg)-airexp1e(trigarg+zeta)*airexp1o(zeta));

        	X=X/(%pi*sqrtz)*eenxk(dz,-1/4);
	end

endfunction

// calculates the Y_i^\pm from the appendix
function Y = Yairy(W0,W1,Delta)

	x=W0/(W1^2);
	Delta=W1*Delta;

	if(abs(Delta)>=abs(x)/2 | abs(x) <= 10) // use airy routines for small arguments or if delta not small compared to x.
        	Y=real(airyaip(x)*airybi(x+Delta)-airyai(x+Delta)*airybip(x));
	elseif x > 0    //else, use asymptotic expansion, which looks like this for possitive x
        	z=x;

	        sqrtz=sqrt(z);
        	zeta=2/3*sqrtz^3;
	        dz=Delta/x;

        	exparg=eenxkm(dz,3/2)*zeta;
	        Y=-exp(exparg)*airexp2(-zeta)*airexp1(exparg+zeta)-exp(-exparg)*airexp2(zeta)*airexp1(-exparg-zeta);

        	Y=Y*eenxk(dz,-1/4)/(2*%pi);
	else            //else, use asymptotic expansion, which looks like this for negative x
        	z=abs(x);

	        sqrtz=sqrt(z);
        	zeta=2/3*sqrtz^3;
	        dz=Delta/x;

	        trigarg=eenxkm(dz,3/2)*zeta;
        	Y=cos(trigarg)*(airexp2e(zeta)*airexp1e(trigarg+zeta)+airexp1o(zeta+trigarg)*airexp2o(zeta))+sin(trigarg)*(airexp1o(zeta+trigarg)*airexp2e(zeta)-airexp1e(trigarg+zeta)*airexp2o(zeta));

	        Y=-Y/%pi*eenxk(dz,-1/4);
	end

	Y=W1*Y;

endfunction

// The above two functions call the helper functions below.



function e = eenxkm(x,k)

	i=1;
	c=k;
	e(i)=c*x;

	while(abs(e(i)/e(1))>2*10*%eps)
        	c=c*(k-i)/(i+1);
	        i=i+1;

        	e(i)=c*x^i;
	end

	e=sum(e);

endfunction

function e = eenxk(x,k)

	e=eenxkm(x,k);
	e=e+1;

endfunction


function expans = airexp1(x)

	k=0;
	expans(1)=1;

	while(abs(expans(k+1)/expans(1))>2*10*%eps)
        	k=k+1;
	        c=gamma(3*k+1/2)/(factorial(k)*54^k*gamma(k+1/2));

        	expans(k+1)=c/(x)^k;
	end

	expans=sum(expans);

endfunction

function expans = airexp1e(x)

	k=0;
	expans(1)=1;

	while(abs(expans(k+1)/expans(1))>2*10*%eps)
        	k=k+1;
	        c=gamma(3*(2*k)+1/2)/(factorial(2*k)*54^(2*k)*gamma((2*k)+1/2));

        	expans(k+1)=(-1)^k*c/(x^(2*k));
	end

	expans=sum(expans);

endfunction

function expans = airexp1o(x)

	k=0;
	//c=gamma(3*(2*k+1)+1/2)/(factorial(2*k+1)*54^(2*k+1)*gamma((2*k+1)+1/2));
	c=3.75/54;
	expans(1)=c/x;

	while(abs(expans(k+1)/expans(1))>2*10*%eps)
        	k=k+1;
	        c=gamma(3*(2*k+1)+1/2)/(factorial(2*k+1)*54^(2*k+1)*gamma(2*k+1+1/2));

        	expans(k+1)=(-1)^k*c/(x^(2*k+1));
	end

	expans=sum(expans);

endfunction

function expans = airexp2(x)

	k=0;
	expans(1)=1;

	while(abs(expans(k+1)/expans(1))>2*10*%eps)
        	k=k+1;
	        c=gamma(3*k+1/2)/(factorial(k)*54^k*gamma(k+1/2));
        	d=-c*(6*k+1)/(6*k-1);

	        expans(k+1)=d/(x)^k;
	end

	expans=sum(expans);

endfunction

function expans = airexp2e(x)

	k=0;
	expans(1)=1;

	while(abs(expans(k+1)/expans(1))>2*10*%eps)
        	k=k+1;
	        c=gamma(3*(2*k)+1/2)/(factorial(2*k)*54^(2*k)*gamma((2*k)+1/2));
        	d=-c*(6*2*k+1)/(6*2*k-1);

	        expans(k+1)=(-1)^k*d/(x)^(2*k);
	end

	expans=sum(expans);

endfunction

function expans = airexp2o(x)

	k=0;
	//c=gamma(3*(2*k+1)+1/2)/(factorial(2*k+1)*54^(2*k+1)*gamma((2*k+1)+1/2));
	c=3.75/54;
	d=-c*7/5;
	expans(1)=d/x;

	while(abs(expans(k+1)/expans(1))>2*10*%eps)
        	k=k+1;
	        c=gamma(3*(2*k+1)+1/2)/(factorial(2*k+1)*54^(2*k+1)*gamma(2*k+1+1/2));
        	d=-c*(6*(2*k+1)+1)/(6*(2*k+1)-1);

	        expans(k+1)=(-1)^k*d/(x)^(2*k+1);
	end

	expans=sum(expans);

endfunction


function ai=airyai(x)

	if(x>0)
		ai=besselk(1/3,2/3*x^(3/2))*sqrt(x/3)/%pi;
	else
		x=abs(x)
		ai=sqrt(x)/3*(besselj(1/3,2/3*x^(3/2))+besselj(-1/3,2/3*x^(3/2)));
	end

endfunction

function bi=airybi(x)

	if(x>0)
        	bi=sqrt(x/3)*(besseli(1/3,2/3*x^(3/2))+besseli(-1/3,2/3*x^(3/2)));
	else
        	x=abs(x)
	        bi=sqrt(x/3)*(-besselj(1/3,2/3*x^(3/2))+besselj(-1/3,2/3*x^(3/2)));
	end

endfunction

function aip=airyaip(x)

	if(x>0)
        	aip=-x/(%pi*sqrt(3))*besselk(2/3,2/3*x^(3/2))
	else
        	x=abs(x)
	        aip=-x/3*(-besselj(2/3,2/3*x^(3/2))+besselj(-2/3,2/3*x^(3/2)));
	end

endfunction


function bip=airybip(x)

	if(x>0)
        	bip=x/sqrt(3)*(besseli(2/3,2/3*x^(3/2))+besseli(-2/3,2/3*x^(3/2)));
	else
        	x=abs(x)
	        bip=x/sqrt(3)*(besselj(2/3,2/3*x^(3/2))+besselj(-2/3,2/3*x^(3/2)));
	end

endfunction
