function Q=numerov(Rgrid,E,mu,wmat_fun,wmat_vars)
//fixed step numerov propagator.
//
//npoints, delta the stepsize, E scattering energy. mu the reduced mass and wmat, the W-matrix at zero scattering energy.
//inputs all wmats (costs some memory, but i might be able to save them more cheaply as a cell array of sparse matrices)

	npoints=length(Rgrid);
	delta=Rgrid(2)-Rgrid(1);

	wmat=wmat_fun(Rgrid(1),wmat_vars);
	Q=zeros(wmat);
	Fmineen=eye(wmat)-delta^2/12*(wmat-2*mu*E*eye(wmat));
	wmat=wmat_fun(Rgrid(2),wmat_vars);
	F=eye(wmat)-delta^2/12*(wmat-2*mu*E*eye(wmat));

	for(i=3:npoints)
		Fmintwee=Fmineen;
		Fmineen=F;
		wmat=wmat_fun(Rgrid(i),wmat_vars);
		F=eye(wmat)-delta^2/12*(wmat-2*mu*E*eye(wmat));
		Q=(12*eye(wmat)-10*Fmineen-Fmintwee*Q)\F;
	end

endfunction

function [Q]=Qsin(Rgrid,E,mu,Q,wmat_fun,wmat_vars)
// comments go here.

	// initialize stuff.
	siz=size(Q(1));
	nchans=siz(1);
	nrgs=length(E);
	npoints=length(Rgrid);

	// loop grid
	for(iR=2:npoints-1)

		// get W-matrix except for scattering energy, at the midpoint of the two-point intervals
		dRp=Rgrid(iR+1)-Rgrid(iR);
		dRm=Rgrid(iR-1)-Rgrid(iR);
		wmatm=wmat_fun((Rgrid(iR)+Rgrid(iR-1))/2,wmat_vars);
		wmatp=wmat_fun((Rgrid(iR)+Rgrid(iR+1))/2,wmat_vars);

		// transform to locally adiabatic basis.
		[Un,Vp]=seig(wmatp);    //diagonalize W on the second interval.
		Vm=diag(Un'*wmatm*Un);  //and assume the same transformation diagonalizes W on the first interval.

		// loop successive energies.
		for(iE=1:nrgs)

			// energy dependent part of potential.
			Wm=Vm-2*mu*E(iE)*ones(Vm);
			Wp=Vp-2*mu*E(iE)*ones(Vp);

			// fix that scilab does not get that these quantities are real
			Wm=real(Wm);
			Wp=real(Wp);

			// get wronskians
			Wronp=-sqrt(abs(Wp));
			Wronm=-sqrt(abs(Wm));

			// two-point interval factors
			posm=Wm>0;
			posp=Wp>0;
			Xm(posm)=sinh(Wronm(posm)*dRm);
			Xp(posp)=sinh(Wronp(posp)*dRp);
			Ym(posm)=-Wronm(posm).*cosh(Wronm(posm)*dRm);
			Yp(posp)=-Wronp(posp).*cosh(Wronp(posp)*dRp);
			Xm(~posm)=sin(Wronm(~posm)*dRm);
			Xp(~posp)=sin(Wronp(~posp)*dRp);
			Ym(~posm)=-Wronm(~posm).*cos(Wronm(~posm)*dRm);
			Yp(~posp)=-Wronp(~posp).*cos(Wronp(~posp)*dRp);
			At=Wronm.*Xp./(Xp.*Ym-Xm.*Yp);
			Ct=Wronp.*Xm./(Xm.*Yp-Xp.*Ym);

			// transform to primitive basis 
			A=Un*diag(At)*Un';
			C=Un*diag(Ct)*Un';

			//propagate Q
			Q(iE)=-(A*Q(iE)+eye(nchans,nchans))\C;

		end //end loop energy
	end //end loop R

endfunction


function [Q]=Qairy(Rgrid,E,mu,Q,wmat_fun,wmat_vars)
// comments go here

	// initialize stuff.
	siz=size(Q(1));
	nchans=siz(1);
	nrgs=length(E);
	npoints=length(Rgrid);

	// loop grid
	for(iR=2:npoints-1)

		// get W-matrix except for scattering energy, at the midpoint of the two-point intervals
		dRp=Rgrid(iR)-Rgrid(iR+1);
		dRm=Rgrid(iR)-Rgrid(iR-1);
		wmatm=wmat_fun((Rgrid(iR)+Rgrid(iR-1))/2,wmat_vars);
		wmatp=wmat_fun((Rgrid(iR)+Rgrid(iR+1))/2,wmat_vars);

		// transform to locally adiabatic basis.
		wmat1=wmat_fun(Rgrid(iR),wmat_vars);
		[Un,wmat1]=seig(wmat1);     //diagonalize W at the center.
		wmat1=diag(Un'*wmat_fun((Rgrid(iR)+Rgrid(iR-1))/2-dRm/(2*sqrt(3)),wmat_vars)*Un);
		wmat2=diag(Un'*wmat_fun((Rgrid(iR)+Rgrid(iR-1))/2+dRm/(2*sqrt(3)),wmat_vars)*Un);
		W1m=(wmat2-wmat1)/((Rgrid(iR)-Rgrid(iR-1))/sqrt(3));
		W1m(abs(W1m)<%eps)=%eps;  //set slope to at least machine precission, to avoid dividing by zero
		V0m=(wmat1+wmat2+W1m*dRm)/2;
		wmat1=diag(Un'*wmat_fun((Rgrid(iR)+Rgrid(iR+1))/2-dRp/(2*sqrt(3)),wmat_vars)*Un);
		wmat2=diag(Un'*wmat_fun((Rgrid(iR)+Rgrid(iR+1))/2+dRp/(2*sqrt(3)),wmat_vars)*Un);
		W1p=(wmat2-wmat1)/((Rgrid(iR)-Rgrid(iR+1))/sqrt(3));
		W1p(abs(W1p)<%eps)=%eps;  //set slope to at least machine precission, to avoid dividing by zero
		V0p=(wmat1+wmat2+W1p*dRp)/2;

		// fix that scilab does not get that these quantities are real
		V0m=real(V0m);
		W1m=real(W1m);
		V0p=real(V0p);
		W1p=real(W1p);

		// loop successive energies.
		for(iE=1:nrgs)

			// energy dependent part of potential.
			W0m=V0m-2*mu*E(iE)*ones(nchans,1);
			W0p=V0p-2*mu*E(iE)*ones(nchans,1);

			// get wronskians
			Wronp=nthroot(W1p,3)/%pi;
			Wronm=nthroot(W1m,3)/%pi;

			// two-point interval factors
			for(i1=1:length(W0m))
				Xm=Xairy(W0m(i1),nthroot(W1m(i1),3),-dRm);
				Xp=Xairy(W0p(i1),nthroot(W1p(i1),3),-dRp);
				Ym=Yairy(W0m(i1),nthroot(W1m(i1),3),-dRm);
				Yp=Yairy(W0p(i1),nthroot(W1p(i1),3),-dRp);
				At(i1)=Wronm(i1)*Xp/(Xp*Ym-Xm*Yp);
				Ct(i1)=Wronp(i1)*Xm/(Xm*Yp-Xp*Ym);
			end

			// transform to primitive basis 
			A=Un*diag(At)*Un';
			C=Un*diag(Ct)*Un';

			//propagate Q
			Q(iE)=-(A*Q(iE)+eye(nchans,nchans))\C;

		end //end loop energy
	end //end loop R

endfunction
