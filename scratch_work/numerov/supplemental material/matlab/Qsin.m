function [Q]=Qsin(Rgrid,E,mu,Q,wmat_fun,wmat_vars)
%
% tkarman 18 mei 2012
% variable stepsize sine/cosine propagator.
%function Q=Qsin(Rgrid,E,mu,wmat,Q)
%
%Rgrid ranges from the last two points done in Numerov, to the
%last done in sines and cosines.
%E is a vector of energies and Q a struct of Q matrices (accordingly).

% initialize stuff.
siz=size(Q{1});
nchans=siz(1);
nrgs=length(E);
npoints=length(Rgrid);

% loop grid
for(iR=2:npoints-1)

	% get W-matrix except for scattering energy, at the midpoint of the two-point intervals
	dRp=Rgrid(iR+1)-Rgrid(iR);
	dRm=Rgrid(iR-1)-Rgrid(iR);
	wmatm=feval(wmat_fun,(Rgrid(iR)+Rgrid(iR-1))/2,wmat_vars);
	wmatp=feval(wmat_fun,(Rgrid(iR)+Rgrid(iR+1))/2,wmat_vars);

	% transform to locally adiabatic basis.
	[Un,Vp]=seig(wmatp);	%diagonalize W on the second interval.
	Vm=diag(Un'*wmatm*Un);	%and assume the same transformation diagonalizes W on the first interval.

	% loop successive energies.
	for(iE=1:nrgs)

		% energy dependent part of potential.
		Wm=Vm-2*mu*E(iE)*ones(nchans,1);
		Wp=Vp-2*mu*E(iE)*ones(nchans,1);

		% get wronskians
		Wronp=-sqrt(abs(Wp));
		Wronm=-sqrt(abs(Wm));

		% two-point interval factors
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
		At=Wronm'.*Xp./(Xp.*Ym-Xm.*Yp);
		Ct=Wronp'.*Xm./(Xm.*Yp-Xp.*Ym);

		% transform to primitive basis 
		A=Un*diag(At)*Un';
		C=Un*diag(Ct)*Un';


		%propagate Q
		Q{iE}=-(A*Q{iE}+eye(nchans))\C;

	end %end loop energy
end %end loop R

end

