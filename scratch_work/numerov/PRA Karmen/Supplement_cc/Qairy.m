function [Q]=Qairy(Rgrid,E,mu,Q,wmat_fun,wmat_vars)
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
	dRp=Rgrid(iR)-Rgrid(iR+1);
	dRm=Rgrid(iR)-Rgrid(iR-1);
	wmatm=feval(wmat_fun,(Rgrid(iR)+Rgrid(iR-1))/2,wmat_vars);
	wmatp=feval(wmat_fun,(Rgrid(iR)+Rgrid(iR+1))/2,wmat_vars);

	% transform to locally adiabatic basis.
	wmat1=feval(wmat_fun,Rgrid(iR),wmat_vars);
	[Un,~]=seig(wmat1);	%diagonalize W at the center.
	wmat1=diag(Un'*feval(wmat_fun,(Rgrid(iR)+Rgrid(iR-1))/2-dRm/(2*sqrt(3)),wmat_vars)*Un);
	wmat2=diag(Un'*feval(wmat_fun,(Rgrid(iR)+Rgrid(iR-1))/2+dRm/(2*sqrt(3)),wmat_vars)*Un);
	W1m=(wmat2-wmat1)/((Rgrid(iR)-Rgrid(iR-1))/sqrt(3));
	W1m(abs(W1m)<eps)=eps;	%set slope to at least machine precission, to avoid dividing by zero
	V0m=(wmat1+wmat2+W1m*dRm)/2;
	wmat1=diag(Un'*feval(wmat_fun,(Rgrid(iR)+Rgrid(iR+1))/2-dRp/(2*sqrt(3)),wmat_vars)*Un);
	wmat2=diag(Un'*feval(wmat_fun,(Rgrid(iR)+Rgrid(iR+1))/2+dRp/(2*sqrt(3)),wmat_vars)*Un);
	W1p=(wmat2-wmat1)/((Rgrid(iR)-Rgrid(iR+1))/sqrt(3));
	W1p(abs(W1p)<eps)=eps;	%set slope to at least machine precission, to avoid dividing by zero
	V0p=(wmat1+wmat2+W1p*dRp)/2;

	% loop successive energies.
	for(iE=1:nrgs)

		% energy dependent part of potential.
		W0m=V0m-2*mu*E(iE)*ones(nchans,1);
		W0p=V0p-2*mu*E(iE)*ones(nchans,1);

		% get wronskians
		Wronp=nthroot(W1p,3)/pi;
		Wronm=nthroot(W1m,3)/pi;

		% two-point interval factors
		for(i1=1:length(W0m))
			Xm=Xairy(W0m(i1),nthroot(W1m(i1),3),-dRm);
			Xp=Xairy(W0p(i1),nthroot(W1p(i1),3),-dRp);
			Ym=Yairy(W0m(i1),nthroot(W1m(i1),3),-dRm);
			Yp=Yairy(W0p(i1),nthroot(W1p(i1),3),-dRp);
			At(i1)=Wronm(i1)*Xp/(Xp*Ym-Xm*Yp);
			Ct(i1)=Wronp(i1)*Xm/(Xm*Yp-Xp*Ym);
		end

		% transform to primitive basis 
		A=Un*diag(At)*Un';
		C=Un*diag(Ct)*Un';

		%propagate Q
		Q{iE}=-(A*Q{iE}+eye(nchans))\C;

	end %end loop energy
end %end loop R

end

