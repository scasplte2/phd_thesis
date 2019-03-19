function w=ScScW2(R,params)
% W-matrix generating function for Scandium-Scandium.
%
%august 20.
	B=params{1}; %0 for fieldfree, contains fieldstrength otherwise.
        mu=params{2};
	Uperm=params{3};
	Vso=diag(params{4});
	Vmagndip=params{5};
        Vcent=diag(params{6});
        Vzeeman=params{7};
	Vint=params{8};
	ktabs=params{9};
	ctabs=params{10};
	Ukc=params{11};
	x=params{12};
	alpha=params{13};
	betas=params{14};
	b=params{15};
	c=params{16};

	[kpoints,~]=size(ktabs);
	Vskkk=potentialSpaceFixedExpansion2(x,alpha,betas,b,c,R,ktabs,ctabs,Ukc);
	Vintsummed=zeros(size(Vso));
	for(ik=1:kpoints)
                Vintsummed=Vintsummed+Vint(:,:,ik)*Vskkk(ik);
        end
	wtemp=sparse(Vcent/(R^2)+2*mu*(Vintsummed+Vso+B*Vzeeman+Vmagndip/(R^3)));
	w=Uperm'*wtemp*Uperm;
	w=full((w+w')/2);

end
