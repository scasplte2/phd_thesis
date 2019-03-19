function [totalElastic,totalInelastic,Egrid,B]=ScSc_field(lmax,M,parity,permsym,B,Ebegin,Eeind,Estappen)
%
% Calculation with magnetic field.
% Loop energy, one field strength
%

tic;

%load stuff.
mu=(44.9559119/5.485799110e-4)/2;
K=315775.04;
gauss=4.254382549059841e-10;
load output_fitShortrange2.mat
toc

%energy (/field) grids
if(Estappen==1)
	Egrid=Ebegin/K;
else
	Egrid=10.^(log10(Ebegin):(log10(Eeind)-log10(Ebegin))/(Estappen-1):log10(Eeind))/K;
end
B=B*gauss;

%shift zeropoint of energy.
%This will be done numerically below.

%define angular basis.
parity
permsym		% permsym corresponds to the total permutation symmetry, not the nuclear permutation symmetry defined in the paper. As derived in the paper, the difference is a factor (-1)^n_e=-1 for scandium
lmax
M
qnums=basis_coupledlmax(lmax,M,parity);
qnumssym=basis_coupledpablmax(lmax,M,parity,permsym);
size(qnums)
size(qnumssym)
toc;

%interaction matrices.
toc;
Vso=diag(spin_orbit(qnums));
toc;
Vmagndip=magndip(qnums);
toc;
[Vint,ktabs]=potentialCoupledMatrixElements(qnums);
toc;
ctabs=generate_ctabs();
Ukc=transformAdiabats2SFExpansion(ktabs,ctabs);
toc;
Vcent=diag(centrifugal(qnums));
toc;
Uperm=adapt2PermutationSymmetry(qnums,qnumssym,permsym);
Uperm=sparse(Uperm);
toc;
Vzeeman=real(zeeman(qnums));
toc;

%wmat_parameters.
wmat_vars1{1}=B;
wmat_vars1{2}=mu;
wmat_vars1{3}=Uperm;
wmat_vars1{4}=Vso;
wmat_vars1{5}=Vmagndip;
wmat_vars1{6}=Vcent;
wmat_vars1{7}=Vzeeman;
wmat_vars1{8}=Vint;
wmat_vars1{9}=ktabs;
wmat_vars1{10}=ctabs;
wmat_vars1{11}=Ukc;
wmat_vars1{12}=x;
wmat_vars1{13}=alpha2;
wmat_vars1{14}=betas;
wmat_vars1{15}=b;
wmat_vars1{16}=c;

%radial grid.
R0=1;
stepsize=1e-2;
R1=20;
Rgrid1=R0:stepsize:R1;
npoints1=length(Rgrid1)
n=4;
Rend=5e4;
Rgrid2 = optLRgrid2(Rgrid1(npoints1-1:npoints1),200,50,Rend);
npoints2=length(Rgrid2)
toc

%decouple Has.
[Udec,qnumsunc]=decouplePermutationAdapted(qnumssym,permsym);
Has=Uperm'*(diag(Vso)+B*Vzeeman)*Uperm;
Hasunc=Udec'*Has*Udec;
%now determine the transformation ja->na jb->nb
for(mja=(-5/2):(5/2))
        Uatoms{mja+7/2}=numDiagHasAtomic(mja,B);
end
Udiag=numDiagHas(qnumsunc,Uatoms,permsym);
Hasdiag=Udiag'*Hasunc*Udiag;
Echan=diag(Hasdiag);
qnumsdiag=qnumsunc;
%determine the zeropoint of energy (energy of entrance channels).
E0=Echan(and(qnumsunc(:,1)==1.5,and(qnumsunc(:,2)==1.5,and(qnumsunc(:,3)==1.5,qnumsunc(:,4)==1.5))));
if(max(abs(E0-mean(E0)))/mean(E0)>3*1.11e-16)
        disp('Warning: ambiguity in zeropoint of energy');
end
E0=mean(E0);

%loop energy, and do all numerov propagation steps.
for(iE=1:length(Egrid))
toc;
E=Egrid(iE)+E0;
Qntmp=numerov6(Rgrid1,E,mu,@ScScW2,wmat_vars1);
Qn{iE}=Qntmp;
end
clear Qntmp;

%propagate longrange parallel (to save time on diagonalizations)
toc;
Qn=Qairy(Rgrid2,Egrid+E0,mu,Qn,@ScScW2,wmat_vars1);

%loop energy again, and match to boundary conditions for crosssections.
for(iE=1:length(Egrid))
toc;
E=Egrid(iE)+E0;
Qnunc=Udec'*Qn{iE}*Udec;
Qndiag=Udiag'*Qnunc*Udiag;
[S,Kopen]=matchBoundaryConditions(Rgrid2,qnumsunc(:,5),mu,Qndiag,E,Echan);
T=eye(size(S))-S;
[totalElastic(iE),totalInelastic(iE)]=crosssectionsUncoupled(T,qnumsunc(Echan<E,:),E,Echan(Echan<E,:),mu);
end

end
