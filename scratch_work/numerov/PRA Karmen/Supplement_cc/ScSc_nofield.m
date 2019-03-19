function [totalElastic,totalInelastic,Egrid,lmax,permsym,parity,Ms]=ScSc_nofield(lmax,parity,permsym,Ebegin,Eeind,Estappen,Rend)
% fieldfree calculation. propagation is carried out per J, and finally we sum all the relevant submatrices.

tic;

%load stuff.
mu=(44.9559119/5.485799110e-4)/2;
K=315775.04;
load output_fitShortrange2.mat
toc

%energy (/field) grids
if(Estappen==1)
        Egrid=Ebegin/K;
else
        Egrid=10.^(log10(Ebegin):(log10(Eeind)-log10(Ebegin))/(Estappen-1):log10(Eeind))/K;
end
B=0;

%shift zeropoint of energy.
A=168.34/219474.63137098*(8/20);       %experimental atomic spin-orbit splitting
E0=-3*A; %B=0; 

%echo input.
parity
permsym	% permsym corresponds to the total permutation symmetry, not the nuclear permutation symmetry defined in the paper. As derived in the paper, the difference is a factor (-1)^n_e=-1 for scandium
lmax
toc;

%radial grid.
R0=1;
stepsize=1e-2;
R1=20;
Rgrid1=R0:stepsize:R1;
npoints1=length(Rgrid1)
Rgrid2 = optLRgrid2(Rgrid1(npoints1-1:npoints1),200,50,Rend);
npoints2=length(Rgrid2)
toc

for(J=max(0,3-lmax):(lmax+5))
disp(strcat('Running calculation for J=',num2str(J)))
qnums{J+1}=basis_coupledlmax_J(J,lmax,parity);
qnumssym{J+1}=basis_coupledpablmax_J(J,lmax,parity,permsym);
size(qnums{J+1})
size(qnumssym{J+1})
toc

%interaction matrices.
Vso=diag(spin_orbit(qnums{J+1}));
Vmagndip=magndip(qnums{J+1});
[Vint,ktabs]=potentialCoupledMatrixElements(qnums{J+1});
ctabs=generate_ctabs();
Ukc=transformAdiabats2SFExpansion(ktabs,ctabs);
Vcent=diag(centrifugal(qnums{J+1}));
Uperm=adapt2PermutationSymmetry(qnums{J+1},qnumssym{J+1},permsym);
disp(Uperm)
Uperm=sparse(Uperm);
disp(Uperm)
disp('Done with angular matrices');
toc;

%wmat_parameters.
wmat_vars1{1}=0;	%magnetic field.
wmat_vars1{2}=mu;
wmat_vars1{3}=Uperm;
wmat_vars1{4}=Vso;
wmat_vars1{5}=Vmagndip;
wmat_vars1{6}=Vcent;
wmat_vars1{7}=zeros(size(Vmagndip)); %Vzeeman.
wmat_vars1{8}=Vint;
wmat_vars1{9}=ktabs;
wmat_vars1{10}=ctabs;
wmat_vars1{11}=Ukc;
wmat_vars1{12}=x;
wmat_vars1{13}=alpha2;
wmat_vars1{14}=betas;
wmat_vars1{15}=b;
wmat_vars1{16}=c;
disp('Done setting up W matrix parameters');
toc;

%loop energy, and do all numerov propagation steps.
for(iE=1:length(Egrid))
E=Egrid(iE)+E0;
Qntmp=numerov6(Rgrid1,E,mu,@ScScW2,wmat_vars1);
Qn{iE}=Qntmp;
end
clear Qntmp;
disp('Done with Numerov');
toc;

%propagate longrange parallel (to save time on diagonalizations)
Qn=Qairy(Rgrid2,Egrid+E0,mu,Qn,@ScScW2,wmat_vars1);
disp('Done with Airy');
toc;

%loop energy again, and match to boundary conditions.
Has{J+1}=Uperm'*diag(Vso)*Uperm;
for(iE=1:length(Egrid))
E=Egrid(iE)+E0;

[S,Kopen]=matchBoundaryConditions(Rgrid2,qnumssym{J+1}(:,4),mu,Qn{iE},E,diag(Has{J+1}));
T{J+1,iE}=eye(size(S))-S;
end
disp('Done matching');
toc;

end %loop J.
disp('Continue with crosssections');
toc;

%finally the cross sections, as a function of M and E.
Ms=(3-lmax):(3+lmax);
for(iM=1:length(Ms))
M=Ms(iM);
for(iE=1:length(Egrid))
E=Egrid(iE)+E0;

%total basis (contributions of all the J's).
qnumstot=[];
Ttot=[];
Hastot=[];
for(J=abs(M):lmax+5)
qnumstot=[qnumstot; qnumssym{J+1}];
Ttot=blkdiag(Ttot,T{J+1,iE});
Hastot=blkdiag(Hastot,Has{J+1});
end
qnumstot(:,2)=M;

%decouple.
[Udec,qnumsunc]=decouplePermutationAdapted(qnumstot,permsym);
Hasunc=Udec'*Hastot*Udec;
Echantot=diag(Hastot);
Echanunc=diag(Hasunc);
Tunc=Udec(Echantot<E,Echanunc<E)'*Ttot*Udec(Echantot<E,Echanunc<E); %transform the open-open block only.

%actual cross sections for zeeman relaxation.
[totalElastic(iE,iM),totalInelastic(iE,iM)]=crosssectionsUncoupled(Tunc,qnumsunc(Echanunc<E,:),E,Echanunc(Echanunc<E,:),mu);

end%E.
end%M.
disp('All done');
toc

end
