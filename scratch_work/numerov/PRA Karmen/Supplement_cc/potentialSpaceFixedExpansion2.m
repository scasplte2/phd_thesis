function [Vskkk]=potentialSpaceFixedExpansion2(x,alpha,betas,b,c,R,ktabs,ctabs,Ukc)
%tijs karman 12 feb 2012
%transforms the effective potential form as given by effectivePotential to the space fixed expansion coefficients (eq. 19 jpcA 108 8943)
%x=[c5 c6000 c6220 c6202 c6222 c6224 c6110 c6112];
%
%differs from version one in the sense that it transforms with a matrix, instead of with lots of sums...
%this should be faster to allow for evaluation at arbitrary R, instead of having to do the entire grid
%simultaneously, in order to keep it feasible.

%evaluate effective potential due to morse.
for(i=1:30)
        Emorse(i)=morse_eval(R,alpha(i),0,b(i),c(i));
end

%order
Vsclambda(1,1,1)=Emorse(4);
Vsclambda(1,2,1)=Emorse(1);
Vsclambda(1,3,1)=Emorse(2);
Vsclambda(1,4,1)=Emorse(5);
Vsclambda(1,5,1)=Emorse(3);
Vsclambda(1,1,2)=Emorse(11);
Vsclambda(1,2,2)=Emorse(13);
Vsclambda(1,3,2)=Emorse(12);
Vsclambda(1,4,2)=Emorse(14);
Vsclambda(1,1,3)=Emorse(21);
Vsclambda(1,2,3)=Emorse(19);
Vsclambda(1,3,3)=Emorse(20);
Vsclambda(1,1,4)=Emorse(25);
Vsclambda(1,2,4)=Emorse(26);
Vsclambda(1,1,5)=Emorse(29);
Vsclambda(2,1,1)=Emorse(6);
Vsclambda(2,2,1)=Emorse(8);
Vsclambda(2,3,1)=Emorse(9);
Vsclambda(2,4,1)=Emorse(7);
Vsclambda(2,5,1)=Emorse(10);
Vsclambda(2,1,2)=Emorse(17);
Vsclambda(2,2,2)=Emorse(15);
Vsclambda(2,3,2)=Emorse(18);
Vsclambda(2,4,2)=Emorse(16);
Vsclambda(2,1,3)=Emorse(22);
Vsclambda(2,2,3)=Emorse(23);
Vsclambda(2,3,3)=Emorse(24);
Vsclambda(2,1,4)=Emorse(28);
Vsclambda(2,2,4)=Emorse(27);
Vsclambda(2,1,5)=Emorse(30);

%transform
[nc,~]=size(ctabs);
for(ic=1:nc)
Vc(ic,1)=Vsclambda(ctabs(ic,1)+1,ctabs(ic,2),abs(ctabs(ic,3))+1);
end

Vskkk=Ukc*Vc;

%find long range contributions.
k000=find(and(ktabs(:,2)==0,and(ktabs(:,3)==0,ktabs(:,4)==0)));
k110=find(and(ktabs(:,2)==1,and(ktabs(:,3)==1,ktabs(:,4)==0)));
k112=find(and(ktabs(:,2)==1,and(ktabs(:,3)==1,ktabs(:,4)==2)));
k220=find(and(ktabs(:,2)==2,and(ktabs(:,3)==2,ktabs(:,4)==0)));
k222=find(and(ktabs(:,2)==2,and(ktabs(:,3)==2,ktabs(:,4)==2)));
k224=find(and(ktabs(:,2)==2,and(ktabs(:,3)==2,ktabs(:,4)==4)));
k202=find(and(ktabs(:,2)==2,and(ktabs(:,3)==0,ktabs(:,4)==2)));
k022=find(and(ktabs(:,2)==0,and(ktabs(:,3)==2,ktabs(:,4)==2)));

%and add.
Vskkk(k000)=Vskkk(k000)-x(1)./(R.^6).*damping(R*betas(1),6);
Vskkk(k220)=Vskkk(k220)-x(2)./(R.^6).*damping(R*betas(1),6);
Vskkk(k202)=Vskkk(k202)-x(4)./(R.^6).*damping(R*betas(2),6);
Vskkk(k022)=Vskkk(k022)-x(4)./(R.^6).*damping(R*betas(2),6);
Vskkk(k222)=Vskkk(k222)-x(6)./(R.^6).*damping(R*betas(2),6);
Vskkk(k110)=Vskkk(k110)-x(3)./(R.^6).*damping(R*betas(1),6);
Vskkk(k112)=Vskkk(k112)-x(7)./(R.^6).*damping(R*betas(2),6);
Vskkk(k224)=Vskkk(k224)+x(9)./(R.^5).*damping(R*betas(3),5)-x(8)./(R.^6).*damping(R*betas(2),6);

%Vskkk(k000)=Vskkk(k000)-x(2)./(R.^6).*damping(R*betas(1),6);
%Vskkk(k220)=Vskkk(k220)-x(3)./(R.^6).*damping(R*betas(1),6);
%Vskkk(k202)=Vskkk(k202)-x(4)./(R.^6).*damping(R*betas(2),6);
%Vskkk(k022)=Vskkk(k022)-x(4)./(R.^6).*damping(R*betas(2),6);
%Vskkk(k222)=Vskkk(k222)-x(5)./(R.^6).*damping(R*betas(2),6);
%Vskkk(k110)=Vskkk(k110)-x(7)./(R.^6).*damping(R*betas(1),6);
%Vskkk(k112)=Vskkk(k112)-x(8)./(R.^6).*damping(R*betas(2),6);
%Vskkk(k224)=Vskkk(k224)+x(1)./(R.^5).*damping(R*betas(3),5)-x(6)./(R.^6).*damping(R*betas(2),6);

end
