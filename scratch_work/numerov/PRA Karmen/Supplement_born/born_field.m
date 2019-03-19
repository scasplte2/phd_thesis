function sigma = born_field(Egrid,Bgrid)
%function sigma = born_field(Egrid,Bgrid)
%
%born approximation crosssections with magnetic fields

tic;
%constants
mub=0.5
alpha=1/137.035999074
mu=(44.9559119/5.485799110e-4)/2
%K=315775.04;
%gauss=4.254382549059841e-10;
gs=2.0023193043718;


%input
LMAX=12
ja=1.5
jb=1.5
ma=1.5
mb=1.5
LA=2;
SA=1/2;
al=(ja*(ja+1)-SA*(SA+1)+LA*(LA+1))/(2*ja*(ja+1));
as=(ja*(ja+1)+SA*(SA+1)-LA*(LA+1))/(2*JA*(JA+1));
gj=al+gs*as
parity=1

mps=[];
for(i1=-ja:ja)
for(i2=-jb:jb)
if(~(i1==ma&&i2==mb))
mps=[mps; i1 i2];
end
end
end
[n,~]=size(mps);

for(iE=1:length(Egrid))
for(iB=1:length(Bgrid))
E=Egrid(iE);
B=Bgrid(iB);

for(i=1:n)
map=mps(i,1);
mbp=mps(i,2);

%auxilliary
dma=map-ma;
dmb=mbp-mb;
dm=dma+dmb;

k=sqrt(2*mu*E);
kp=sqrt(2*mu*(E-mub*gj*(map+mbp-ma-mb)*B));

summation=0;
for(l=0:LMAX)
for(lp=abs(l-2):(l+2))
for(ml=-l:l)
if((l~=0 || lp~= 0) && (-1)^l == parity && (-1)^lp == parity)
summation=summation+(2*l+1)*(2*lp+1)*(k/kp)^(2*l-1)*(ff_3jm([1 1 2; dma dmb -dm])*ff_3jm([ja 1 ja; -map dma ma])*ff_3jm([jb 1 jb; -mbp dmb mb])*ff_3jm([l 2 lp; 0 0 0])*ff_3jm([l 2 lp; ml dm -(ml+dm)])*gamma((l+lp)/2)/(gamma((lp-l+3)/2)*gamma(l+3/2))*hypergeom([(l-lp-1)/2 (l+lp)/2],l+3/2,(k/kp)^2))^2;
end
end
end
end

tmp(i)=summation*15/2*pi^3*mu^2*(mub*alpha*gj)^4*(2*ja+1)*(2*jb+1)*ja*(ja+1)*jb*(jb+1);

end

toc
sigma(iE,iB)=sum(tmp)

end
end

end
