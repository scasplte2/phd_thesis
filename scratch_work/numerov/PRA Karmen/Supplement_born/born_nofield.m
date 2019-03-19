%born approximation crosssections for field-free case

%constants
mub=0.5;
alpha=1/137.035999074;
mu=(44.9559119/5.485799110e-4)/2;
gs=2.0023193043718;

%input
LMAX=12;
ja=1.5;
jb=ja;
ma=1.5;
mb=1.5;
%map=0.5
%mbp=0.5
%gj=0.8;
LA=2;
SA=1/2;
al=(ja*(ja+1)-SA*(SA+1)+LA*(LA+1))/(2*ja*(ja+1));
as=(ja*(ja+1)+SA*(SA+1)-LA*(LA+1))/(2*ja*(ja+1));
%al=(3/2*(3/2+1)-1/2*(1/2+1)+2*(2+1))/(2*3/2*(3/2+1));
%as=(3/2*(3/2+1)+1/2*(1/2+1)-2*(2+1))/(2*3/2*(3/2+1));
gj=al+gs*as;
parity=1;

%determine final states
mps=[];
for(i1=-ja:ja)
for(i2=-jb:jb)
if(~(i1==ma&&i2==mb))
mps=[mps; i1 i2];
end
end
end
%mps=[3/2 3/2] %for elastic (odd parity).
[n,~]=size(mps);

%actual calculation
for(i=1:n)
map=mps(i,1);
mbp=mps(i,2);

%auxilliary
dma=map-ma;
dmb=mbp-mb;
dm=dma+dmb;

summation=0;
for(l=0:LMAX)
for(lp=abs(l-2):(l+2))
for(ml=-l:l)
if((l~=0 || lp~= 0) && (-1)^l == parity && (-1)^lp == parity)
summation=summation+(2*l+1)*(2*lp+1)*(ff_3jm([1 1 2; dma dmb -dm])*ff_3jm([ja 1 ja; -map dma ma])*ff_3jm([jb 1 jb; -mbp dmb mb])*ff_3jm([l 2 lp; 0 0 0])*ff_3jm([l 2 lp; ml dm -(ml+dm)])*cos(pi/2*(lp-l))/((lp-l-1)*(lp-l+1)*(l+lp)*(l+lp+2)))^2;
end
end
end
end

%summation
sigma(i)=summation*1920*pi*mu^2*(mub*alpha*gj)^4*(2*ja+1)*(2*jb+1)*ja*(ja+1)*jb*(jb+1);

end

sum(sigma,2)
