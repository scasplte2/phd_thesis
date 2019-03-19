function []=generate_18j()
%tkarman 1 maart 2012
%since calculating the 18-j symbols is a bottleneck, this script will generate all 18-js i will possibly need.
%symmetry properties are used to speed up the process a little.
%this will be somewhat messy.
%the biggest gain is really in the fact that the 18-j s dont depend on the total angular momentum, nor on the end over end rotation, so that we dont have to repeat the same thing over-and-over by tabulating  the relevant values once and for all.
%we return nothing, but save a struct eightteenj in eightteenj.mat 
%we assume La=Lb=Lap=Lbp=2 and same for S=1/2. (scandium, or any other ^2D_g)
%eightteenj{s,k1,k2,k,ja,jb,jap,jbp,j,jp}

La=2;
Lb=2;
Sa=1/2;
Sb=1/2;
Lap=2;
Lbp=2;
Sap=1/2;
Sbp=1/2;

eightteenj=zeros(2,5,5,9,2,2,2,2,6,6);

for(S=abs(Sa-Sb):(Sa+Sb))
for(k1=0:(2*La))
for(k2=0:(2*Lb)) %for(k2=k1:(2*Lb))
for(k=0:2:2*(La+Lb))
for(ja=abs(La-Sa):(La+Sa))
for(jap=abs(Lap-Sap):(Lap+Sap))
for(jb=abs(Lb-Sb):(Lb+Sb))
for(jbp=abs(Lbp-Sbp):(Lbp+Sbp))
for(j=abs(ja-jb):(ja+jb))
for(jp=abs(jap-jbp):(jap+jbp))

if(k<abs(k1-k2) || k > k1+k2) %special cases in which the eighteen j is zero

sommatie=0;

else	%when it isnt nescesarily zero.

sommatie=0;
for(f=abs(La-Lb):(La+Lb))
for(fp=abs(Lap-Lbp):(Lap+Lbp))
sommatie=sommatie+(2*f+1)*(2*fp+1)*(-1)^(k1+k2+Lap+Lb+Sb+Sap-jb-jap+j+jp)*ff_6j([j jp k; fp f S])*ff_9j([La Lap k1; Lb Lbp k2; f fp k])*ff_9j([Sa La ja; Sb Lb jb; S f j])*ff_9j([Sbp Lbp jbp; Sap Lap jap; S fp jp]);
end
end

end
eightteenj(S+1,k1+1,k2+1,k+1,ja-1/2,jb-1/2,jap-1/2,jbp-1/2,j+1,jp+1)=sommatie;
%eightteenj(S+1,k2+1,k1+1,k+1,jbp-1/2,jap-1/2,jb-1/2,ja-1/2,jp+1,j+1)=sommatie;
%cell array is slower, somehow.
%eightteenj{2*S+1,k1+1,k2+1,k+1,2*ja+1,2*jb+1,2*j+1,2*jap+1,2*jbp+1,2*jp+1}=sommatie;
%eightteenj{2*S+1,k2+1,k1+1,k+1,2*jbp+1,2*jap+1,2*jp+1,2*jb+1,2*ja+1,2*j+1}=sommatie;

end
end
end
end
end
end
end
end
end
end

save eightteenj.mat eightteenj

end
