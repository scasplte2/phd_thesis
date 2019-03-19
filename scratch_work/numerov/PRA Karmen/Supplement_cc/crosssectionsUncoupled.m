function [totalElastic,totalInelastic] = crosssectionsUncoupled(T,qnums,E,Echan,mu)
%tkarman 21 maart 2012
%obtain crosssections from the T matrix in uncoupled basis
%yields total elastic and total inelastic
%
%function [totalElastic,totalInelastic] = crosssectionsUncoupled(T,qnums,E,Echan,mu)

%entrance channels are those in the lowest spin orbit manifold, with maximally stretched projections |1.5 1.5>|1.5 1.5>|l ml>
[n,~]=size(T);
entrance=and(and(and(qnums(:,1)==1.5,qnums(:,2)==1.5),qnums(:,3)==1.5),qnums(:,4)==1.5);
k=sqrt(2*mu*(E-Echan));

totalElastic=0;
totalInelastic=0;
for(i1=1:n)
for(i2=1:n)

	if(entrance(i1))
	if(entrance(i2))	%elastic.
		totalElastic=totalElastic+abs(T(i2,i1)/k(i1))^2;
	else		%inelastic.
		totalInelastic=totalInelastic+abs(T(i2,i1)/k(i1))^2;
	end
	end

end
end
totalElastic=totalElastic*pi;
totalInelastic=totalInelastic*pi;

end
