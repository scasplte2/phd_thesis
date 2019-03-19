function [U]=numDiagHasAtomic(mj,B)
%3 april tkarman
%function [U]=numDiagHasAtomic(mj,B) 
%diagonalizes the spin-orbit and zeeman interaction for one atom
%with projection number mj in a field of strength B.
%assumes mj is a valid projection.
%returns on the basis {|j=3/2 mj>, |j=5/2 mj>} versus {|n=3/2 mj>, |n=5/2 mj>}.

js=zeros(0);
for(j=3/2:5/2)
if(j >= abs(mj))
	js=[js; j];
end
end

gs=2.0023193043718;
A=168.34/219474.63137098*(8/20);
[n,~]=size(js);
for(i1=1:n)
for(i2=1:n)
	j1=js(i1);
	j2=js(i2);

	Has(i1,i2)=B/2*(-1)^(j1-mj+1)*sqrt((2*j1+1)*(2*j2+1))*ff_3jm([j1 1 j2; -mj 0 mj])*((-1)^(j2+1/2)*sqrt(30)*ff_6j([j1 1 j2; 2 1/2 2])+gs*(-1)^(j1+1/2)*sqrt(3/2)*ff_6j([j1 1 j2; 1/2 2 1/2]));
	if(j1==j2)
		Has(i1,i2)=Has(i1,i2)+A/2*(j1*(j1+1)-27/4);
	end
end
end

if(size(Has)==[1 1])
	U=[0 0; 0 1];
	E=Has;
else
	[U,E]=seig(Has);
end

end
