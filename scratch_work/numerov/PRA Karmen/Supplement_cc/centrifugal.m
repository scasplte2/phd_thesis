function [V]=centrifugal(qnums)
%tkarman 2 maart 2012
%function [V]=centrifugal(qnums)
%centrifugal barrier term. (matrix of \hat{l}^2, really).

[n,m]=size(qnums);
V=zeros(n,n);

for(i1=1:n)
V(i1,i1)=qnums(i1,4)*(qnums(i1,4)+1);
end


end
