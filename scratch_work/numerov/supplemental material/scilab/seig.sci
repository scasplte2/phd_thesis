function [V,lambda]=seig(A)

[V,lambda]=spec(A);
[lambda,k]=gsort(diag(lambda),'r','i');
V=V(:,k);

endfunction
