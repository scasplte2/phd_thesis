function [V,lambda]=seig(A);
% [V,D]=seig(A);
% eigenvector and eigenvalues, sorted
%

%[V,D]=eig(A);
%[lambda,k]=sort(diag(D));
%V=V(:,k);

[lambda,V]=ff_dsyev(A);
