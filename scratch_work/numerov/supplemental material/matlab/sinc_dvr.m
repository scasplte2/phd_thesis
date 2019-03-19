function [U,lambda]=sinc_dvr(mu,delta,V);
% [U,lambda]=sinc_dvr(mu,delta,V)
% 1-D Sinc-function dvr calculation
% mu:    reduced mass
% delta: grid spacing
% V:     potential
% (H - lambda[i]) U[:,i] =0
%

% build DVR matrix
n=length(V);
fac=pi*pi/(6*mu*delta*delta);
H=diag(V)+fac*eye(n);
fac=1/(mu*delta*delta);
for i=2:n
for j=1:i-1
  H(i,j)=fac*(-1)^(i-j)/(i-j)^2;
  H(j,i)=H(i,j);
end
end

% diagonalize

[U,lambda]=seig(H);

% make first element of eigenvector positive
[dummy,i1]=max(abs(U));
i2=n*(0:n-1)+i1;	% position of largest elements
ind=find(U(i2)<0);
U(:,ind)=-U(:,ind);

end
