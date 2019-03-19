function [U,E]=qqcoupled(lambda)
%function [U,E]=qqcoupled(lambda)
%
%per lambda, as provided in input, this function calculates eigenvectors and splitting due to quadrupole-quadrupole coupling
%U is given in a basis |L lambda> with L running from abs(lambda) to 4 (maximum L for D-state scandium)

Ls=abs(lambda):4;
n =length(Ls);
H= zeros(n,n);

K = 4;

for i=1:n
  Li = Ls(i);
  for j=1:n
    Lj = Ls(j);
    nj = [2 2 2 ; 2 2 2; Li Lj K];
    H(i,j) = ff_9j(nj)*sqrt((2*Li+1)*(2*Lj+1))*(-1)^(Lj-lambda)*...
             c_gm([Li Lj K; lambda -lambda 0]);
  end
end

[U,E] = seig(H);

end
