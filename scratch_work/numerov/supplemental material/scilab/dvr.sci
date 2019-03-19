function [lambda, U] = sinc_dvr(mu, delta, V, nmax)
//[lambda, U] = sinc_dvr(mu, delta, V, nmax)
//  1-D Sinc-function dvr calculation
//  mu:    reduced mass
//  delta: grid spacing
//  V:     potential
//  nmax:  (optional) return at most the lowest nmax energies and eigenvectors
//  (H - lambda(i)) U(:,i) = 0

//  GCG, 20-Mar-2003; 15-May-2003

  nrhs = argn(2)
  if nrhs<4, nmax=[], end
  if nmax==[], nmax=%inf, end

  //  build DVR matrix
  n    = length(V)
  fac  = %pi*%pi/(6*mu*delta*delta)
  H    = diag(V) + fac*eye(n, n)
  fac  = 1/(mu*delta*delta)
  for i=2:n
    for j=1:i-1
      H(i,j) = fac*(-1)^(i-j)/(i-j)^2
      H(j,i) = H(i,j)
    end
  end

  // diagonalize

  nmax = min(n, nmax)
  ind  = 1:nmax;

  if argn(1)<=1
    lambda = seig(H)
    lambda = lambda(ind)
  else
    [U,lambda] = seig(H)
    lambda     = lambda(ind)
    U          = U(:, ind)
    // make first element of eigenvector positive
    [dummy,i1] = max(abs(U),'r')
    i2         = n*(0:nmax-1)+i1         // position of largest elements
    ind        = find(U(i2)<0)
    U(:,ind)   = -U(:,ind)
  end
endfunction

