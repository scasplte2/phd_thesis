function [Vy, alpha, lambda]=Vy()
  delta = 0.3;
  y     = (-6:delta:6)';
  V     = 0.5*y.^2;
  muy   = 1;
  [U,lambda] = sinc_dvr(muy, delta, V);
  nchan = 6;
  U     = U(:,1:nchan);
  lambda= lambda(1:nchan);
  
  A     = 41000;
  alpha = 0.3;

  vvy   = A*exp(alpha*y);
  Vy    = U'*diag(vvy,0)*U;
end

