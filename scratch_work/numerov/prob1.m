function [absR, absS] = prob1(E)
% Problem specifics
m = 1822;
%E = 4.1e-4;

V0 = 8e-4;
a  = 3;

x = linspace(-5, 10, 1e3)';

% Functions which influence propagation
kFree = @(x) sqrt(2*m*E);
kReal = @(x) sqrt(2*m*(E - V0));
kImag = @(x) -1i*sqrt(2*m*(V0 - E));
if E > V0; kReg2 = kReal; else kReg2 = kImag; end  

% Define incoming and outgoing asymptotic solutions
incPlane = @(A, k, x) A*exp(1i*k.*x);
outPlane = @(B, k, x) B*exp(-1i*k.*x);

dIncPlane = @(A, k, x)  1i*k*incPlane(A, k, x);
dOutPlane = @(B, k, x) -1i*k*outPlane(B, k, x);

mL = @(x) [ incPlane(1, kFree(x), x)  outPlane(1, kFree(x), x) 
            dIncPlane(1, kFree(x), x) dOutPlane(1, kFree(x), x) ];
mR = @(x) [ incPlane(1, kReg2(x), x)  outPlane(1, kReg2(x), x) 
            dIncPlane(1, kReg2(x), x) dOutPlane(1, kReg2(x), x) ];
 
M = mR(a)*(mR(0)\mL(0));
T = [ incPlane(1, kFree(a), a)  
      dIncPlane(1, kFree(a), a) ];
  
SRcoeffs = M\T;
S = 1/SRcoeffs(1);
R = S*SRcoeffs(2);

absR = abs(R).^2;
absS = abs(S).^2;

ABcoeffs = (mR(0)\mL(0))*[1; R];
A = ABcoeffs(1);
B = ABcoeffs(2);

psi1 = @(x) incPlane(1, kFree(x), x) + outPlane(R, kFree(x), x);
psi2 = @(x) incPlane(A, kReg2(x), x) + outPlane(B, kFree(x), x);
psi3 = @(x) incPlane(S, kFree(x), x);

psi = [ psi1(x( x<=0 ))' psi2(x( x>0 & x<a  ))' psi3(x( x>=a ))' ];
V   = zeros(1, length(x)); V(x>0 & x<a ) = V0;

A = 2*V0/max(psi);

plot(x, A*psi, x, V);
end