function prob5
% Problem specifics
m = 1822;
E = 1.7e-4;

V0 = 8e-4;
a  = 1;

x = linspace(-10, 10, 5e3)';
d = mean(diff(x));

% Functions which influence propagation
V   = @(x) V0*sech(a*x).^2;
k   = @(x) sqrt(2*m*(E - V(x)));

% Define incoming and outgoing asymptotic solutions
incPlane = @(A, x) A*exp(1i*k(x).*x);
outPlane = @(B, x) B*exp(-1i*k(x).*x);

R = nan(length(x),2);
R([1 2],:) = [ incPlane(1, x(1:2)) outPlane(1, x(1:2)) ];

for i = 3:length(x)
    kStep = k(x([i-2 i-1 i])).^2;
    rStep = R([i-2 i-1], :);
    R(i, :) = [ numerov(d, kStep, rStep(:,1)) numerov(d, kStep, rStep(:,2)) ];
end

% Find coefficeints
C = [incPlane(1, x([end-1 end])) outPlane(1, x([end-1 end]))] \ R([end-1 end], 1);

Sco = C(1) - abs(C(2))^2/conj(C(1))
Rco = -C(2)/conj(C(1))

S2 = abs(Sco)^2
R2 = abs(Rco)^2

% Plotting
figure;
A = 2*max(V(x))/max(max(R));
plot(x, A*R, x, V(x))

end

function out = numerov(d, k, R)
% General Numerov propagator (used for linear second order differential equations)
% INPUTS
%   d: grid spacing (grid must be equidistanct)
%   k: 3 element squared wavevector, [k_n-2, k_n-1, k_n] where k_i is the position at a specified index i, and i = n is the current index
%   R: 2 element value vector, [R_n-2, R_n-1] where R(z) is the scaled radial wavefunction being solved for

out = ((2 - (10*d^2/12)*k(2))*R(2) - (1 + (d^2/12)*k(1))*R(1))/(1 + (d^2/12)*k(3));
end