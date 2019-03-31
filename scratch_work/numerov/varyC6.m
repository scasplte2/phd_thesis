function [C6Out, aOut] = varyC6
% Currently using free vary coeffs

% Free vary
%cRange = [1.52657 1.526615];
%cRange = [1.523 1.53];

% recommended
cRange = [1.523945 1.52399];
%cRange = [1.52 1.527];

ang2Bohr = 1/0.52917721067;
cm2H     = 4.556335270832134e-06;

c = linspace(cRange(1), cRange(2), 10);
b = zeros(1,length(c));

for i = 1:length(c)
    b(i) = halo86(c(i));
end

figure;
plot(c*1e7*cm2H*(ang2Bohr)^6, b,...
    [c(1) c(end)]*1e7*cm2H*(ang2Bohr)^6, [83 83],...
    [c(1) c(end)]*1e7*cm2H*(ang2Bohr)^6, [83.2 83.2],...
    [c(1) c(end)]*1e7*cm2H*(ang2Bohr)^6, [82.8 82.8] )
axis tight

P0  = polyfit(b-83, c, 1);
C60 = P0(2)*1e7*cm2H*ang2Bohr^6;

Pp1  = polyfit(b-83.2, c, 1);
C6p1 = Pp1(2)*1e7*cm2H*ang2Bohr^6;

C6Out = [C60 C6p1-C60];

a0 = calcScatterLengths(P0(2));
a0 = a0.sr86(2);

ap1 = calcScatterLengths(Pp1(2));
ap1 = ap1.sr86(2);

aOut = [a0 a0-ap1];

end