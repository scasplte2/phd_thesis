function [V11, V22] = srPEC_3P1plus1S0(r, isotope)
% All values given in atomic units from Borkowski et al 2014 except where noted
%
%%%%%%  ^3\Sigma_u^+  %% ^3\Pi_u  %%%%%%
A0  = [  1.29406314e2   5.78723038e6  ];
A1  = [ -7.90551852e1  -3.46113235e6  ];
A2  = [  1.87863441e1   7.79019763e5  ];
A3  = [ -1.96979418    -7.85317879e4  ];
A4  = [  7.88636443e-2  3.01833743e3  ];

gam   = [ 7.61382806e-2  1.34967817e-3 ];
beta  = [ 1              1.03238202    ];

switch isotope
    case 84
        alpha = [ 0.045485       1.9893     ]; % 84 value (from Reschovsky thesis 2017)
        %alpha = [ 0.045301189    1.99037413 ];
    case 86
        alpha = [ 0.045690735    1.99188286 ];
    case 88
        alpha = [ 0.045648282    1.99168225 ];
    otherwise
        error('Invalid isotope choice. Please choose either 84, 86, or 88')
end

C12 = [ -5.31841848e9 -1.06415514e10 ];
C10 = [ 2.20495e8      5.24064e7     ];
C8  = [ 2.3574797e6    3.4156471e5   ];
C6  = [ 4.3015063e3    3.8683912e3   ];
C3  = 1.52356615e-2;

% bsxfun requires a specific orientation
if ~iscolumn(r); r = r'; end

dampSum  = @(r, beta, n) bsxfun(@times, bsxfun(@power, beta.*r', (0:n)'), 1./factorial(0:n)');
funcDamp = @(r, beta, n) 1 - exp(-beta.*r).*sum( dampSum(r, beta, n) )';

funcV = @(r, gam, beta, alpha, A, C) exp(-alpha*r - gam*r.^2).*( A(1) + A(2)*r + A(3)*r.^2 + A(4)*r.^3 + A(5)*r.^4 ) ...
                                     - C(1)*funcDamp(r, beta, 12).*r.^(-12) - C(2)*funcDamp(r, beta, 10).*r.^(-10) ...
                                     - C(3)*funcDamp(r, beta, 8).*r.^(-8)   - C(4)*funcDamp(r, beta, 6).*r.^(-6)   ;
      
V_3sig = @(r) funcV(r, gam(1), beta(1), alpha(1), [A0(1) A1(1) A2(1) A3(1) A4(1)], [C12(1) C10(1) C8(1) C6(1)]);
V_3pi  = @(r) funcV(r, gam(2), beta(2), alpha(2), [A0(2) A1(2) A2(2) A3(2) A4(2)], [C12(2) C10(2) C8(2) C6(2)]);

V11 = V_3pi(r) - C3./r.^3;
V22 = (1/2)*( V_3pi(r) + V_3sig(r) ) + C3./(2*r.^3);

end