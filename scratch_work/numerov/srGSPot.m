function V = srGSPot(r)
% From Analytical potential energy curve for the X1Sigma+g state of Sr2 
% calculated with freely varied long range coefficients Ci
% Tiemann et. al (2010)

% Tiemann values are listed, but converted from their quoted values to
% atomic units using the factors below

% Conversion factors
H_cm1  = 219474.6305; % cm^{-1}, 1Ha = 27.21eV = 219474.6305cm^{-1}
eV_cm1 = 8065.54;     % cm^{-1}, 1eV = 8065.54cm-1
Ang    = 1e-10;       % m, 1 angstrom = 1e-10m
aBohr  = 0.528;       % Ang, 1a0 = 0.528 Ang

% To avoid any ambiguity, these functions are named to say what is
% converting to what
conv.cm1ToH    = (1/H_cm1);
conv.angToBohr = (1/aBohr);

% Range parameters
R_i = 3.963 * conv.angToBohr; % aBohr (7.5 Bohr)
R_a = 10.50 * conv.angToBohr; % aBohr (19.9 Bohr)

% Vectorized method for getting the potential
% 1. Break input r into pieces (using logical indexing)
% 2. Evaluate appropriate V for each section
% 3. Reassemble V(r) vector for output

r_inner   = r( r <= R_i );
r_central = r( r > R_i & r < R_a );
r_long    = r( r >= R_a );

V = [ V_inner(r_inner, conv) V_central(r_central, conv)' V_longRange(r_long, conv) ];
end

function V = V_central(r, conv)
% a values given in cm^-1, converted in vector
a1  = 6.9e-3;           
a2  = 1.5938695e4;      
a3  = -2.9646123e4;    
a4  = -6.270500e3;      
a5  = 4.4954658e4;      
a6  = 8.705240e3;       
a7  = -1.0055114e5;     
a8  = 5.94781848e5;     
a9  = -9.95237400e5;    
a10 = -1.14496527e7;    
a11 = 4.606466552e7;    
a12 = 3.74666948e7;     
a13 = -5.439156789e8;   
a14 = 9.364833437e8;    
a15 = 1.387879473e9;    
a16 = -8.4009060525e9;  
a17 = 1.5781751108e10;  
a18 = -1.5721038639e10; 
a19 = 8.376043926e9;    
a20 = -1.88984087e9;    

a_i = [a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 a12 a13 a14 a15 a16 a17 a18 a19 a20]'*conv.cm1ToH;

b   = -0.17;                        % dimensionless
R_m = 4.6719018  * conv.angToBohr;  % aBohr
T_m = -1081.6350 * conv.cm1ToH;     % Ha

V = [];
if r
    x = (r - R_m)./(r + b*R_m);
    V = T_m + bsxfun(@power, x, (1:20)')'*a_i;
end
end

function V = V_inner(r, conv)
n   = 9.639536875;                                    % dimensionless
A   = -1.7126341e3 * conv.cm1ToH;                     % Ha
B   = 1.00304862e9 * (conv.cm1ToH*conv.angToBohr^n);  % Ha * (aBohr)^n

V = [];
if r; 
    V = A + B./r.^n; 
end
end

function V = V_longRange(r, conv)
C6  = 1.52701e7 * (conv.cm1ToH*conv.angToBohr^6);  % Ha * (aBohr)^6
C8  = 5.0654e8  * (conv.cm1ToH*conv.angToBohr^8);  % Ha * (aBohr)^8
C10 = 1.9752e10 * (conv.cm1ToH*conv.angToBohr^10); % Ha * (aBohr)^10

V = [];
if r
    V = -C6./r.^6 - C8./r.^8 - C10./r.^10;
end
end