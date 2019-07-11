function V = srPEC_1S0plus1S0(r, model, varargin)
% Function to calculate the analytical potential energy curve for the X1Sigma+g state of Sr2 
% calculated using the piece-wise model and long range coefficients from Tiemann et. al (2010)
%
% Tiemann values are given as states, but converted from their quoted values to
% atomic units using the factors below

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
%--------------------------------------------------------------------------
%   Name    | Input Unit  | Description
%-----------|-------------|------------------------------------------------
%   r       | bohr        | Interparticle distance
%
% Optional inputs: (for use with varargin) valid input are one of 
% {'Tiemann - Free', 'Tiemann - Rec',  'Aman - Free', 'Aman - Rec', 'varyC6'}
% If 'varyC6' is used, follow this input with the intended 
% model and the value to be used for C6 in place of the model values.
% ex. {'varyC6','free vary','1.525
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
%--------------------------------------------------------------------------
%   Name    | Input Unit  | Description
%-----------|-------------|-------------------------------------------------
%   V       | Hartree     | potential energy curve evaluated at V(r)

%% Set defaults and Evaluate variable arguements
if ~exist('model','var'); model = ''; end
printDefModelMsg = 1;

for i = 1:length(varargin);
    switch varargin{i}
        case 'noPrint'
            printDefModelMsg = 0;
        case {'free', 'rec'}
            m = i;
    end
end

%% Decide model to use
% Choose whether to use the recommended model or the freely varied model
% Aman '18 value from halo binding energy, Tiemann 2010 value from fourier spec.
switch model
    case {'Tiemann - Free', 'TF'}
        coeffsToUse = {'free vary', 'Tiemann'};
    case {'Tiemann - Rec', 'TR'}
        coeffsToUse = {'recommended', 'Tiemann'};
    case {'Aman - Free', 'AF'}
        coeffsToUse = {'free vary', 'Aman'};
    case 'varyC6'
        switch varargin{m} % uses m as determined from loop above
            case {'free'}
                coeffsToUse = {'free vary', 'other', varargin{m + 1}};
            case {'rec'}
                coeffsToUse = {'recommended', 'other', varargin{m + 1}};
            otherwise
                error('\nYou specified varyC6 but did not specify which model to use.')
        end
    otherwise %default to 'Aman - Rec'
        if printDefModelMsg;
            fprintf('\n Using defaults - Model: Recommended, C6: Aman\n')
        end
        coeffsToUse = {'recommended', 'Aman'};
end

%% Conversion factors
H_cm1  = 219474.6305;   % cm^{-1}, 1Ha = 27.21eV = 219474.6305cm^{-1}
aBohr  = 0.52917721067; % Ang, 1a0 = 0.528 Ang

% To avoid any ambiguity, these functions are named to say what is
% converting to what
conv.cm1ToH    = (1/H_cm1);
conv.angToBohr = (1/aBohr);
        
%% Setup Model
% Range parameters
R_i = 3.963 * conv.angToBohr; % aBohr (7.5 Bohr)
R_a = 10.50 * conv.angToBohr; % aBohr (19.9 Bohr)

% Vectorized method for getting the potential
% 1. Break input r into pieces (using logical indexing)
% 2. Evaluate appropriate V for each section
% 3. Reassemble V(r) vector for output

% r needs to be a column vector for the matrix algebra to work
if ~iscolumn(r); r = r'; end

r_inner   = r( r < R_i );
r_central = r( r >= R_i & r <= R_a );
r_long    = r( r > R_a );

V = [ V_inner(r_inner,     conv, coeffsToUse{1}) 
      V_central(r_central, conv, coeffsToUse{1}) 
      V_longRange(r_long,  conv, coeffsToUse) ];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Central potential
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function V = V_central(r, conv, coeffsToUse)
    % a values given in cm^-1, converted in vector
    a_i = getAVec(coeffsToUse)*conv.cm1ToH;

    % define common coefficients
    b   = -0.17;                        % dimensionless
    R_m = 4.6719018  * conv.angToBohr;  % aBohr

    switch coeffsToUse
        case 'free vary'
            T_m = -1081.6350 * conv.cm1ToH;     % Ha
        case 'recommended'
            T_m = -1081.6384 * conv.cm1ToH;     % Ha
        otherwise
            error('Invalid coefficient option')
    end

    V = [];
    if r
        x = (r - R_m)./(r + b*R_m);
        V = T_m + bsxfun(@power, x, (1:20))*a_i;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inner potential
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function V = V_inner(r, conv, coeffsToUse)
    switch coeffsToUse
        case 'free vary'
            n   = 9.639536875;                                    % dimensionless
            A   = -1.7126341e3 * conv.cm1ToH;                     % Ha
            B   = 1.00304862e9 * (conv.cm1ToH*conv.angToBohr^n);  % Ha * (aBohr)^n
        case 'recommended'
            n   =  12.362;                                           % dimensionless
            A   = -1.3328825e3    * conv.cm1ToH;                     % Ha
            B   =  3.321662099e10 * (conv.cm1ToH*conv.angToBohr^n);  % Ha * (aBohr)^n
        otherwise
            error('Invalid coefficient option')
    end

    V = [];
    if nonzeros(r); 
        V = A + B./r.^n; 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Long-range potential (ground state asymptote defines E = 0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function V = V_longRange(r, conv, coeffsToUse)
    switch coeffsToUse{1}
        case 'free vary'
            switch coeffsToUse{2}
                case 'Tiemann'
                    C6  = 1.52701 * 1e7 * (conv.cm1ToH*conv.angToBohr^6);  % Ha * (aBohr)^6
                case 'Aman'
                    C6  = 1.526593 * 1e7 * (conv.cm1ToH*conv.angToBohr^6);  % Ha * (aBohr)^6
                case 'other'
                    C6 = coeffsToUse{3} * 1e7 * (conv.cm1ToH*conv.angToBohr^6); % Ha * (aBohr)^6
            end
            
            C8  = 5.0654e8  * (conv.cm1ToH*conv.angToBohr^8);  % Ha * (aBohr)^8
            C10 = 1.9752e10 * (conv.cm1ToH*conv.angToBohr^10); % Ha * (aBohr)^10
            
        case 'recommended'
            switch coeffsToUse{2}
                case 'Tiemann'
                    C6  = 1.525 * 1e7 * (conv.cm1ToH*conv.angToBohr^6);  % Ha * (aBohr)^6
                case 'Aman'
                    C6  = 1.523967 * 1e7 * (conv.cm1ToH*conv.angToBohr^6);  % Ha * (aBohr)^6
                case 'other'
                    C6 = coeffsToUse{3} * 1e7 * (conv.cm1ToH*conv.angToBohr^6); % Ha * (aBohr)^6
            end
            
            C8  = 5.159e8 * (conv.cm1ToH*conv.angToBohr^8);  % Ha * (aBohr)^8
            C10 = 1.91e10 * (conv.cm1ToH*conv.angToBohr^10); % Ha * (aBohr)^10
            
        otherwise
            error('Invalid coefficient option')
    end

    V = [];
    if r
        V = -C6./r.^6 - C8./r.^8 - C10./r.^10;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Central part coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a_i = getAVec(coeffsToUse)
% a values given in cm^-1
    switch coeffsToUse
        case 'free vary'
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

        case 'recommended'
            a1  = -6.5e-2;
            a2  =  1.5939056e4;
            a3  = -2.9646778e4;
            a4  = -6.269777e3;
            a5  =  4.4952358e4;
            a6  =  8.709016e3;
            a7  = -1.0054929e5;
            a8  =  5.94784152e5;
            a9  = -9.95239126e5;
            a10 = -1.14496717e7;
            a11 =  4.606463055e7;
            a12 =  3.74666573e7;
            a13 = -5.439157146e8;
            a14 =  9.364833940e8;
            a15 =  1.387879748e9;
            a16 = -8.4009054730e9;
            a17 =  1.5781752106e10;
            a18 = -1.5721037673e10;
            a19 =  8.376043061e9;
            a20 = -1.88984880e9;

        otherwise
            error('Invalid coefficient option')
    end

    a_i = [a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 a12 a13 a14 a15 a16 a17 a18 a19 a20]';
end