function genExPEC
% Function to calculate the excited state potential and save to output files

%% Problem specifics
% select isotope for mass and quantum defect to use
% NOTE: quantum defects are only given for homonuclear combinations in Borkowski
isotope = 84; % valid options [88, 86, 84]

% specify units for saving the output
units = 'cgs'; % valid options are ['cgs', 'au']

% define region of interest (always specified in atomic units)
% NOTE: Length input for the PEC evaluation is expected in AU
r_roi  = (0.2:0.005:400)'; 

% Base output path and filenames
basePath  = 'Z:\';

filenames = {'srVA.dat' 'srVB.dat'}; % 84
%filenames = {'srVF.dat' 'srVR.dat'}; % 86
%filenames = {'srVY.dat' 'srVZ.dat'}; % 88

%% physical constants
ang2a0 = 1.889726125457828;
H2cm   = 219475;
amu2Au = 1.822889e3;

% isotope of interest (this is also implicitly assumed in the PEC function
% as the quantum defect value (alpha) is programmed for the 86 value. See Borkowski 2014 
% for the other values.
m84 = 83.913425*amu2Au;
m86 = 85.9092607309*amu2Au;
m88 = 87.9056122571*amu2Au;

% centrifugal coupling term
funcB  = @(r, mu) 1/(2*mu*r.^2);
funcMu = @(m1, m2) prod([m1 m2])/sum([m1 m2]);

%% Make selections based on otpions above
% calculate reduced mass
switch isotope
    case 84
        mu = funcMu(m84, m84);
    case 86
        mu = funcMu(m86, m86);
    case 88
        mu = funcMu(m88, m88);
end

switch units
    case 'au'
        % Default units are au
        lenConv = 1;
        enConv  = 1;
    case 'cgs'
        lenConv = ang2a0;
        enConv  = H2cm;
end

%% Solve for the potential
% This gets the interaction potential using Hunds case (a) potentials 
[V11, V22] = srPEC_3P1plus1S0(r_roi*lenConv, isotope);
B          = funcB(r_roi*lenConv, mu);

% Loop through and calculate the adiabatic potentials at each point in space
% Quick and dirty solution, this could probably be done more quickly using linear algebra
[V1, V2] = deal(zeros(length(r_roi), 1));
for i = 1:length(r_roi);
    % Explicitly calculating for the J = 1 case
    V = [ V11(i) + 4*B(i)    -sqrt(8)*B(i)
          -sqrt(8)*B(i)      V22(i) + 2*B(i)  ];
      
    [~, eigVal] = eig(V);
    
    % Remember that up to this point everything is still in atomic units
    V1(i) = eigVal(1,1)*enConv;
    V2(i) = eigVal(2,2)*enConv;
end

%% Save output
% Save files and convert from atomic units to chemist units of Angstrom and eV
fnv1 = fopen([basePath filenames{1}], 'w+'); 
fprintf(fnv1, '%24.15f  %24.15f\n', [r_roi V1]'); 
fclose(fnv1);

fnv2 = fopen([basePath filenames{2}], 'w+'); 
fprintf(fnv2, '%24.15f  %24.15f\n', [r_roi V2]'); 
fclose(fnv2);

end