function genExPEC
% Function to calculate the excited state potential and save to output files

%% Problem specifics
% select isotope for mass and quantum defect to use
% NOTE: quantum defects are only given for homonuclear combinations in Borkowski
isotope = 84; % valid options [88, 86, 84]

% specify units for saving the output
units = 'cgs'; % valid options are ['cgs', 'au']

% define region of interest
% units are determined by that option [either bohr (au) or angstroms (cgs)]
range  = (2:0.01:300)'; 

% Base output path and filenames
basePath  = 'Z:\';
filenames = {'srVA.dat' 'srVB.dat'};

%% physical constants
aBohr  = 0.52917721067;
cm2H   = 4.556335270832134e-06;
amu2Au = 1.822889e+03;

% isotope of interest (this is also implicitly assumed in the PEC function
% as the quantum defect value (alpha) is programmed for the 86 value. See Borkowski 2014 
% for the other values.
m86 = 85.9092607309;
m84 = 83.913425;
m88 = 87.9056122571;

% centrifugal coupling term
B      = @(r, mu) (2*mu*r.^2).^(-1);
funcMu = @(m1, m2) prod([m1 m2])/sum([m1 m2]);

%% Make selections based on otpions above
% calculate reduced mass
switch isotope
    case 84
        mu = funcMu(m84, m84)*amu2Au;
    case 86
        mu = funcMu(m86, m86)*amu2Au;
    case 88
        mu = funcMu(m88, m88)*amu2Au;
end

switch units
    case 'au'
        % Default units are au
        lenConv = 1;
        enConv  = 1;
    case 'cgs'
        lenConv = aBohr;
        enConv  = 1/cm2H;     
end

% Length input for the PEC evaluation is expected in AU so the lenConv 
% uses the appropriate factor to correctly evaluate the potential
r_eval   = range*lenConv;

%% Solve for the potential
% This gets the interaction potential using Hunds case (a) potentials 
[V11, V22] = srPEC_3P1plus1S0(r_eval, isotope);
Bvec = B(r_eval, mu);

% Loop through and calculate the adiabatic potentials at each point in space
% Quick and dirty solution, this could probably be done more quickly using linear algebra
[VF, VR] = deal(zeros(length(range), 1));
for i = 1:length(range);
    % Explicitly calculating for the J = 1 case
    V = [ V11(i) + 4*Bvec(i)    -sqrt(8)*Bvec(i)
          -sqrt(8)*Bvec(i)      V22(i) + 2*Bvec(i)  ];
      
    [~, eigVal] = eig(V);
    
    % Remember that at this point everything is still in atomic units
    VF(i) = eigVal(1);
    VR(i) = eigVal(4);
end

%% Save output
% Save files and convert from atomic units to chemist units of Angstrom and eV
fnvF = fopen([basePath filenames{1}], 'w+'); 
fprintf(fnvF, '%24.15f  %24.15f\n', [range (VF*enConv)]'); 
fclose(fnvF);

fnvR = fopen([basePath filenames{2}], 'w+'); 
fprintf(fnvR, '%24.15f  %24.15f\n', [range (VR*enConv)]'); 
fclose(fnvR);

end