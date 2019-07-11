function genGsPEC
% Function to calculate the ground state potential and save to output files
% This script assumes as el=0 interaction (this only saves the electronic potential and neglects any centrifugal contribution)

%% Problem specifics
% specify units for saving the output
units = 'cgs'; % valid options are ['cgs', 'au']

% define region of interest (always specified in atomic units)
% NOTE: Length input for the PEC evaluation is expected in AU
r_roi  = (0.01:0.01:5300)'; 

% Base output path and filenames
basePath  = 'Z:\';
%basePath  = 'C:\Users\jaa2\Desktop\phd_thesis\scratch_work\WavePacket-6.0.2\strontium\ground_state\'

filenames = 'srGS_halo.dat';
%pot_1.dat

%% physical constants
ang2a0 = 1.889726125457828; % a_0/Ang
H2cm   = 219475; % cm/H

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
% Must query potential in atomic units
% input distance in bohr, return energy in hartree
V = srPEC_1S0plus1S0(r_roi*lenConv, 'C6', 1.526593545404851*1e7)*enConv;

%% Save output
% Save files and convert from atomic units to chemist units of Angstrom and eV
fnV = fopen([basePath filenames], 'w+'); 
fprintf(fnV, '%24.15f  %24.15f\n', [r_roi V]'); 
fclose(fnV);

end