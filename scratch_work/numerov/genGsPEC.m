% This script assumes as el=0 interaction (this only saves the electronic potential and neglects any centrifugal contribution)

% physical constants
aBohr  = 0.52917721067;
cm2H   = 4.556335270832134e-06;
amu2Au = 1.822889e+03;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For use with WavePacket - export in atomic units
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define region of interest
r_bohr = (5:0.01:3000)';

% This gets the interaction potential using Hunds case (a) potentials 
V = srPEC_1S0plus1S0(r_bohr);

% Save files and convert from atomic units to chemist units of Angstrom and eV
fnvF = fopen('C:\Users\jaa2\Desktop\phd_thesis\scratch_work\WavePacket-6.0.2\strontium\ground_state\pot_1.dat', 'w+');
fprintf(fnvF, '%24.15f  %24.15f\n', [r_bohr (V)]'); 
fclose(fnvF);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For use with pydiatomic - export angstroms and cm-1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define region of interest
x_ang  = (2:0.005:250)';
r_bohr = x_ang/aBohr;

% This gets the interaction potential using Hunds case (a) potentials 
V = srPEC_1S0plus1S0(r_bohr);

% Save files and convert from atomic units to chemist units of Angstrom and eV
fnvF = fopen('Z:\srGs.dat', 'w+'); 
fprintf(fnvF, '%24.15f  %24.15f\n', [x_ang (V/cm2H)]'); 
fclose(fnvF);