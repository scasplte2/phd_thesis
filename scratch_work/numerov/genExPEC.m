% physical constants
aBohr2Ang = 0.52917721067;
cm2H      = 4.556335270832134e-06;
amu2Au    = 1.822889e+03;

% isotope of interest (this is also implicitly assumed in the PEC function
% as the quantum defect value (alpha) is programmed for the 86 value. See Borkowski 2014 
% for the other values.
m86 = 85.9092607309*amu2Au;
m84 = 83.913425*amu2Au;
mu86 = prod([m86 m86])/sum([m86 m86]);
mu84 = prod([m84 m84])/sum([m84 m84]);

% centrifugal coupling term
B     = @(r, mu) (2*mu*r.^2).^(-1);

% define region of interest
x_ang  = (0.2:0.02:500)';
r_bohr = x_ang*aBohr2Ang;

% This gets the interaction potential using Hunds case (a) potentials 
[V11, V22] = srPEC_3P1plus1S0(r_bohr);
Bvec = B(r_bohr, mu84);

% Loop through and calculate the adiabatic potentials at each point in space
% Quick and dirty solution, this could probably be done more quickly using linear algebra
[VF, VR] = deal(zeros(length(x_ang), 1));
for i = 1:length(x_ang);
    % Explicitly calculating for the J=1 case
    V = [ V11(i) + 4*Bvec(i)    -sqrt(8)*Bvec(i)
          -sqrt(8)*Bvec(i)      V22(i) + 2*Bvec(i)  ];
      
    [~, eigVal] = eig(V);
    
    % Remember that at this point everything is still in atomic units
    VF(i) = eigVal(1);
    VR(i) = eigVal(4);
end


% Bound state energies in eV
bnd86 = [1.633 44.246 159.984 348.751]*1e6*4.13567e-15;

% Save files and convert from atomic units to chemist units of Angstrom and eV
fnvF = fopen('Z:\srVA.dat', 'w+'); 
fprintf(fnvF, '%24.15f  %24.15f\n', [x_ang (VF/cm2H)]'); 
fclose(fnvF);

fnvR = fopen('Z:\srVB.dat', 'w+'); 
fprintf(fnvR, '%24.15f  %24.15f\n', [x_ang (VR/cm2H)]'); 
fclose(fnvR);