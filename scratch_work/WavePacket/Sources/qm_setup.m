%------------------------------------------------------------------------------
%
% Clears the workspace to avoid side-effects from multiple WavePacket runs.
% Also closes/recreates log files etc. Call this function whenever you need a
% fresh workspace.
%
% https://sourceforge.net/p/wavepacket/matlab/wiki/Reference.Programs.qm_setup
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%               2008 Burkhard Schmidt
%               2011 Ulf Lorenz, Boris Schaefer-Bung
%               2012 Jeremy Rodriguez, Ulf Lorenz
%               2016 Burkhard Schmidt
%
% see the README file for license details.

function qm_setup (toggle)

% Clear variables/functions and close files (from previous runs);
clear global atomic bilinear balanced correlate control dim_red expect hamilt info plots psi space time uncert 
clear functions;
fclose('all');

% Initialize plotting of densities, expectation values, and spectra
% You do not want to have to set all these parameters by hand
init.plots();

%% Converting between SI units and atomic units
global atomic

% Constants in SI units (taken from Wikipedia, 2016)
m_e = 9.10938215E-31;                    % electron rest mass [kg]
q_e = 1.60217662E-19;                    % elementary charge [C]
h_b = 1.05457186E-34;                    % Planck's constant [Js] 
e_c = 8.9875517873682E+9;                % electric constant [(kg*m^3)/(s^2*C^2)]
n_a = 6.0221408E23;                      % Avogadro's number [1/mol]
c_v = 299792458;                         % vacuum speed of light [m/s]
k_b = 1.380648E-23;                      % Boltzmann's constant [J/K]

% Atomic unit of length
atomic.mass.KiloGram = m_e;                  % rest mass of electron [kg]
atomic.mass.Gram = m_e * 1e3;                % same, but in gram 
atomic.mass.Dalton =atomic.mass.Gram * n_a;  % same, but in g/mol 

% Atomic unit of length
atomic.length.Meter = h_b^2 / (e_c * m_e * q_e^2);     % Bohr radius [m]
atomic.length.nanoMeter = atomic.length.Meter * 1e09;
atomic.length.Angstrom =  atomic.length.Meter * 1e10;
atomic.length.picoMeter = atomic.length.Meter * 1e12;

% Atomic unit of energy (and equivalents)
atomic.energy.Joule = m_e * q_e^4 * e_c^2 / h_b^2;     % Hartree energy [J]
atomic.energy.JoulePerMole =       atomic.energy.Joule * n_a;
atomic.energy.KiloJoulePerMole =   atomic.energy.JoulePerMole / 1e3;
atomic.energy.CaloriePerMole =     atomic.energy.JoulePerMole / 4.184;
atomic.energy.KiloCaloriePerMole = atomic.energy.CaloriePerMole / 1e3;
atomic.energy.ElectronVolt =       atomic.energy.Joule / q_e;
atomic.energy.MilliElectronVolt =  atomic.energy.ElectronVolt * 1e3;
atomic.energy.MicroElectronVolt =  atomic.energy.ElectronVolt * 1e6;

% Atomic unit of temperature
atomic.temperature.Kelvin =      atomic.energy.Joule / k_b;
atomic.temperature.MilliKelvin = atomic.temperature.Kelvin * 1e3;
atomic.temperature.MicroKelvin = atomic.temperature.Kelvin * 1e6;

% Atomic unit of wavenumbers
atomic.wavenumber.PerMeter      = atomic.energy.Joule / (2*pi*h_b*c_v);         
atomic.wavenumber.PerCentiMeter = atomic.wavenumber.PerMeter / 100;         

% Atomic unit of time
atomic.time.Second = h_b / atomic.energy.Joule; 
atomic.time.NanoSecond =   atomic.time.Second * 1e09; 
atomic.time.PicoSecond =   atomic.time.Second * 1e12; 
atomic.time.FemtoSecond =  atomic.time.Second * 1e15; 
atomic.time.AttoSecond =   atomic.time.Second * 1e18; 

% Atomic unit of (angular!) frequency
atomic.frequency.Hertz =     atomic.energy.Joule / h_b; 
atomic.frequency.KiloHertz = atomic.frequency.Hertz / 1e03; 
atomic.frequency.MegaHertz = atomic.frequency.Hertz / 1e06; 
atomic.frequency.GigaHertz = atomic.frequency.Hertz / 1e09; 
atomic.frequency.TeraHertz = atomic.frequency.Hertz / 1e12; 
atomic.frequency.PetaHertz = atomic.frequency.Hertz / 1e15; 

% Atomic unit of dipole moment
atomic.dipole.CoulombMeter     = atomic.length.Meter * q_e; 
atomic.dipole.Debye            = atomic.dipole.CoulombMeter * c_v * 1e21; 

% Atomic unit of electric field
atomic.efield.VoltPerMeter     = atomic.energy.Joule / (q_e * atomic.length.Meter); 
atomic.efield.MegaVoltPerMeter = atomic.efield.VoltPerMeter / 1e06; 
atomic.efield.GigaVoltPerMeter = atomic.efield.VoltPerMeter / 1e09; 
atomic.efield.VoltPerAngstrom  = atomic.efield.VoltPerMeter / 1e10; 

% Atomic unit of light intensity
atomic.intensity.WattPerMeter2          = atomic.efield.VoltPerMeter^2*c_v/(8*pi*e_c); 
atomic.intensity.WattPerCentiMeter2     = atomic.intensity.WattPerMeter2 / 1e04;
atomic.intensity.MegaWattPerCentiMeter2 = atomic.intensity.WattPerCentiMeter2 / 1e06;
atomic.intensity.GigaWattPerCentiMeter2 = atomic.intensity.WattPerCentiMeter2 / 1e09;
atomic.intensity.TeraWattPerCentiMeter2 = atomic.intensity.WattPerCentiMeter2 / 1e12;
atomic.intensity.PetaWattPerCentiMeter2 = atomic.intensity.WattPerCentiMeter2 / 1e15;

% Output
if nargin==0
    toggle = false;
end
if toggle
    format longE
    util.disp ('   ')
    util.disp ('Using atomic units throughout WavePacket!')
    util.disp ('   ')
    util.disp (atomic.mass)
    util.disp (atomic.length)
    util.disp (atomic.energy)
    util.disp (atomic.temperature)
    util.disp (atomic.wavenumber)
    util.disp (atomic.time)
    util.disp (atomic.frequency)
    util.disp (atomic.dipole)
    util.disp (atomic.efield)
    util.disp (atomic.intensity)
    format short
end
