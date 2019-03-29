%------------------------------------------------------------------------------
%
% If "state" is an object of class "wave":
% Fully quantum mechanical dynamics:
% Solves the time-dependent Schroedinger equation by 
% propagating (coupled) wave function(s) in time.
%
% If "state" is an object of class "traj":
% Mixed quantum-mechanical / classical dynamics:
% Solves the quantum-classical Liouville equation by 
% propagating (coupled) phase-space density(s) in time.
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-20.. Burkhard Schmidt
%               2007-2011 Ulf Lorenz
%
% see the README file for license details.

function qm_propa ( state )

global time;

% Initializes general information and sets up log files.
log.init (mfilename('fullpath'));

% Initialize spatial discretization
grid.init (state);

% Initialize Hamiltonian operator
init_ham (state);

% Initialize temporal discretization
init (time.steps);

% Initialize the electric field
efi.init;

% Initialize wave functions / densities
init_obj (state);

% Transform to adiabatic representation (if desired)
adiabatic ( state, -1 );

% Initialize expectation values and uncertainties
exp.init;

%% Main loop over time steps (step=1 is the initial step)
for step = 1 : time.steps.m_number
    
    % Numerical propagation using pde solvers, possibly with absorbing boundary conditions
    propagate ( state, step );
    
    % Transform to adiabatic representation (if desired)
    adiabatic ( state, step, 'dia2adi' );
    
    % Expectation values and uncertainties
    observe ( state, step );
    
    % Logging and plot title
    exp.log ( step );
        
    % Show visualization of densities and expectation values
    vis.show ( state, step );
    
    % Transform back to diabatic representation
    adiabatic ( state, step, 'adi2dia' );
        
    % Store the wave function.
    i_o.save ( state, step );
    
    % End of main loop
end

% Output clock/date/time
log.clock;

end
