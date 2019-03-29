%--------------------------------------------------------------------------
%
% Creates movie output from previous qm_bound or qm_propa calculations
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2010-2011 Ulf Lorenz
%
% see the README file for license details.

function qm_movie(state)

global time;

% Initializes general information and sets up log files.
log.init (mfilename('fullpath'));

log.disp ('***************************************************************');
log.disp ('Make movie from saved calculation in file        ');
log.disp ('https://sourceforge.net/p/wavepacket/matlab/wiki/Reference.Programs.qm_movie')
log.disp ('***************************************************************');
log.disp (' ');

% Load the general data of the saved calculation.
% Sets most global variables, except for: plots
% which should be set before running this function
i_o.load_0 (state, true);

%% Main loop over time steps (step=1 is the initial step)
for step  = 1:time.steps.m_number
    
    % Load states from data files
    i_o.load_n ( state, step );
    
    % Transform to adiabatic representation (if desired)
    % adiabatic ( state, 'dia2adi' );
    
    % Logging and plot title
    exp.log ( step );
    
    % Show visualization of densities and expectation values
    vis.show ( state, step );
    
    % Transform back to diabatic representation
    % adiabatic ( state, 'adi2dia' );
    
end

% Output clock/date/time
log.clock;

end
