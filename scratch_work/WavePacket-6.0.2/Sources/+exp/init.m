%------------------------------------------------------------------------------
%
% Initialize expectation values and their uncertainties
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2011 Ulf Lorenz
%               2008 Burkhard Schmidt
%
% see the README file for license details.

function init
global expect space hamilt time

% Population threshold for logging, plotting
expect.min_pop = 10^-3;

% Population
expect.pop = exp.generic ('pop');

% Additional multiplicative operators
if isfield(hamilt, 'amo')
    for p = 1:length(hamilt.amo)
        if ~isempty (hamilt.amo{p})
            expect.amo{p} = exp.generic('amo');
            expect.amo{p}.ind = p;
        end
    end
end

% DVR and FBR for each spatial dimension
for k = 1:space.n_dim
    expect.pos{k} = exp.generic('pos');
    expect.pos{k}.ind = k;
    expect.mom{k} = exp.generic('mom');
    expect.mom{k}.ind = k;
end

% Potential/kinetic energy
expect.pot = exp.generic ('pot');
expect.kin = exp.generic ('kin');

% Total energy
expect.total = zeros ( time.steps.m_number, 1 );
