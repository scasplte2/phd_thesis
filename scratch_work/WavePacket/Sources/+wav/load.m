%------------------------------------------------------------------------------
%
% Loads the initial wavefunction from a file generated
% in a previous run of qm_propa or qm_bound. The user has to
% specify the directory, filename, and (time step) index,
%
% It takes those arguments
% psi.init.corr.dir  - the directory where the savefile resides
% psi.init.corr.file - the file where the saved wavefunction is stored
% psi.init.corr.index- if set, gives the index of the timestep for which
%                      we want to load the wavefunction.
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2008-2010 Ulf Lorenz
%               2008 Burkhard Schmidt
%
% see the README file for license details.

function load()

global psi space hamilt

%% Start with output

util.disp('**************************************************')
util.disp('Loading initial wavefunction from file            ')
util.disp('**************************************************')
util.disp( ['Directory name : ' psi.init.corr.dir ] )
util.disp( ['File name      : ' psi.init.corr.file] )
if isfield(psi.init.corr, 'index') && ~isempty(psi.init.corr.index)
    util.disp( ['time step or index    : ' num2str(psi.init.corr.index)] )
else
    psi.init.corr.index = [];
    util.disp('Loading last index');
end
if isfield(psi.init.corr, 'state') && ~isempty(psi.init.corr.state)
    util.disp( ['state(s) to load      : ' num2str(psi.init.corr.state)])
else
    psi.init.corr.state = 1:hamilt.coupling.n_eqs;
end

% Load and check the grid data from the other calculation
setup = ket.load(psi.init.corr.dir, psi.init.corr.file);

if size(setup.space.dvr.grid_ND{1}) ~= size(space.dvr.grid_ND{1})
    util.error('Loaded initial wavefunction has wrong dimensions')
end

% Determine the index to load.
my_index = psi.init.corr.index;
if isempty(my_index)
    my_index = setup.time.main.n;
end

% Load the given index
setup = ket.load(setup, my_index);
wavefunc = setup.wf.grid_ND;

% Handle this gracefully!
if setup.hamilt.coupling.n_eqs ~= hamilt.coupling.n_eqs
    if setup.hamilt.coupling.n_eqs ~= 1
        util.error('Error: Loaded initial wavefunction has different number of electronic states.')
    else
        util.disp('WARNING: Initial wavefunction has just one electronic state!')
        util.disp('         Will place it on the first surface or according to ')
        util.disp('         initial coefficients if set.                       ')
        psi.dvr.grid_ND = wavefunc{1};
        for m = 2:hamilt.coupling.n_eqs
            psi.dvr.grid_ND{m} = zeros(size(psi.dvr.grid_ND{1}));
        end
    end
else
    for mystate = psi.init.corr.state
        psi.dvr.grid_ND{mystate} = wavefunc{mystate};
    end
end
