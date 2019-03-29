%--------------------------------------------------------------------------
%
% Initialize (classical) Hamiltonian function
% for use with densities sampled by trajectories
%
% Notice that for the while being this function is still limited
%
% - no truncation of energies
% - no other than Cartesian coordinates
% - no "extra" kinetic operators
% - no negative imaginary potential
% - no coupling to external fields
% - no system-bath coupling
%
% but these features will be (hopefully) implemented in future versions
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-20.. Burkhard Schmidt's group
%               2007-2008 Ulf Lorenz
%               2011 Boris Schaefer-Bung, Ulf Lorenz
%
% see the README file for license details.

function init_ham (state)
global hamilt space time

%% Close coupling scheme
init (hamilt.coupling, state);
disp (hamilt.coupling);

%% Optionally truncate all energies below/above thresholds
if space.n_dim<3
    
    if ~isfield (hamilt,'truncate') % If not already existing
        hamilt.truncate = ham.truncate; % Construct object
    end
    init (hamilt.truncate);
    disp (hamilt.truncate);
    
end

%% Grid representation of kinetic energy (diabatic: scalar)
if space.n_dim<3
    
    % Required for some reason inside grid_.../init_kin
    time.steps.s_delta = 1e-10; % set to dummy value
    
    % Kinetic operators associated with each type of DVR
    % So far, the classical codes are restricted to Cartesian coordinates
    for k = 1:space.n_dim
        if isa (space.dof{k},'grid.fft')
            init_kin (space.dof{k}, 1);
        else
            log.error('No code for other than Cartesian coordinates available')
        end
    end
    
    % "Extra" kinetic operators
    if isfield(hamilt, 'kin')
        log.error('No code for "extra" kinetic operators available')
    end
    
end

%% Set up potential energy (diabatic: matrix representation)
if ~isfield (hamilt, 'pot')
    hamilt.pot = cell(hamilt.coupling.n_eqs);
end
for m = 1:hamilt.coupling.n_eqs
    for n = m:hamilt.coupling.n_eqs
        log.disp('***************************************************************')
        log.disp(['Potential energy for channels (' hamilt.coupling.labels{m} ', ' hamilt.coupling.labels{n} '):'])
        if isempty (hamilt.pot{m,n})
            hamilt.pot{m,n} = pot.generic;
        end
        hamilt.pot{m,n}.row = m; % Tell each POT its channels
        hamilt.pot{m,n}.col = n; % Tell each POT its channels
        init_pot( hamilt.pot{m,n} );
        disp    ( hamilt.pot{m,n} ); log.disp ( ' ' )
        if space.n_dim<3
            grid_pot( hamilt.pot{m,n} );
        end
    end
end

%% Set up additional multiplicative operators
if isfield(hamilt, 'amo')
    for p = 1:length(hamilt.amo)
        if ~isempty (hamilt.amo{p})
            log.disp('***************************************************************')
            log.disp([ int2str(p) '-th additional multiplicative operators for channels (' hamilt.coupling.labels{m} ', ' hamilt.coupling.labels{n} '):'])
            if isempty (hamilt.amo{p})
                hamilt.amo{p} = amo.generic;
            end
            hamilt.amo{p}.ind = p; % Tell each AMO its index
            init_amo( hamilt.amo{p} );
            disp    ( hamilt.amo{p} ); log.disp ( ' ' )
        end
    end
else
    log.disp('***************************************************************')
    log.disp('No additional multiplicative operators available  ')
    log.disp('***************************************************************')
    log.disp('   ')
end

%% Set up system-bath coupling (diabatic: matrix representation)
if isfield(hamilt, 'sbc')
    log.warning('No TRAJ code for system-bath coupling available')
end

%% Set up negative imaginary potential (diabatic: vector)
if isfield(hamilt, 'nip')
    log.warning('No TRAJ code for negative imaginary potentials available')
end

%% Set up dipole moments (diabatic: matrix representation)
if isfield(hamilt, 'dip')
    log.warning('No TRAJ code for dipole moments available')
end

%% Set up polariabilities (diabatic: matrix representation)
if isfield(hamilt, 'pol')
    log.warning('No TRAJ code for polarizabilities available')
end

%% Truncating potential and kinetic energy
if space.n_dim<3
    trunc_pot_kin (hamilt.truncate)
end

