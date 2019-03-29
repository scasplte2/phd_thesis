%--------------------------------------------------------------------------
%
% Initialize (quantum) Hamiltonian operator
% for use with wave functions represented on grids
%
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
if ~isfield (hamilt,'truncate') % If not already existing
    hamilt.truncate = ham.truncate; % Construct object
end
init (hamilt.truncate);
disp (hamilt.truncate);

%% Grid representation of kinetic energy (diabatic: scalar)

% Required for some reason inside grid_.../init_kin
time.steps.s_delta = 1e-10; % set to dummy value

% Kinetic operators associated with each type of DVR
for k = 1:space.n_dim
    init_kin (space.dof{k}, 1);
end

% "Extra" kinetic operators
if isfield(hamilt, 'kin')
    for k = 1:length(hamilt.kin)
        init_kin (hamilt.kin{k}, 1);
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
        grid_pot( hamilt.pot{m,n} );
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
            grid_amo( hamilt.amo{p} );
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
    for m = 1:hamilt.coupling.n_eqs
        for n = m:hamilt.coupling.n_eqs
            log.disp('***************************************************************')
            log.disp(['System-bath coupling for channels (' hamilt.coupling.labels{m} ', ' hamilt.coupling.labels{n} '):'])
            if isempty (hamilt.sbc{m,n})
                hamilt.sbc{m,n} = sbc.generic;
            end
            hamilt.sbc{m,n}.row = m; % Tell each SBC its channels
            hamilt.sbc{m,n}.col = n; % Tell each SBC its channels
            init_sbc( hamilt.sbc{m,n} );
            disp    ( hamilt.pot{m,n} ); log.disp ( ' ' )
            grid_sbc( hamilt.sbc{m,n} );
        end
    end
else
    log.disp('***************************************************************')
    log.disp('No system-bath coupling available                 ')
    log.disp('***************************************************************')
    log.disp(' ')
end

%% Set up negative imaginary potential (diabatic: vector)
if isfield(hamilt, 'nip')
    for m = 1:hamilt.coupling.n_eqs
        log.disp('***************************************************************')
        log.disp(['Negative imaginary potential for channel (' hamilt.coupling.labels{m} '):'])
        if isempty (hamilt.nip{m})
            hamilt.nip{m} = nip.generic;
        end
        hamilt.nip{m}.ind = m; % Tell each NIP its channel
        init_nip( hamilt.nip{m} );
        disp    ( hamilt.nip{m} ); log.disp ( ' ' )
        grid_nip( hamilt.nip{m} );
    end
else
    log.disp('***************************************************************')
    log.disp('No negative imaginary potentials available        ')
    log.disp('***************************************************************')
    log.disp(' ')
end

%% Set up dipole moments (diabatic: matrix representation)
if isfield(hamilt, 'dip')
    for p = 1:length(hamilt.dip)
        if ~isempty (hamilt.dip{p})
            for m = 1:hamilt.coupling.n_eqs
                for n = m:hamilt.coupling.n_eqs
                    log.disp('***************************************************************')
                    log.disp([ int2str(p) '-th component of dipole moment for channels (' hamilt.coupling.labels{m} ', ' hamilt.coupling.labels{n} '):'])
                    if isempty (hamilt.dip{p}{m,n})
                        hamilt.dip{p}{m,n} = dip.generic;
                    end
                    hamilt.dip{p}{m,n}.pol = p; % Tell each DIP its polarization
                    hamilt.dip{p}{m,n}.row = m; % Tell each DIP its channels
                    hamilt.dip{p}{m,n}.col = n; % Tell each DIP its channels
                    init_dip( hamilt.dip{p}{m,n} );
                    disp    ( hamilt.dip{p}{m,n} );  log.disp ( ' ' )
                    grid_dip( hamilt.dip{p}{m,n} );
                end
            end
        end
    end
else
    log.disp('***************************************************************')
    log.disp('No dipole moments available                       ')
    log.disp('***************************************************************')
    log.disp('   ')
end

%% Set up polariabilities (diabatic: matrix representation)
if isfield(hamilt, 'pol')
    for p = 1:size(hamilt.pol,1)
        for q = p:size(hamilt.pol,2)
            if ~isempty (hamilt.pol{p,q})
                for m = 1:hamilt.coupling.n_eqs
                    for n = m:hamilt.coupling.n_eqs
                        log.disp('***************************************************************')
                        log.disp([ int2str(p) '-' int2str(q) '-th component of polarizability for channels (' hamilt.coupling.labels{m} ', ' hamilt.coupling.labels{n} '):'])
                        if isempty (hamilt.pol{p,q}{m,n})
                            hamilt.pol{p,q}{m,n} = pol.generic;
                        end
                        hamilt.pol{p,q}{m,n}.p_1 = p; % Tell each POL its polarization
                        hamilt.pol{p,q}{m,n}.p_2 = q; % Tell each POL its polarization
                        hamilt.pol{p,q}{m,n}.row = m; % Tell each POL its channels
                        hamilt.pol{p,q}{m,n}.col = n; % Tell each POL its channels
                        init_pol( hamilt.pol{p,q}{m,n} );
                        disp    ( hamilt.pol{p,q}{m,n} ); log.disp ( ' ' )
                        grid_pol( hamilt.pol{p,q}{m,n} );
                    end
                end
            end
        end
    end
else
    log.disp('***************************************************************')
    log.disp('No polarizabilities available                     ')
    log.disp('***************************************************************')
    log.disp('   ')
end

%% Truncating potential and kinetic energy
trunc_pot_kin (hamilt.truncate)
