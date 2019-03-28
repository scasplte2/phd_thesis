%-------------------------------------------------------------------------
%
% Use dressed states. This function basically takes a setup for naked states
% and replaces them by an expansion in dressed states. Note that the code
% has the following limitations:
%
% a) It is not possible to expand in dressed states of more than one laser
%    pulse.
% b) Only adiabatic Floquet theory is done; the wave function adiabatically
%    adapts to a change of the amplitude, and the central frequency is even
%    kept constant.
% c) As usual, the propagation is always performed in the diabatic frame of the
%    new dressed basis (usually called "dressed states"). You can produce plots
%    and calculate expectation values in the adiabatic basis (called
%    "light-induced potentials" or "Floquet states"; Cohen-Tanoudji, however,
%    calls these the "dressed states") in the usual way by setting
%    hamilt.coupling.representation to "adi", though.
%
%-------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008 Ulf Lorenz
%               2011 Ulf Lorenz
%
% see the README file for license details.

function floquet
global hamilt time

if time.efield.n_pulse>1
    util.error('Dressed states not yet implemented for more than one pulse')
end

if numel(time.efield.photons) ~= hamilt.coupling.n_eqs
    util.error('time.efield.photons needs one entry per electronic state.');
end

util.disp ('   ')
util.disp ('******************************************************')
util.disp (' Dressed (Floquet) states  ')
util.disp ('******************************************************')
util.disp ('   ')

% Save "naked" state properties
hamilt.coupling.naked = hamilt.coupling.n_eqs; % Number of naked states
naked_labels       = hamilt.coupling.labels;   % Labels of naked states
hamilt.pot.naked   = hamilt.pot.grid_ND;       % Naked potential functions
hamilt.d_x.naked   = hamilt.d_x.grid_ND;       % Naked x-dipole moments
hamilt.d_y.naked   = hamilt.d_y.grid_ND;       % Naked y-dipole moments

% Extended Hilbert space comprising dressed states
hamilt.coupling.n_eqs = numel(cell2mat(time.efield.photons));  % Number of dressed states
hamilt.coupling.labels = cell(hamilt.coupling.n_eqs,1);        % Labels of dressed states
hamilt.pot.grid_ND  = cell(hamilt.coupling.n_eqs);             % Dressed potential functions
hamilt.d_x.grid_ND  = cell(hamilt.coupling.n_eqs);             % Dressed x-dipole moments
hamilt.d_y.grid_ND  = cell(hamilt.coupling.n_eqs);             % Dressed y-dipole moments


%% Labels of "dressed" states
util.disp('Dressed states')
util.disp(' ')

% Note here and in the following: The states are ordered by the number of
% photons varying fastest, and the electronic state increasing slower.
index=0;
for m = 1:hamilt.coupling.naked
    for ph = 1:numel(time.efield.photons{m})
        index = index+1;
        num_photons = time.efield.photons{m}(ph);
        label = naked_labels{m};

        if num_photons<-1
            hamilt.coupling.labels{index} = [label ' - ' int2str(abs(num_photons)) '\omega'];
        elseif num_photons==-1
            hamilt.coupling.labels{index} = [label ' - \omega'];
        elseif num_photons==0
            hamilt.coupling.labels{index} = label;
        elseif num_photons==1
            hamilt.coupling.labels{index} = [label ' + \omega'];
        else
            hamilt.coupling.labels{index} = [label ' + ' int2str(num_photons) '\omega'];
        end

        util.disp(hamilt.coupling.labels{index})
    end
end

util.disp(' ')
util.disp(['Total number of naked   states: ' int2str(hamilt.coupling.naked)])
util.disp(['Total number of dressed states: ' int2str(hamilt.coupling.n_eqs)])
util.disp(' ')


%% Elements of the field-free Hamiltonian in Floquet space: Dressed potentials and
%  diabatic couplings.
index = 0;
for m = 1:hamilt.coupling.naked
    for ph = 1:numel(time.efield.photons{m})
        index = index + 1;
        num_photons = time.efield.photons{m}(ph);

        hamilt.pot.grid_ND{index,index} = hamilt.pot.naked{m,m} ...
                + num_photons * time.efield.frequ(1);
    end
end


%% Diabatic coupling can only take place between dressed states with the same
%number of photons.
indexleft = 0;
for m = 1:hamilt.coupling.naked
    for phm = 1:numel(time.efield.photons{m})
        indexleft = indexleft + 1;
        indexright = 0;

        photonsleft = time.efield.photons{m}(phm);

        for n = 1:hamilt.coupling.naked
            for phn = 1:numel(time.efield.photons{n})
                indexright = indexright+1;
                photonsright = time.efield.photons{n}(phn);

                % continue until we have reached the upper right triangle of the matrix
                if indexright <= indexleft
                    continue;
                end

                % Couple if there is some diabatic coupling and the same number
                % of dressing photons.
                if ~isempty(hamilt.pot.naked{m,n}) && photonsleft == photonsright
                    hamilt.pot.grid_ND{indexleft, indexright} = hamilt.pot.naked{m,n};
                end
            end
        end
    end
end


%% Dipole elements: only allowed if the states are coupled and differ by +/- one photon.
%  otherwise, this is very similar to the diabatic coupling case.
indexleft = 0;
for m = 1:hamilt.coupling.naked
    for phm = 1:numel(time.efield.photons{m})
        indexleft = indexleft + 1;
        indexright = 0;

        photonsleft = time.efield.photons{m}(phm);

        for n = 1:hamilt.coupling.naked
            for phn = 1:numel(time.efield.photons{n})
                indexright = indexright + 1;
                photonsright = time.efield.photons{n}(phn);

                % as with the other matrices, we only fill states {i,j} with j >= i
                % This test automatically implies that n >= m
                if indexright < indexleft
                    continue;
                end

                % Coupling along x
                if ( ~isempty(hamilt.d_x.naked{m,n}) ) && ( abs(photonsleft-photonsright) == 1 )
                    hamilt.d_x.grid_ND{indexleft, indexright} = hamilt.d_x.naked{m, n};
                end

                % Coupling along y
                if ( ~isempty(hamilt.d_y.naked{m,n}) ) && ( abs(photonsleft-photonsright) == 1 )
                    hamilt.d_y.grid_ND{indexleft, indexright} = hamilt.d_y.naked{m, n};
                end
            end
        end
    end
end


% Calculate which transitions are allowed and print them
allow_x = zeros(hamilt.coupling.n_eqs);
allow_y = zeros(hamilt.coupling.n_eqs);

for m = 1:hamilt.coupling.n_eqs;
    for n = 1:hamilt.coupling.n_eqs;
        if ~isempty(hamilt.d_x.grid_ND{m,n}) || ~isempty(hamilt.d_x.grid_ND{n,m})
            allow_x(m,n) = 1;
        end

        if ~isempty(hamilt.d_y.grid_ND{m,n}) || ~isempty(hamilt.d_y.grid_ND{n,m})
            allow_y(m,n) = 1;
        end
    end
end

% Display coupling schemes
util.disp (' ')
if any(allow_x(:))
    util.disp ('Allowed transitions for x-polarization')
    util.disp (' ')
    util.disp (num2str(allow_x));
else
    util.disp ('No allowed transitions for x-polarization')
end

util.disp (' ')
if any(allow_y(:))
    util.disp ('Allowed transitions for y-polarization')
    util.disp (' ')
    util.disp (num2str(allow_y))
else
    util.disp ('No allowed transitions for y-polarization')
end
util.disp (' ')
