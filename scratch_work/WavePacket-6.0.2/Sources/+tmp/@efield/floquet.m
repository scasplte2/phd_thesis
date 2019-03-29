%-------------------------------------------------------------------------
%
% Use dressed states. This function basically takes a setup for bare states
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
%    hamilt.coupling.represent to "adi", though.
%
%-------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2008 Ulf Lorenz
%               2011 Ulf Lorenz
%
% see the README file for license details.

function floquet (obj)
global hamilt time

if ~isfield (hamilt,'dip')
    log.error ('Dipole interaction required for Floquet!')
end

if isfield (hamilt,'pol')
    log.error ('Induced dipole interaction not yet availalable in Floquet!')
end

if obj.n_pulse>1
    log.error('Dressed states not yet implemented for more than one pulse')
end

if length(obj.photons) ~= hamilt.coupling.n_eqs
    log.error('obj.photons needs one entry per (coupled) state.');
end

log.disp ('***************************************************************')
log.disp (' Dressed (Floquet) states  ')
log.disp ('***************************************************************')
log.disp ('   ')

% Save bare state properties
bare.coupling.n_eqs  = hamilt.coupling.n_eqs;    % Number of bare states
bare.coupling.labels = hamilt.coupling.labels;   % Labels of bare states

bare.pot = cell(hamilt.coupling.n_eqs);          % Potentials for bare states
for m=1:hamilt.coupling.n_eqs
    for n = m:hamilt.coupling.n_eqs
        if ~isempty ( hamilt.pot{m,n}.dvr )
            bare.pot{m,n} = hamilt.pot{m,n};     % Objects of same class
        end
    end
end
for p = 1:length(hamilt.dip)                     % Dipole moments for bare states
    if ~isempty (hamilt.dip{p})
        bare.dip{p} = cell(hamilt.coupling.n_eqs);
        for m = 1:hamilt.coupling.n_eqs
            for n = m:hamilt.coupling.n_eqs
                if ~isempty (hamilt.dip{p}{m,n}.dvr)
                    bare.dip{p}{m,n} = hamilt.dip{p}{m,n}; % Objects of same class
                end
            end
        end
    end
end

% Note here and in the following: The states are ordered by the number of
% photons varying fastest, and the electronic state increasing slower.
index=0;
hamilt.coupling.n_eqs = 0;
for m = 1:bare.coupling.n_eqs
    hamilt.coupling.n_eqs = hamilt.coupling.n_eqs + length(obj.photons{m});
    for ph = 1:length(obj.photons{m})
        index = index+1;
        num_photons = obj.photons{m}(ph);
        label = [int2str(index) ' : ' bare.coupling.labels{m}];
        
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
        
        log.disp(hamilt.coupling.labels{index})
    end
end

log.disp(' ')
log.disp(['Total number of bare    states: ' int2str(  bare.coupling.n_eqs)])
log.disp(['Total number of dressed states: ' int2str(hamilt.coupling.n_eqs)])
log.disp(' ')


%% Elements of the field-free Hamiltonian in Floquet space: Dressed potentials
hamilt.pot = cell(hamilt.coupling.n_eqs);
index = 0;
for m = 1:bare.coupling.n_eqs
    for ph = 1:length(obj.photons{m})
        index = index + 1;
        num_photons = obj.photons{m}(ph);
        
%        hamilt.pot{index,index}      = bare.pot{m,m}; % Objects of same class: troubles with pointers!!!
        hamilt.pot{index,index}.dvr = bare.pot{m,m}.dvr ... % not an object any more!
            + num_photons * time.pulse{1}.frequ;
    end
end

%% Diabatic coupling only between dressed states with the same number of photons
indexleft = 0;
for m = 1:bare.coupling.n_eqs
    for phm = 1:length(obj.photons{m})
        indexleft = indexleft + 1;
        photonsleft = obj.photons{m}(phm);
        
        indexright = 0;
        for n = 1:bare.coupling.n_eqs
            for phn = 1:length(obj.photons{n})
                indexright = indexright+1;
                photonsright = obj.photons{n}(phn);
                
                % continue until we have reached the upper right triangle of the matrix
                if indexright <= indexleft
                    continue;
                end
                
                % Potential coupling only if the number of photons is the same
                if ~isempty(bare.pot{m,n}) && photonsleft == photonsright
%                     hamilt.pot{indexleft, indexright} = bare.pot{m,n};    % Objects of same class : troubles with pointers!!!
                    hamilt.pot{indexleft, indexright}.dvr = bare.pot{m,n}.dvr;    % not an object any more!
                else
                    hamilt.pot{indexleft, indexright}.dvr = [];
                end
            end
        end
    end
end


%% Dipole elements: only allowed if the states are coupled and differ by +/- one photon
allow = cell(length(bare.dip),1);
for p = 1:length(bare.dip)
    if ~isempty (bare.dip{p})
        hamilt.dip{p}  = cell(hamilt.coupling.n_eqs);
        allow{p} = zeros(hamilt.coupling.n_eqs);
    end
end

indexleft = 0;
for m = 1:bare.coupling.n_eqs
    for phm = 1:length(obj.photons{m})
        indexleft = indexleft + 1;
        photonsleft = obj.photons{m}(phm);
        
        indexright = 0;
        for n = 1:bare.coupling.n_eqs
            for phn = 1:length(obj.photons{n})
                indexright = indexright + 1;
                photonsright = obj.photons{n}(phn);
                
                % continue until we have reached the diagonal of the matrix
                if indexright < indexleft
                    continue;
                end
                
                % Dipole coupling only if the number of photons differ by one
                for p = 1:length(hamilt.dip)
                    if ~isempty (hamilt.dip{p})
                        if ~isempty(bare.dip{p}{m,n}) && abs(photonsleft-photonsright) == 1
%                             hamilt.dip{p}{indexleft,indexright} = bare.dip{p}{m,n}; % Objects of same class: troubles with pointers!!
                            hamilt.dip{p}{indexleft,indexright}.dvr = bare.dip{p}{m,n}.dvr; % not an object any more!
                            allow{p}(indexleft, indexright) = 1;
                            allow{p}(indexright, indexleft) = 1;
                        else
                            hamilt.dip{p}{indexleft, indexright}.dvr = [];
                        end
                    end
                end
                
                
            end
        end
    end
end


%% Display coupling schemes
for p = 1:length(bare.dip)
    if ~isempty (bare.dip{p})
        
        % Display coupling schemes
        log.disp (' ')
        if any(allow{p}(:))
            log.disp (['Allowed transitions with polarization along ' int2str(p)])
            log.disp (' ')
            log.disp (num2str(allow{p}));
        else
            log.disp (['No allowed transitions with polarization along ' int2str(p)])
        end
        
        log.disp (' ')
        
    end
end

end
