%--------------------------------------------------------------------------
%
% Matrix (eigen) representations of important operators using
% output from previous qm_bound calculation
%
% This function calculates the eigenenergies, matrix elements of dipole
% and/or polarizability along x and/or y (if applicable) and optionally 
% also matrix elements of the system/bath coupling and/or additional multi-
% plicative operators. These vectors/matrices are combined in structure  
% tise and are written to file "tise.mat" in the current working directory. 
%
% An optional second parameter can be supplied to define a cut-off; matrix
% elements are set to zero if their absolute values are below the cut-off, 
% which keeps the output files more readable (and could be used for sparse
% representations later on).
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2014 - 2017 Burkhard Schmidt
%               2011 Ulf Lorenz, Boris Schaefer-Bung, Burkhard Schmidt
%               2012 Jeremy Rodriguez, Burkhard Schmidt, Ulf Lorenz
%               
%
% see the README file for license details.

function qm_matrix ( state, cutoff )

global control hamilt time

% Initializes general information and sets up log files.
log.init (mfilename('fullpath'));

% Calculations of matrix elements for wavefunctions only
if ~isa(state,'wave') 
    log.error ('Calculations of matrix elements for wavefunctions only')
end

% Provide default values for missing input arguments
if nargin<2
    cutoff=0;
end

log.disp ('***************************************************************')
log.disp ('Matrix representations of important operators     ')
log.disp ('***************************************************************')
log.disp (' ')
log.disp ('by numerical integration using DVR (quadrature) schemes     ')
log.disp ('https://sourceforge.net/p/wavepacket/matlab/wiki/Reference.Programs.qm_matrix ');
log.disp (' ')

%% Load the input
i_o.load_0 (state, true);

%% Preallocate fields (ham, dip, pol, sbc, amo) of structure 'tise' (if applicable)
tise.ham = zeros(time.steps.m_number, 1);
log.disp('Energies: \delta_ij E_i = <i|H|j> ')
log.disp(' ')

if isfield(hamilt,'dip')
    for p = 1:length(hamilt.dip)
        if ~isempty(hamilt.dip{p}{1,1}.dvr)
            tise.dip{p} = zeros(time.steps.m_number);
            log.disp(['Dipole moments: \mu_{p,ij} = <i|\mu_p|j> with p=' int2str(p)])
            log.disp(' ')
        else
            tise.dip{p} = [];
        end
    end
end

if isfield(hamilt,'pol')
    for p = 1:size(hamilt.pol,1)
        for q = p:size(hamilt.pol,2)
            if ~isempty(hamilt.pol{p,q}{1,1}.dvr)
                tise.pol{p,q} = zeros(time.steps.m_number);
                log.disp(['Polarizabilities: \alpha_{pq,ij} = <i|\alpha_pq|j> with pq='  int2str(p)  int2str(q)])
                log.disp(' ')
            else
                tise.pol{p,q} = [];
            end
        end
    end
end

if isfield(hamilt,'sbc')
    if ~isempty(hamilt.sbc{1,1}.dvr)
        tise.sbc = zeros(time.steps.m_number);
        log.disp('System-bath coupling: \chi_{ij} = <i|\chi|j> ')
        log.disp(' ')
    end
end

if isfield(hamilt, 'amo') && strcmpi(control.observe.types,'amo')
    tise.amo = cell (length(hamilt.amo),1);
    if ~isempty (hamilt.amo{p})
        for p = 1:length(hamilt.amo)
            if ~isempty(hamilt.amo{p}.dvr)
                tise.amo{p} = zeros(time.steps.m_number);
                log.disp(['Additional multiplicative operators: O_{p,ij} = <i|O_p|j> with O being: ' hamilt.amo{p}.label ])
                log.disp(' ')
            end
        end
    end
end

%% Calculate all matrix elements

% Outer loop: bra-states
for bra_index = 1:time.steps.m_number
	i_o.load_n ( state, bra_index );
	bra_state = state.dvr;
    
    % calculate the eigenenergies
    apply_ham( state, [0 0], 0 ); % provides H|psi> in psi.new
    for m = 1:hamilt.coupling.n_eqs
        tise.ham(bra_index) = wave.braket(bra_state, state.new);
    end
    
	% Inner loop: ket-states 
	for ket_index = 1:time.steps.m_number
		i_o.load_n( state, ket_index );
		ket_state = state.dvr;

        % dipole moment
        if isfield (tise, 'dip')
            for p = 1:length(hamilt.dip)
                if ~isempty(hamilt.dip{p}{1,1}.dvr)
                    tise.dip{p}(bra_index, ket_index) = wave.sandwich(bra_state, hamilt.dip{p}, ket_state);
                end
            end
        end

        % polarization
        if isfield (tise, 'pol')
            for p = 1:size(hamilt.pol,1)
                for q = p:size(hamilt.pol,2)
                    if ~isempty(hamilt.pol{p,q}{1,1}.dvr)  
                        tise.pol{p,q}(bra_index, ket_index) = wave.sandwich(bra_state, hamilt.pol{p,q}, ket_state);
                    end
                end
            end
        end
        
        % system-bath coupling
        if isfield (tise, 'sbc')
            tise.sbc(bra_index, ket_index) = wave.sandwich(bra_state, hamilt.sbc, ket_state);
        end
        
        % additional multiplicative operators
        if isfield(tise, 'amo')
            for p = 1:length(hamilt.amo)
                if ~isempty (hamilt.amo{p})
                    if ~isempty(hamilt.amo{p}.dvr)
                        tise.amo{p}(bra_index, ket_index) = wave.sandwich(bra_state, hamilt.amo{p}, ket_state);
                    end
                end
            end
        end
        
    end
end

%% Optionally truncate matrix elements
if cutoff>0
    if isfield (tise, 'dip')
        for p = 1:length(tise.dip)
            tise.dip{p}(abs(tise.dip{p}) < cutoff) = 0;
        end
    end
    
    if isfield (tise, 'pol')
        for p = 1:size(tise.pol,1)
            for q = p:size(tise.pol,2)
                tise.pol{p,q}(abs(tise.pol{p,q}) < cutoff) = 0;
            end
        end
    end
    
    if isfield (tise, 'sbc')
        tise.sbc(abs(tise.sbc) < cutoff) = 0;
    end
    
    if isfield(tise, 'amo')
        for p = 1:length(hamilt.amo)
            if ~isempty (hamilt.amo{p})
                tise.amo{p}(abs(tise.amo{p}) < cutoff) = 0;
            end
        end
    end
    
end

%% Matrix/vector representations and labels of observables (to be used in qm_abncd)
tise.lab = cell(length(control.observe.choices),1);
tise.obs = control.observe.types;

switch lower(control.observe.types)
    
    % Additional multiplicative operators
    case 'amo'
        if ~isfield (tise, 'amo')
            log.error ('No additional multiplicative operators defined')
        end
        
        % Only observables corresponding to *one* operator
        for len=1:length (control.observe.choices)
            if length(control.observe.choices{len})>1
                log.error('Combinations of more than one observables not *yet* implemented')
            end
        end

        % If not specified otherwise, observables will be labeled like operators
        if ~isfield (control.observe,'labels')
            for len=1:length (control.observe.choices)
                control.observe.labels{len} = hamilt.amo{len}.label;
            end
        end
            
        % Set labels and observables ==> structure "tise"
        tise.mat = cell(length(control.observe.choices),1);
        for len=1:length (control.observe.choices)
            tise.lab{len} = control.observe.labels{len};
            log.disp (['Observable ' int2str(len) ': Additional multiplicative operators: ' tise.lab{len}])
            tise.mat{len} = tise.amo{control.observe.choices{len}};
        end
        
    % Populations as projectors onto eigenstates
    case 'prj'
        tise.mat = cell(length(control.observe.choices),1);
        for len=1:length(control.observe.choices)
            tise.lab{len} = control.observe.labels{len};
            log.disp (['Observable ' int2str(len) ': Populations of eigenstates: ' tise.lab{len}])
            tise.mat{len} = zeros (time.steps.m_number);
            for m=1:length(control.observe.choices{len})
                ii=control.observe.choices{len}(m)+1;
                tise.mat{len}(ii,ii) = 1;
            end
        end
        
    % Populations from overlaps with eigenstates
    case 'ovl'
        tise.vec = cell(length(control.observe.choices),1);
        for len=1:length(control.observe.choices)
            tise.lab{len} = control.observe.labels{len};
            log.disp (['Observable ' int2str(len) ': Overlaps with eigenstates: ' tise.lab{len}])
            tise.vec{len} = zeros (time.steps.m_number,1);
            for m=1:length(control.observe.choices{len})
                ii=control.observe.choices{len}(m)+1;
                tise.vec{len}(ii,1) = 1;
            end
        end
        
    otherwise
        log.error('Wrong choice of observable types')
end

%% Save structure 'tise' to disk
log.disp (' ')
log.disp ('Saving all matrix representations in file: tise.mat')
save('tise','tise');

%% Plot matrices

% Open new figure
figure(11);
clf;
thisplot = vis.styles; % Construct object
show_logo (thisplot)

% Energy level diagram (diagonal matrix)
subplot(2,2,1)
plot(0:time.steps.m_number-1,real(tise.ham),'o')
h = gca;
h.LineWidth  = thisplot.l_thick;
h.FontName   = thisplot.f_name;
h.FontSize   = thisplot.f_large;
h.FontWeight = thisplot.f_heavy;
xlabel('n')
ylabel('E(n)')
title('Energy level diagram')

% System-bath coupling
if isfield (tise, 'sbc')
    subplot(2,2,2)
    surf(tise.sbc)
    h = gca;
    h.LineWidth  = thisplot.l_thick;
    h.FontName   = thisplot.f_name;
    h.FontSize   = thisplot.f_large;
    h.FontWeight = thisplot.f_heavy;    
    xlabel('n')
    ylabel('m')
    zlabel('\chi(n,m)')
    title('System-bath coupling')
end

% Dipole moments (along x and/or y)
if isfield (tise, 'dip')
    for p=1:length(tise.dip)
        if ~isempty (tise.dip{p})
            subplot(2,2,2+p)
            surf(tise.dip{p})
            h = gca;
            h.LineWidth  = thisplot.l_thick;
            h.FontName   = thisplot.f_name;
            h.FontSize   = thisplot.f_large;
            h.FontWeight = thisplot.f_heavy;            
            xlabel('n')
            ylabel('m')
            switch p
                case 1
                    zlabel('\mu_x(n,m)')
                    title('Dipole moments (x)')
                case 2
                    zlabel('\mu_y(n,m)')
                    title('Dipole moments (y)')
            end
        end
    end
end


% Output clock/date/time
log.clock;

end



