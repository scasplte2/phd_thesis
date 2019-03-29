%--------------------------------------------------------------------------
%
% Initialize Hamiltonian operator
% 
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2007-2008 Ulf Lorenz
%               2011 Boris Schaefer-Bung, Ulf Lorenz
%
% see the README file for license details.

function hamilt 
global hamilt space time

if ~isfield(hamilt, 'coupling')
    hamilt.coupling = [];
end

% Number of (coupled) equations
if ~isfield (hamilt.coupling, 'n_eqs')
    hamilt.coupling.n_eqs=1;
end

% Labels of (coupled) states
if ~isfield(hamilt.coupling, 'labels')
    for m = 1:hamilt.coupling.n_eqs
        hamilt.coupling.labels{m} = int2str(m);
    end
end

% Diabatic/adiabatic representation
if ~isfield (hamilt.coupling, 'representation')
    hamilt.coupling.representation='dia';
end
if hamilt.coupling.n_eqs == 1
    hamilt.coupling.representation='adi';
end

% Required for initialization: set to dummy value 
time.sub.delta = 1e-10;

%% Grid representation of kinetic energy (diabatic: scalar)
% Set up and output the kinetic energy operators, in case they are not set up later
% (ket.splitting will override this behaviour)

% Kinetic energy will be truncated at this value
if ~isfield(hamilt,'truncate')
    hamilt.truncate = [];
end
if ~isfield(hamilt.truncate,'min')
    hamilt.truncate.min = -Inf;
end
if ~isfield(hamilt.truncate,'max')
    hamilt.truncate.max = +Inf;
end
hamilt.truncate.delta = hamilt.truncate.max - hamilt.truncate.min;
for k = 1:space.size.n_dim
    space.dof{k} = init_kin(space.dof{k}, 1);
end

if isfield(hamilt, 'kin')
    for n = 1:numel(hamilt.kin)
        hamilt.kin{n} = init(hamilt.kin{n}, 1);
    end
end

%% Grid representation of potential energy (diabatic: matrix representation)
if isfield(hamilt, 'pot') && isfield(hamilt.pot, 'handle')
    hamilt.pot.grid_ND = cell(hamilt.coupling.n_eqs);
    feval ( hamilt.pot.handle );
else
    util.disp('************************************')
    util.disp('No potential energy (free particle) ')
    util.disp('************************************')
    util.disp('   ')

    for m = 1:hamilt.coupling.n_eqs
        hamilt.pot.grid_ND{m,m} = zeros(size(space.dvr.grid_ND{1}));
    end
end

%% Grid representation of negative imaginary potential (diabatic: scalar)
if isfield(hamilt, 'nip') && isfield(hamilt.nip, 'handle')
    feval ( hamilt.nip.handle );
else
    util.disp('*************************************')
    util.disp('no absorbing boundary conditions     ')
    util.disp('*************************************')
    util.disp('   ')
    hamilt.nip.grid = [];
end
    
%% Grid representation of spatial projection function (diabatic: scalar)
if isfield(space, 'prj') && isfield(space.prj, 'handle')
    feval ( space.prj.handle );
else
    util.disp('************************************')
    util.disp('No projection calculated')
    util.disp('*************************************')
    util.disp('   ')
    space.prj.grid_ND = [];
end

%% Grid representation of system-bath coupling (diabatic: matrix representation)
if isfield(hamilt, 'sbc') && isfield(hamilt.sbc, 'handle')
    hamilt.sbc.grid_ND = cell(hamilt.coupling.n_eqs);
    feval ( hamilt.sbc.handle );
else
    util.disp('************************************')
    util.disp('No system-bath coupling calculated')
    util.disp('*************************************')
    util.disp('   ')
    
end

%% Grid representation of dipole moments (diabatic: matrix representation)
if isfield(hamilt, 'dip') && isfield(hamilt.dip, 'handle') 
    % Skip preallocation because often only one polarization direction x|y
    % Needs smarter solution in future versions: see Ticket #39
	% Perhaps hamilt.dip{i} which would allow arbitrary number of components
    % hamilt.d_x.grid_ND = cell(hamilt.coupling.n_eqs);
    % hamilt.d_y.grid_ND = cell(hamilt.coupling.n_eqs);
    feval ( hamilt.dip.handle );
end

%% Initialize electric field
if ~isfield(time, 'efield')
	time.efield.n_pulse = 0;
end
if ~isfield(time.efield, 'n_pulse')
     time.efield.n_pulse = size(time.efield.shape, 1);
end

% For every polarization direction (x|y) a corresponding dipole is required
% Needs smarter solution in future versions
if time.efield.n_pulse>0
    if ~isfield(hamilt, 'd_x') || ~isfield(hamilt.d_x, 'grid_ND') || isempty(hamilt.d_x.grid_ND)
        util.error('In the presence of an electric field, you need to set up a dipole moment!');
    end

    % Initialize dressed (Floquet) states if desired
	if ~isfield(time.efield, 'dressed')
		time.efield.dressed = false;
	end

    if time.efield.dressed
        ket.floquet;
    end
end

%% Truncations

% Manually truncate potential energy
for m = 1:hamilt.coupling.n_eqs
    hamilt.pot.grid_ND{m,m}(hamilt.pot.grid_ND{m,m}>hamilt.truncate.max) = hamilt.truncate.max;
    hamilt.pot.grid_ND{m,m}(hamilt.pot.grid_ND{m,m}<hamilt.truncate.min) = hamilt.truncate.min;
end

% Estimate spectral range of Hamiltonian (neglect off-diagonal coupling)
kin_max = 0;
for k = 1:space.size.n_dim
    kin_max = kin_max + min(space.dof{k}.kin_max, hamilt.truncate.delta);
end

% Estimate spectral range of additional kinetic operators
if isfield(hamilt, 'kin')
    for k = 1:length(hamilt.kin)
        kin_max = kin_max + min(hamilt.kin{k}.kin_max, hamilt.truncate.delta);
    end
end

hamilt.pot.min = +realmax; % Largest positive floating-point number
hamilt.pot.max = -realmax;
for m=1:hamilt.coupling.n_eqs % Perhaps CELL2MAT could be used to avoid loop
    hamilt.pot.min = min ( hamilt.pot.min, min ( hamilt.pot.grid_ND{m,m}(:) ) );
    hamilt.pot.max = max ( hamilt.pot.max, max ( hamilt.pot.grid_ND{m,m}(:) ) );
end

if ~isfield(hamilt, 'min')
    hamilt.min     = hamilt.pot.min;
end
if ~isfield(hamilt, 'max')
    hamilt.max = kin_max + hamilt.pot.max;
end

hamilt.delta   = hamilt.max - hamilt.min;

%% Output

util.disp ('*******************************************************')
util.disp ('Initialize grid representation of Hamiltonian operator')
util.disp ('*******************************************************')
util.disp ( [ 'Diabatic/adiabatic representation  : ' hamilt.coupling.representation ] )
util.disp (' ')
util.disp ( [ 'Lower truncation of energies  : ' num2str(hamilt.truncate.min) ] )
util.disp ( [ 'Upper truncation of energies  : ' num2str(hamilt.truncate.max) ] )
util.disp ( ' ' )
util.disp ( [ 'Maximum of kinetic energy     : ' num2str(kin_max) ] )
util.disp ( [ 'Minimum of potential energy   : ' num2str(hamilt.pot.min) ] )
util.disp ( [ 'Maximum of potential energy   : ' num2str(hamilt.pot.max) ] )
util.disp ( ' ' )
util.disp ( [ 'Spectral range of Hamiltonian : ' num2str(hamilt.delta) ] )

if hamilt.coupling.n_eqs > 1
    util.disp ( ' ' )
    util.disp ( [ 'Number of (coupled) equations      : ' int2str(hamilt.coupling.n_eqs) ] )
    for m=1:hamilt.coupling.n_eqs
        util.disp ( [ int2str(m) ': ' hamilt.coupling.labels{m} ] )
    end
end
util.disp ( ' ' )
