% Initial state of a trajectory object
function init_obj (obj)

global hamilt space time

% Preallocate
obj.pos = cell(space.n_dim,1);
obj.mom = cell(space.n_dim,1); 
obj.frc = cell(space.n_dim,1);
for k = 1:space.n_dim
    obj.pos{k} = zeros (obj.n_p, 1); % positions of particles
    obj.mom{k} = zeros (obj.n_p, 1); % momenta of particles
    obj.frc{k} = zeros (obj.n_p, 1); % forces acting on particles
end
obj.pot = zeros (obj.n_p, 1); % potential energy
obj.kin = zeros (obj.n_p, 1); % kinetic energy

% Initialize positions/momenta as a direct product of 1D functions
if space.n_dim>1
    log.disp ('***************************************************************')
    log.disp ('Initialize trajectories using a direct product of 1-dim states ')
    log.disp ('***************************************************************')
    log.disp (' ')
end

if isfield(time, 'corr')
    log.error ('Fully correlated initial state for trajectories not implemented')
end

if ~isfield(time, 'dof')
    log.error ('No information on initial (product!) state found')
end

if length(time.dof) ~= space.n_dim
    log.error ('Wrong dimensionality of initial product state')
end

for k = 1:space.n_dim
    log.disp ('***************************************************************')
    log.disp (['Initial state along degree of freedom : ' int2str(k)])
    if ~isempty(time.dof{k})
        time.dof{k}.dof = k; % Tell each d-o-f about its index
        init_dof ( time.dof{k} );
        disp     ( time.dof{k} ); log.disp ( ' ' )
        traj_dof ( time.dof{k}, obj)
    else
        log.error ('No specification for this d-o-f found')
    end
end

% Population of coupled states
obj.cha = int32 ( ones (obj.n_p, 1) ); % single channel
if hamilt.coupling.n_eqs>1
    
    % Channels initially populated according to individual coefficients
    if ~isempty(hamilt.coupling.ini_coeffs)
        
        % Check number of coefficients
        if length(hamilt.coupling.ini_coeffs) ~= hamilt.coupling.n_eqs
            log.error('Wrong number of initial coefficients')
        end
        
        % Convert populations to number of trajectories
        n_t = round ( obj.n_p * hamilt.coupling.ini_coeffs.^2 );
        if sum(n_t)~=obj.n_p
            log.error ('Inconsistent number of trajectories')
        end
        
        % No transformation required if:
        % Adiabatic propagation with adiabatic initial populations OR
        %  Diabatic propagation with  diabatic initial populations
        if strcmpi(hamilt.coupling.represent,hamilt.coupling.ini_rep)
            where = 0;
            for m=1:hamilt.coupling.n_eqs
                howmany = n_t (m);
                obj.cha(where+1:where+howmany) = int32 ( m*ones (howmany, 1) );
                where = where + howmany;
            end
            
        else
            
            % Get position vectors
            r = cell(space.n_dim,1);
            for d=1:space.n_dim
                r{d} = obj.pos{d}(:);
            end
            
            % Set up potential matrices
            for p = 1:hamilt.coupling.n_eqs % diagonal
                obj.pot_mat(p,p,:) = V ( hamilt.pot{p,p}, r );
                for q = p+1:hamilt.coupling.n_eqs % off-diagonal
                    V_pq = V ( hamilt.pot{p,q}, r );
                    obj.pot_mat(p,q,:) = V_pq;
                    obj.pot_mat(q,p,:) = V_pq;
                end
            end
            
            % Loop over all trajectories
            for t = 1:obj.n_p
                
                % Diagonalize diabatic potential matrix
                [U,~] = eig(obj.pot_mat(:,:,t));
                
                % Transform from diabatic to adiabatic representation:
                % Adiabatic propagation with diabatic initial populations
                if strcmpi(hamilt.coupling.represent,'adi')
                    obj.cha(t) = randomPop ( U' * hamilt.coupling.ini_coeffs' );
                    
                    % Transform from adiabatic to diabatic representation
                    % Diabatic propagation with adiabatic initial populations
                elseif strcmpi(hamilt.coupling.represent,'dia')
                    obj.cha(t) = randomPop ( U * hamilt.coupling.ini_coeffs' );
                end
                
            end
            
        end
        
    else
        log.disp ('***************************************************************')
        log.disp ('Initial populations taken from initial function.  ')
        log.disp ('***************************************************************')
        log.disp (' ')
        log.error ('Code missing')
    end
    
    % If adiabatic representation is chosen, we'll need the (numeric)
    % eigenvalues/eigenvectors of diabatic potential energy matrices 
    if strcmpi(hamilt.coupling.represent,'adi')
        obj.U_new   = zeros (hamilt.coupling.n_eqs,hamilt.coupling.n_eqs,obj.n_p);
        obj.D_new   = zeros (hamilt.coupling.n_eqs,                      obj.n_p);
        obj.pot_mat = zeros (hamilt.coupling.n_eqs,hamilt.coupling.n_eqs,obj.n_p);
        
        obj.frc_mat = cell(space.n_dim,1);
        for d=1:space.n_dim
            obj.frc_mat{d} = zeros (hamilt.coupling.n_eqs,hamilt.coupling.n_eqs,obj.n_p);
        end
        
    end
    
    % For  quantum/classical surface hopping
    if isfield (time,'hop')
        
        % Hamiltonian matrices attached to trajectories
        obj.ham       = cell(hamilt.coupling.n_eqs);
        for m=1:hamilt.coupling.n_eqs
            for n=1:hamilt.coupling.n_eqs
                obj.ham {n,m} = zeros(obj.n_p,1);
            end
        end
        
        % First-order non-adiabatic coupling tensor
        if strcmpi(hamilt.coupling.represent,'adi')
            obj.nac_1 = cell(space.n_dim,1);
            for d = 1:space.n_dim
                obj.nac_1{d} = cell (hamilt.coupling.n_eqs);
                for m=1:hamilt.coupling.n_eqs
                    for n=m+1:hamilt.coupling.n_eqs
                        obj.nac_1{d}{n,m} = zeros(obj.n_p,1);
                        obj.nac_1{d}{m,n} = zeros(obj.n_p,1);
                    end
                end 
            end
        end
        
    end
    
end


function ind = randomPop(coeffs)
% randomly choose index of a state to be populated
% for a given vector of coefficients of the states

% Normalize coefficient vector
coeffs = coeffs / norm (coeffs);

% Single uniformly distributed random number in the interval (0,1).
xi = rand;

% Find interval 
summ = 0;
for mm=1:length(coeffs)
    prob = coeffs(mm)^2;
    if xi < summ+prob 
        break
    else
        summ = summ + prob;
    end
end

% Return index
ind = mm;

end

end

