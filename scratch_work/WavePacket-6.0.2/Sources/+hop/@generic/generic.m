%--------------------------------------------------------------------------
%
% Generic properties of all "surface hopping" class definitions
%
% For use in mixed quantum-classical dynamics where the trajectories
% may undergo transitions between different (diabatic or adiabatic)
% states of the quantum system in a stochastic manner.
%
% Note that the method traj_hop defined here is calling methods prob_hop 
% and tidy_hop which have to be implemented in each of the sub-classes.
% In addition, the sub-classes may overwrite the methods init_hop as well
% as prep_hop, thus providing the necessary flexibility in setting up 
% implementations of different variants of surface hopping.
% 
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt's group
%
% see the README file for license details.

classdef generic < handle
    
    properties (Access = public)
        
        Allowed     % Indices of trajectories where hopping is allowed
        Forbidden   % Indices of trajectories where hopping is forbidden
        
        rescale     % Toggle momentum rescaling upon surface hopping
        sca_nac     % Scaling along non-adiabatic coupling vector
        
        verbose     % Toggle diagnostic extra output

    end
    
    methods (Access = public)
        
        % Constructor
        function obj = generic
            obj.verbose = false;         % Diagnostic extra output

        end

        % Initialization of hopping
        function init_hop (obj)
            global hamilt
            
            % Indices of trajectories where hopping is allowed/forbidden
            obj.Allowed   = cell(hamilt.coupling.n_eqs);
            obj.Forbidden = cell(hamilt.coupling.n_eqs);
                        
            % Toggle momentum rescaling upon surface hopping
            if isempty(obj.rescale)
                if strcmpi(hamilt.coupling.represent,'dia')
                    obj.rescale = false;
                elseif strcmpi(hamilt.coupling.represent,'adi')
                    obj.rescale  = true;
                end
            end
            
            % Scaling along non-adiabatic coupling vectors
            % Else, scaling is along the momentum vectors
            if obj.rescale
                if isempty(obj.sca_nac)
                    obj.sca_nac = false;
                end
            end
            
            % Scaling along NAC vectors only in adiabatic representation
            if strcmpi(hamilt.coupling.represent,'dia') && obj.sca_nac
               log.error ('Scaling along NAC vectors only in adiabatic representation') 
            end
            
        end
        
        % Display surface hopping details, overloading default disp method
        function disp(obj)
            
            % Logfile output
            log.disp ('***************************************************************')
            log.disp ('Surface hopping trajectories')
            log.disp ('***************************************************************')
            log.disp (' ')
            log.disp (['Rescaling momenta upon surface hopping : ' int2str(obj.rescale)])
            if obj.sca_nac
                log.disp ('Rescaling along non-adiabatic coupling vectors')
            else
                log.disp ('Rescaling along momentum vectors')
            end
            log.disp (' ')

        end
        
        % Preprocessing: before hopping
        function prep_hop ( obj,state,first_call )
            global hamilt
            
            % Reset statistics for allowed/forbidden hopping events
            % So far, this is used in scatter plots only
            if first_call
                obj.Allowed   = cell (hamilt.coupling.n_eqs);
                obj.Forbidden = cell (hamilt.coupling.n_eqs);
            end
            
        end
        
        % Perform the trajectory hopping
        function traj_hop (obj, state)
            global hamilt space
            
            % Indices of trajectories about to hop in current time step
            ind_h = cell (hamilt.coupling.n_eqs);
            
            % Current quantum state "m"
            for m = 1:hamilt.coupling.n_eqs
                
                % Find indices of trajectories initially in "active" channel "m"
                ind_m = find (state.cha==m);
                if ~isempty (ind_m)
                    
                    % Uniform random numbers in the interval (0,1)
                    zeta = rand(size(ind_m));
                    
                    % Initialize summations of hopping probabilities
                    summ = zeros(size(ind_m));
                    
                    % Loop over all other quantum states "n"
                    for n = 1:hamilt.coupling.n_eqs
                        
                        % Adiabatic representation: neighboring states only
                        % Diabatic representation: All other states
                        if n~=m && (n==m+1 || n==m-1 || strcmpi(hamilt.coupling.represent,'dia'))
                            
                            % Get probabilities of hopping (from sub-classes)
                            probable = prob_hop (obj,state,m,n,ind_m);
                            
                            % If probabilities are negative, set them to zero
                            probable ( probable<0 ) = 0;
                            
                            % Summing up probabilities, see e.g. Eq. (10) of doi:10.1063/1.5000843
                            prev = summ;
                            summ = summ + probable;
                            
                            % Find indices of hopping trajectories by comparison with
                            % zeta, a uniform random number in the interval (0,1)
                            ind_h{n,m} = ind_m (zeta>prev & zeta<summ );
                            
                        end
                    end
                end
            end
            
            % Double loop: Hopping from current quantum state "m" to other quantum states "n"
            for m = 1:hamilt.coupling.n_eqs
                for n = 1:hamilt.coupling.n_eqs
                    
                    % If probabilities sufficient
                    if ~isempty (ind_h{n,m})
                        
                        % Number of hopping trajectories
                        n_hop = length (ind_h{n,m});
                        if obj.verbose
                            log.disp (['Hopping from state ' int2str(m) ...
                                ' to ' int2str(n) ...
                                ': ' int2str(n_hop) ...
                                ' trajectories'])
                        end
                        
                        pos = zeros(n_hop , space.n_dim);
                        mom = zeros(n_hop , space.n_dim);
                        
                        for d = 1:space.n_dim
                            % Position/momentum while hopping
                            pos(:,d) = state.pos{d}(ind_h{n,m});
                            mom(:,d) = state.mom{d}(ind_h{n,m});
                        end
                        
                        % Potential energy before|after hopping
                        pot_m = state.ham{m,m}(ind_h{n,m});
                        pot_n = state.ham{n,n}(ind_h{n,m});
                        
                        % Kinetic energy before|after hopping from energy conservation
                        kin_m = zeros(n_hop,1);
                        for d = 1:space.n_dim
                            kin_m = kin_m + ( mom(:,d).^2 ./ (2*space.dof{d}.mass));
                        end
                        kin_n = kin_m - (pot_n - pot_m);
                        
                        % Verbose
                        if obj.verbose
                            for j=1:n_hop
                                if kin_n(j)<0
                                    doit = 'FORBIDDEN ';
                                else
                                    doit = 'allowed   ';
                                end
                                log.disp ( [ doit...
                                    int2str(ind_h (j)) ': ' ...
                                    num2str(pos (j,:)) ' | ' ...
                                    num2str(pot_m (j)) ' | ' ...
                                    num2str(pot_n (j)) ' # ' ...
                                    num2str(mom (j,:)) ' | ' ...
                                    num2str(kin_m (j)) ' | ' ...
                                    num2str(kin_n (j)) ] )
                            end
                            log.disp ('   ')
                        end
                        
                        % Energy conservation; allowed and frustated hops
                        if obj.rescale
                            
                            % Energy conservation by incrementing momenta in the direction of d_mn
                            if obj.sca_nac
                                
                                d_mn = zeros (n_hop, space.n_dim);
                                A    = zeros (n_hop, 1          );
                                B    = zeros (n_hop, 1          );
                                
                                for d = 1:space.n_dim
                                    d_mn(:,d) = state.nac_1 {d}{n,m}(ind_h{n,m});
                                    A = A + 0.5 / space.dof{d}.mass * d_mn(:,d).^2;
                                    B = B + 1.0 / space.dof{d}.mass * d_mn(:,d).*mom(:,d);
                                end
                                
                                C = pot_n - pot_m;
                                D = B ./ (2*A);
                                E = C ./ A;
                                
                                ind_allowed = (D.^2 - E >= 0);
                                allowed   = ind_h{n,m}(D.^2 - E >= 0);
                                forbidden = ind_h{n,m}(D.^2 - E < 0);
                                
                                if ~isempty(allowed)
                                    ratio = - D + (-1).^(D<0) .* sqrt( (D.^2 - E) );
                                    for d = 1:space.n_dim
                                        state.mom{d}(allowed) = state.mom{d}(allowed) + ratio(ind_allowed) .* d_mn(ind_allowed,d);
                                    end
                                end
                                
                            % Energy conservation by simple rescaling of momenta
                            else
                                
                                ind_allowed = (kin_n>=0);
                                allowed   = ind_h{n,m}(kin_n>=0);
                                forbidden = ind_h{n,m}(kin_n<0);
                                
                                if ~isempty(allowed)
                                    ratio = sqrt (kin_n./kin_m);
                                    for d = 1:space.n_dim
                                        state.mom{d}(allowed) = state.mom{d}(allowed) .* ratio(ind_allowed);
                                    end
                                end
                                
                            end
                            
                        % if rescaling not desired
                        else
                            allowed   = ind_h{n,m}(kin_n>=0);
                            forbidden = ind_h{n,m}(kin_n<0);
                        end
                        
                        % Perform hopping
                        state.cha(allowed) = n;
                        
                        % Collecting allowed/forbidden jumps during every MAIN time step
                        a1 = obj.Allowed  {n,m};
                        f1 = obj.Forbidden{n,m};
                        
                        % Unique values in array
                        obj.Allowed  {n,m} = unique ( [a1; allowed  ] );
                        obj.Forbidden{n,m} = unique ( [f1; forbidden] );
                        
                        % Tidying up after hopping (from sub-classes)
                        tidy_hop (obj,state,m,n,allowed);
                        
                    end
                end  
            end    
        end

    end
        
end
    
