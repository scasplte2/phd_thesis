%--------------------------------------------------------------------------
%
% Evaluate potential energy and forces which
% is needed for all types of fully classical
% or quantum-classical trajectory propagations.
% 
% Using either diabatic or adiabatic representation.
% In the latter case, all the diabatic potential and
% force matrices are not only used for calculating 
% adiabatic potentials and forces but they are also
% saved for use in tdse_ham.m (SH and FSSH only)
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt and Leonardo Cancissu-Araujo
%
% see the README file for license details.

function eval_V_F (obj)
global hamilt space 

switch hamilt.coupling.n_eqs
    case 1 % Single channel
        if ~isa (hamilt.pot,'pot.generic')
            obj.pot = V ( hamilt.pot{1,1}, obj.pos );
            obj.frc = F ( hamilt.pot{1,1}, obj.pos );
        else
            obj.pot = zeros ( size ( obj.pos{1} ) );
            for d=1:space.n_dim
                obj.frc{d} = zeros ( size ( obj.pos{1} ) );
            end
        end
        
    case 2 % Two channels
        
        % Get potentials/forces in diabatic representation
        if strcmpi(hamilt.coupling.represent,'dia') 
            
            % Get position vectors of trajectories in first channels
            r_1 = cell(space.n_dim,1);
            for d=1:space.n_dim
                r_1{d} = obj.pos{d}(obj.cha==1);
            end
            
            % Get corresponding potentials and forces
            if ~isempty(r_1)
                V_1 = V ( hamilt.pot{1,1}, r_1 );
                F_1 = F ( hamilt.pot{1,1}, r_1 );
                obj.pot (obj.cha==1) = V_1;
                for d=1:space.n_dim
                    obj.frc{d}(obj.cha==1) = F_1 {d};
                end         
            end
 
            % Get position vectors of trajectories in second channels
            r_2 = cell(space.n_dim,1);
            for d=1:space.n_dim
                r_2{d} = obj.pos{d}(obj.cha==2);
            end
            
            % Get corresponding potentials and forces
            if ~isempty(r_2)
                V_2 = V ( hamilt.pot{2,2}, r_2 );
                F_2 = F ( hamilt.pot{2,2}, r_2 );
                obj.pot (obj.cha==2) = V_2;
                for d=1:space.n_dim
                    obj.frc{d}(obj.cha==2) = F_2 {d};
                end
            end

        % Get potentials/forces in adiabatic representation 
        else
            
            % Get position vectors of trajectories in all(!) channels
            r = cell(space.n_dim,1);
            for d=1:space.n_dim
                r{d} = obj.pos{d}(:);
            end
            
            % Get diabatic potential and force matrices of all trajectories
            V_11 = V ( hamilt.pot{1,1}, r );
            F_11 = F ( hamilt.pot{1,1}, r );
            V_22 = V ( hamilt.pot{2,2}, r );
            F_22 = F ( hamilt.pot{2,2}, r );
            V_12 = V ( hamilt.pot{1,2}, r );
            F_12 = F ( hamilt.pot{1,2}, r );
            
            % Save diabatic potential & force matrices for use in tdse_ham
            obj.pot_mat(1,1,:) = V_11;
            obj.pot_mat(2,2,:) = V_22;
            obj.pot_mat(1,2,:) = V_12;
            obj.pot_mat(2,1,:) = V_12;
            for d=1:space.n_dim
                obj.frc_mat{d}(1,1,:) = F_11{d};
                obj.frc_mat{d}(2,2,:) = F_22{d};
                obj.frc_mat{d}(1,2,:) = F_12{d};
                obj.frc_mat{d}(2,1,:) = F_12{d};
            end
            
            % Indices of (coupled) channels
            ind_1 = obj.cha==1;
            ind_2 = obj.cha==2;
                        
            % Adiabatic potential energy curves|surfaces
            % Analytic solutions for a 2-state problem
            % See e.g. Eq. (4.8) in doi:10.1063/1.1522712
            dlt = (V_11 - V_22)/2;
            eta = (V_11 + V_22)/2;
            rho = sqrt ( dlt.^2 + V_12.^2 );
            
            obj.D_new   = zeros (2,obj.n_p);            
            obj.D_new(1,:) = eta - rho; % lower adiabat
            obj.D_new(2,:) = eta + rho; % upper adiabat

            obj.pot        = zeros(obj.n_p,1);
            obj.pot(ind_1) = obj.D_new(1,ind_1);
            obj.pot(ind_2) = obj.D_new(2,ind_2);
       
            % Forces along adiabatic potential curves|surfaces
            % Analytic solutions for a 2-state problem
            obj.frc = cell(size(F_11));
            for d=1:space.n_dim
                flt = (F_11{d} - F_22{d})/2;
                fta = (F_11{d} + F_22{d})/2;
                fra = ( dlt.* flt + V_12.*F_12{d} ) ./ rho;
                obj.frc{d}        = zeros(obj.n_p,1);
                obj.frc{d}(ind_1) = fta(ind_1) - fra(ind_1); % lower adiabat
                obj.frc{d}(ind_2) = fta(ind_2) + fra(ind_2); % lower adiabat
            end
            
            % Compute eigenvectors with explicit formulas (only for SH and FSSH)
            if ~strcmpi ( obj.q_c, 'lz_2' )

                ind_n = V_12~=0;
                n21   = (V_11 > V_22)*1.0;
                n22   = (V_11 < V_22)*1.0;
                
                b21          = V_11(ind_n)-V_22(ind_n)+2*rho(ind_n);
                b22          = 2*V_12(ind_n);
                n22(ind_n)   = b22 ./ sqrt(b21.^2 + b22.^2);
                n21(ind_n)   = b21 ./ sqrt(b21.^2 + b22.^2);

                % Save new eigenvectors
                obj.U_new = zeros (2,2,obj.n_p);
                obj.U_new(1,1,:) = - n22;
                obj.U_new(2,1,:) =   n21;
                obj.U_new(1,2,:) =   n21;
                obj.U_new(2,2,:) =   n22;
            end
            
        end
        
        
    otherwise % Multiple channels
        
        % Get potentials/forces in diabatic representation
        if strcmpi(hamilt.coupling.represent,'dia')
            
            % Loop over coupled channels
            r_m = cell(space.n_dim,1);
            for m = 1:hamilt.coupling.n_eqs

                % Get position vectors of trajectories in respective channels
                for d=1:space.n_dim
                    r_m{d} = obj.pos{d}(obj.cha==m);
                end
                
                % Get corresponding potentials and forces
                if ~isempty(r_m)
                    V_m = V ( hamilt.pot{m,m}, r_m );
                    F_m = F ( hamilt.pot{m,m}, r_m );
                    obj.pot (obj.cha==m) = V_m;
                    for d=1:space.n_dim
                        obj.frc{d}(obj.cha==m) = F_m{d};
                    end
                end
                
            end
            
        % Get potentials/forces in adiabatic representation
        else
        
            % Get position vectors of trajectories in all(!) channels
            r = cell(space.n_dim,1);
            for d=1:space.n_dim
                r{d} = obj.pos{d}(:);
            end
            
            % Calculate and save all diabatic potentials and forces in advance
            for p = 1:hamilt.coupling.n_eqs % diagonal
                for q = p:hamilt.coupling.n_eqs % off-diagonal
                    V_pq = V ( hamilt.pot{p,q}, r );
                    F_pq = F ( hamilt.pot{p,q}, r );
                    obj.pot_mat(p,q,:) = V_pq;
                    obj.pot_mat(q,p,:) = V_pq;
                    for d=1:space.n_dim
                        obj.frc_mat{d}(p,q,:) = F_pq{d};
                        obj.frc_mat{d}(q,p,:) = F_pq{d};
                    end
                end
            end
            
            % Save eigenvector matrices from previous time step (only temporary)
            U_old = obj.U_new;
            
            % Diagonalize diabatic potential matrices for each trajectory
            for t = 1:obj.n_p
                % Entries of vector D are eigenvalues
                % Columns of matrix U are right eigenvectors
                [obj.U_new(:,:,t),obj.D_new(:,t)] = eig(obj.pot_mat(:,:,t),'vector');
            end
            
            % Continuity of eigenvectors U: avoid sudden sign changes of columns
            if ~strcmpi ( obj.q_c, 'lz_2' )
                for p = 1:hamilt.coupling.n_eqs
                    minus = zeros(obj.n_p,1);
                    minus(:) = sum(obj.U_new(:,p,:) .* U_old(:,p,:)) < 0;
                    ind_minus = find(minus==1);
                    
                    if( isempty(ind_minus) == 0 )
                        obj.U_new(:,p,ind_minus) = - obj.U_new(:,p,ind_minus);
                    end
                end
            end
            
            % Loop over (coupled) channels m
            for m = 1:hamilt.coupling.n_eqs
                
                % Indices, number of trajectories in that channel
                ind_m = find(obj.cha==m);
                n_traj = nnz (obj.cha==m);
                if n_traj>0

                    % Save adiabatic potential energies
                    obj.pot(ind_m) = obj.D_new(m,ind_m);
                    
                    % Calculate and save adiabatic forces
                    % Eq. (7) from doi:10.1007/3-540-45054-8_36
                    for d=1:space.n_dim
                        
                        % Computation of the adiabatic force: F_adi_d = u_m^T * F_d * u_m
                        F_d     = obj.frc_mat{d}(:,:,ind_m);
                        u_m_F_d = reshape( sum( obj.U_new(:,m,ind_m) .* F_d , 1 ) , [hamilt.coupling.n_eqs,n_traj] );
                        u_m     = squeeze(obj.U_new(:,m,ind_m));
                        
                        obj.frc{d}(ind_m) = sum( u_m_F_d .* u_m ,1);
                        
                    end 
                end
            end
                        
        end
        
end

end