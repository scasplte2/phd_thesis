%--------------------------------------------------------------------------
%
% Evaluate Hamiltonian matrices attached
% to ALL trajectories in ALL channels
% for use in surface hopping methods.
%
% Depending on the surface hopping variant,
% these matrices are used for the following:
% - Calculating energy gaps and/or
% - Integrating the TDSEs attached
%
% In the DIABATIC case, this function provides the following:
%
% trajectories.ham (full potential matrix)
% 
% In the ADIABATIC case, this function provides the following:
%
% trajectories.ham   (diagonal only: adiabatic potential energy surfaces)
% trajectories.nac_1 (off-diagonal only: non-adiabatic coupling, 1st order)
% where the latter ones are optional (for example, not needed in LZ_2)
% 
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2018-.... Burkhard Schmidt and Leonardo Cancissu-Araujo
%
% see the README file for license details.

function eval_ham (obj, nac)
global hamilt space

switch hamilt.coupling.n_eqs
    case 1 % Single channel
        log.error ('No surface hopping for single channel')
        
    case 2 % Two channels
        
        if strcmpi(hamilt.coupling.represent,'adi') % adiabatic representation
            
            % Lower/upper adiabatic potential (ground/excited state)
            obj.ham{1,1}(:) = obj.D_new(1,:);
            obj.ham{2,2}(:) = obj.D_new(2,:);
            
            % First order non-adiabatic coupling: not needed for LZ_2
            if nac
                V_11 = squeeze(obj.pot_mat(1,1,:));
                V_22 = squeeze(obj.pot_mat(2,2,:));
                V_12 = squeeze(obj.pot_mat(1,2,:));
                
                F_11 = cell(space.n_dim,1);
                F_22 = cell(space.n_dim,1);
                F_12 = cell(space.n_dim,1);
                
                for d=1:space.n_dim
                    F_11{d} = squeeze(obj.frc_mat{d}(1,1,:));
                    F_22{d} = squeeze(obj.frc_mat{d}(2,2,:));
                    F_12{d} = squeeze(obj.frc_mat{d}(1,2,:));
                end
                
                % First order non-adiabatic coupling vectors
                % Analytic solutions for a 2-state problem
                % See e.g. Eq. (4.9) in doi:10.1063/1.1522712
                dlt = (V_11 - V_22)/2;
                rho2 = dlt.^2 + V_12.^2;
                for d=1:space.n_dim
                    flt = (F_11{d} - F_22{d})/2;
                    C_12_d = ( - dlt.*F_12{d}/2 + V_12.*flt/2 ) ./ rho2;
                    obj.nac_1{d}{1,2} = - C_12_d;
                    obj.nac_1{d}{2,1} = + C_12_d;  % anti-symmetric 
                end
            end
            
        else % diabatic representation
            V_11 = V ( hamilt.pot{1,1}, obj.pos );
            V_22 = V ( hamilt.pot{2,2}, obj.pos );
            V_12 = V ( hamilt.pot{1,2}, obj.pos );
            
            obj.ham{1,1} = V_11;
            obj.ham{2,2} = V_22;
            obj.ham{1,2} = V_12; 
            obj.ham{2,1} = V_12; % symmetric
        end
        
    otherwise % Multiple channels
        
        if strcmpi(hamilt.coupling.represent,'adi') % Adiabatic representation

            % Adiabatic potential energy surfaces
            for m = 1:hamilt.coupling.n_eqs
                obj.ham{m,m}(:) = obj.D_new(m,:);
            end
            
            % First order non-adiabatic coupling: not needed for LZ_2
            if nac
                
                % Reuse eigenvectors D, eigenvalues U from eval_V_F.m
                U = obj.U_new;
                D = obj.D_new;
                
                % Calculate non-adiabatic couplings, using also the diabatic(!) forces
                for d=1:space.n_dim
                    for m = 1:hamilt.coupling.n_eqs
                        for n = m+1:hamilt.coupling.n_eqs
                            
                            % Computation of u_m^T * F_d * u_n
                            F_d = obj.frc_mat{d}(:,:,:);
                            u_m_F_d = reshape( sum( U(:,m,:) .* F_d , 1 ) , [hamilt.coupling.n_eqs,obj.n_p] );
                            u_n = squeeze(U(:,n,:));
                            u_m_F_d_u_n = sum( u_m_F_d .* u_n ,1);
                            
                            % Computation of the NAC-vector
                            C_mn_d = u_m_F_d_u_n ./ (D(m,:) - D(n,:));
                            C_mn_d = C_mn_d'; % switch to a column vector
                            
                            obj.nac_1{d}{m,n} = - C_mn_d ;
                            obj.nac_1{d}{n,m} = + C_mn_d ;
                        end
                    end
                end
 
            end
            
        else % Diabatic formulation
            
            for m = 1:hamilt.coupling.n_eqs
                obj.ham{m,m} = V ( hamilt.pot{m,m}, obj.pos );
                for n = m+1:hamilt.coupling.n_eqs
                    V_mn = V ( hamilt.pot{m,n}, obj.pos );
                    obj.ham{m,n} = V_mn;
                    obj.ham{n,m} = V_mn;
                end
            end
        end
        
end


end

