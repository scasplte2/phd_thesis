%--------------------------------------------------------------------------
%
% Provide input for bilinear control problems
%
% Matrix   *A*  is created from energies (and coupling to bath for LvNE)
% Vectors  *B* are created from transition dipole moments
% Matrices *N* are created from transition dipole moments
% Vectors  *C* are created from (linear) observables (LvNE)
% Matrices *D* are created from (quadratic) observables (TDSE)
% Vector *x.initial* is the initial state
% Vector *y.initial* is the corresponding output
% Vector *x.equilib* is the equilibrium state (fix point of A)
% Vector *y.equilib* is the corresponding output
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2014-17 Burkhard Schmidt
%               2012 Boris Schaefer-Bung, Burkhard Schmidt, 
%                    Ulf Lorenz, Jeremy Rodriguez
%
% see the README file for license details.
%

function qm_abncd(eom)

global control bilinear 

% Initializes general information and sets up log files.
log.init (mfilename('fullpath'));

% Load the energies, dipole, and (system-bath) coupling matrices
load ('tise');
dim=size(tise.ham, 1);

log.disp ('***************************************************************')
log.disp ('Matrices A, N, B and C, D and vectors x_i, x_e     ')
log.disp ('***************************************************************')
log.disp (' ')
log.disp ('for use in a bilinear control problem                       ')
log.disp ('https://sourceforge.net/p/wavepacket/matlab/wiki/Reference.Programs.qm_abncd/')
log.disp ('                                                            ')
log.disp ('         d                                                  ')
log.disp ('control: -- x(t) = ( A + iu(t)N ) x(t) + iu(t)B             ')
log.disp ('         dt                                                 ')
log.disp ('                          T                                 ')
log.disp ('observe: y(t) = C x(t) + x (t) D x(t)                       ')
log.disp ('                                                            ')
log.disp ('see B. Schaefer-Bung, C. Hartmann, B. Schmidt, Ch. Schuette ')
log.disp ('J. Chem. Phys. 135, 014112-1-13 (2011)                      ')

switch lower(eom)
    case 'lvne'  
        % Map matrices onto vectors using the so-called tetradic (Liouville) 
        % convention introduced by Mukamel's group
        %
        % We use columnwise ordering of the matrix elements, e.g. for
        % the density rho = (rho_00, rho_10, ..., rho_01, rho_11, ...)
        % Note that this is the standard ordering of Fortran/Matlab anyway
        
        log.disp (' ')
        log.disp (' coming from the quantum Liouville-von Neumann equation')
        log.disp (' using atomic units throughout, i.e. hbar = m_e = e = 1')
        log.disp (' for quantum systems interacting with electrical fields')
        log.disp ('                                                       ')
        log.disp ('  d             i                                      ')
        log.disp (' -- rho(t) = - ---- [H - F(t) mu, rho(t)] +L [rho(t)]  ')
        log.disp (' dt            hbar   0                    D           ')
        log.disp ('                                                       ')
        log.disp (' H_0 is a diagonal matrix for the unperturbed system,  ')
        log.disp (' mu is the Hermitian matrix with dipole moments        ')
        log.disp (' and F(t) is the electric field.                       ')
        log.disp ('                                                       ')
        log.disp (' with L [rho] = ... Lindblad dissipation/dephasing ... ')
        log.disp ('       D                                               ')
        log.disp ('                                                       ')
        log.disp (' <C>(t) = tr( C rho(t) )                               ')
        log.disp (' ')
        log.disp ('Parameters of dissipative Liouvillian: Lindblad form   ')
        log.disp ('                                                       ')
        
        % Initial and equilibrium state vectors/densities for propagation
        oct.lvne (real(tise.ham))

        % set temperature; generate title string for plots
        if ~isfield(control.lvne,'temperature')
            control.lvne.temperature = 0;
        end
        log.disp (['temperature * k_B = ' num2str(control.lvne.temperature)]) 
        log.disp (' ') 
        kBT = control.lvne.temperature;
        bilinear.title = ['LvNE: \Theta =' num2str(control.lvne.temperature) ', '];
  
        % Construct omega-matrix (Bohr frequencies)
        w=zeros(dim);
        for k=1:dim
            for l=1:dim
                w(k,l)=tise.ham(k)-tise.ham(l);
                w = real(w); % to be done: delete?
            end
        end
        
        % Set up matrix A; coherent evolution from omega-matrix
        bilinear.A = diag( reshape(-1i*w, dim^2, 1) );

        % Lindblad operators for population relaxation
        if ~isfield(control,'relax')
            control.relax = [];
        end
        if ~isfield(control.relax,'model')
            control.relax.model = 'none';
        end
        
        % Construct Gamma matrix
        % convention used throughout: Gam(k,l) = Gamma_{k <- l}
        Gam = zeros(dim);
        
        switch lower(control.relax.model)
            
            % Relaxation rates from Fermi's golden rule, see Andrianov & Saalfrank 2006
            case 'fermi'
                log.disp ('Relaxation rates from Fermi''s golden rule')
                log.disp (['relaxation rate = ' num2str(control.relax.rate)])
                log.disp (['==> upper state = ' int2str(control.relax.upper)])
                log.disp (['==> lower state = ' int2str(control.relax.lower)])
                bilinear.title = [bilinear.title '\Gamma =' num2str(control.relax.rate)];
                bilinear.title = [bilinear.title ' (' int2str(control.relax.upper)];
                bilinear.title = [bilinear.title '->' int2str(control.relax.lower) '), '] ;

                ratios = abs( triu(tise.sbc,1) ).^2 / abs( tise.sbc(control.relax.lower+1,control.relax.upper+1) )^2;
                
                if kBT == 0 % for zero temperature: only downward relaxation
                    for k=1:dim % upper right triangle of Gamma matrix
                        for l=k+1:dim % l \geq k
                            Gam(k,l)=ratios(k,l) ...
                                *w(control.relax.upper+1,control.relax.lower+1) / w(l,k) ...
                                * control.relax.rate;
                        end
                    end
                    
                else % for finite temperature
                    for k=1:dim % upper right triangle: downward relaxation
                        for l=k+1:dim % l \geq k
                            Gam(k,l) = ratios(k,l) ...
                                * w(control.relax.upper+1,control.relax.lower+1) / w(l,k)...
                                * (1-exp(-w(control.relax.upper+1,control.relax.lower+1)/kBT)) ...
                                / (1-exp(-w(l,k)/kBT)) ...
                                * control.relax.rate;
                            % upward transitions from detailed balance, see Eq. (4)
                            Gam(l,k) = Gam(k,l) * exp(-w(l,k)/kBT);
                        end
                    end
                end
                
            % Relaxation rates from Einstein's spontaneous emission
            case 'einstein'
                log.disp ('Relaxation rates from Einstein''s spontaneous emission')
                log.disp (['relaxation rate = ' num2str(control.relax.rate)])
                log.disp (['==> upper state = ' int2str(control.relax.upper)])
                log.disp (['==> lower state = ' int2str(control.relax.lower)])
                bilinear.title = [bilinear.title '\Gamma =' num2str(control.relax.rate)];
                bilinear.title = [bilinear.title ' (' int2str(control.relax.upper)];
                bilinear.title = [bilinear.title '->' int2str(control.relax.lower) '), '] ;
                
                ratios = abs( triu(tise.dip{1},1) ).^2 / abs( tise.dip{1}(control.relax.lower+1,control.relax.upper+1) )^2;
                
                for k=1:dim % upper right triangle: downward relaxation
                    for l=k+1:dim % l \geq k
                        Gam(k,l) = ratios(k,l) ...
                            * (w(l,k) / w(control.relax.upper+1,control.relax.lower+1))^3 ...
                            * control.relax.rate;
                        % upward transitions from detailed balance, see Eq. (4)
                        Gam(l,k) = Gam(k,l) * exp(-w(l,k)/kBT);
                    end
                end
                
            % Constant relaxation rates (for all downward transitions)
            case 'constant'
                log.disp ('Constant relaxation rates')
                log.disp (['relaxation rate = ' num2str(control.relax.rate)])
                bilinear.title = [bilinear.title '\Gamma =' num2str(control.relax.rate)];
                bilinear.title = [bilinear.title ' (const.), '] ;

                for k=1:dim % upper right triangle: downward relaxation
                    for l=k+1:dim % l \geq k
                        Gam(k,l)=control.relax.rate;
                        if kBT>0 % upward transitions from detailed balance
                            Gam(l,k) = Gam(k,l) * exp(-w(l,k)/kBT);
                        end
                    end
                end
                
            % Relaxation rates from data file
            case 'datafile'
                log.disp ('Reading relaxation rates from file')
                bilinear.title = [bilinear.title '\Gamma from file, '] ;

                data = load ('relax');
                Gam = data.relax;
                
            case 'none'
                log.disp ('Not accounting for relaxation rates')

            otherwise
                log.error ('Invalid choice of relaxation model')
                
        end
        
        if ~strcmpi(control.relax.model,'none')
            
            % Construct total dephasing rate, i. e. gamma matrix from Eq. (5)
            GamSummed = sum(Gam, 1);
            gamma.r = zeros(dim);
            for k=1:dim
                for l=1:dim
                    gamma.r(k,l)=0.5*(GamSummed(k)+GamSummed(l));
                end
            end
            
            % population gain, similar to Eq. (A2), but for columnwise order
            for l=0:dim-1
                ll = 1 + l*dim + l;
                for k=0:dim-1
                    kk = 1 + k*dim + k;
                    bilinear.A(ll, kk) = bilinear.A(ll, kk) + Gam(l+1,k+1);
                end
            end
            
            % Include population loss and total dephasing in matrix A
            bilinear.A = bilinear.A + diag( reshape(-gamma.r, dim^2, 1) );
            
            % Find extrema of gamma.r
            min_gr = abs(gamma.r(1,2)); min_lo=1-1; min_up=2-1;
            max_gr = abs(gamma.r(1,2)); max_lo=1-1; max_up=2-1;
            for k=1:dim % upper right triangle with diag: downward relaxation
                for l=k:dim % l \gt k
                    abs_gr = abs(gamma.r(k,l));
                    if abs_gr < min_gr
                        min_gr = abs_gr;
                        min_lo = k-1;
                        min_up = l-1;
                    end
                    if abs_gr > max_gr
                        max_gr = abs_gr;
                        max_lo = k-1;
                        max_up = l-1;
                    end
                end
            end

            log.disp (['min. relax. dephasing = ' num2str(min_gr)])
            log.disp (['==> upper state       = ' int2str(min_up)])
            log.disp (['==> lower state       = ' int2str(min_lo)])
            log.disp (['==> Bohr frequency    = ' num2str(w(min_up+1,min_lo+1))])
            log.disp (['max. relax. dephasing = ' num2str(max_gr)])
            log.disp (['==> upper state       = ' int2str(max_up)])
            log.disp (['==> lower state       = ' int2str(max_lo)])
            log.disp (['==> Bohr frequency    = ' num2str(w(max_up+1,max_lo+1))])
                 
        end
        log.disp (' ')      
        
        % Lindblad operator for pure dephasing
        if ~isfield(control,'depha')
            control.depha = [];
        end
        if ~isfield(control.depha,'model')
            control.depha.model = 'none';
        end
        
        switch lower (control.depha.model)
            
            % Dephasing rates from stochastic Gaussian model (quadratic energy gap dependence)
            case 'gauss'
                log.disp ('Dephasing rates from stochastic Gaussian model')
                log.disp (['pure dephasing rate = ' num2str(control.depha.rate)])
                log.disp (['==> upper state     = ' int2str(control.depha.upper)])
                log.disp (['==> lower state     = ' int2str(control.depha.lower)])
                log.disp (['==> Bohr frequency  = ' num2str(w(control.depha.upper+1,control.depha.lower+1))])
                bilinear.title = [bilinear.title '\Gamma^* =' num2str(control.depha.rate)];
                bilinear.title = [bilinear.title ' (' int2str(control.depha.upper)];
                bilinear.title = [bilinear.title '->' int2str(control.depha.lower) '), '] ;

                gamma.d = (w/w(control.depha.upper+1,control.depha.lower+1)).^2 * control.depha.rate;
                
                % Constant dephasing rates
            case 'constant'
                log.disp ('Constant pure dephasing rates')
                log.disp (['pure dephasing rate = ' num2str(control.depha.rate)])
                bilinear.title = [bilinear.title '\Gamma^* =' num2str(control.depha.rate)];
                bilinear.title = [bilinear.title ' (const.), '] ;

                gamma.d = ones(size(w)) * control.depha.rate;
                
                % Read dephasing rates from data file
            case 'datafile'
                log.disp ('Reading pure dephasing data from file')
                bilinear.title = [bilinear.title '\Gamma^* from file, '] ;

                data = load ('depha');
                gamma.d = data.depha;
                
            case 'none'
                log.disp ('Not accounting for pure dephasing rates')
                
            otherwise
                log.error ('Invalid choice of pure dephasing model')
                
        end
        
        if ~strcmpi(control.depha.model,'none')
            
            % Contribution to A-matrix
            bilinear.A = bilinear.A + diag( reshape(-gamma.d, dim^2, 1) );
            
            % Find extrema
            min_gd = abs(gamma.d(1,2)); min_lo=1-1; min_up=2-1;
            max_gd = abs(gamma.d(1,2)); max_lo=1-1; max_up=2-1;
            for k=1:dim-1 % upper right triangle: downward relaxation
                for l=k+1:dim % l \geq k
                    abs_gd = abs(gamma.d(k,l));
                    if abs_gd < min_gd
                        min_gd = abs_gd;
                        min_lo = k-1;
                        min_up = l-1;
                    end
                    if abs_gd > max_gd
                        max_gd = abs_gd;
                        max_lo = k-1;
                        max_up = l-1;
                    end
                end
            end
            
            log.disp (['min. pure dephasing = ' num2str(min_gd)])
            log.disp (['==> upper state     = ' int2str(min_up)])
            log.disp (['==> lower state     = ' int2str(min_lo)])
            log.disp (['==> Bohr frequency  = ' num2str(w(min_up+1,min_lo+1))])
            log.disp (['max. pure dephasing = ' num2str(max_gd)])
            log.disp (['==> upper state     = ' int2str(max_up)])
            log.disp (['==> lower state     = ' int2str(max_lo)])
            log.disp (['==> Bohr frequency  = ' num2str(w(max_up+1,max_lo+1))])
            
        end
        log.disp (' ')
        
        % Coupling(s) to control field(s)
        if isfield (tise, 'dip')
            for p = 1:length(tise.dip)
                if ~isempty (tise.dip{p})
                    bilinear.N{p}=zeros(dim^2);
                    
                    % Set up N matrices, similar to  Eq. (A3), but for columnwise order
                    for l=0:dim-1
                        for m=0:dim-1
                            index_lm = 1 + l + m*dim;
                            
                            for k=0:dim-1
                                index_km = 1 + k + m*dim;
                                index_lk = 1 + l + k*dim;
                                
                                bilinear.N{p}(index_lm,index_km) = bilinear.N{p}(index_lm,index_km) + tise.dip{p}(l+1,k+1);
                                bilinear.N{p}(index_lm,index_lk) = bilinear.N{p}(index_lm,index_lk) - tise.dip{p}(k+1,m+1);
                                
                            end
                            
                        end
                    end
                    
                    % Set up B vectors, similar to Eq. (A5), but for columnwise order
                    bilinear.B{p} = bilinear.N{p}*bilinear.x.equilib;
                    
                end
            end
        end
        
        
        
        % Choice of observables is given in control.observe.targets
        % Set up C vectors: columnwise ordering of matrices
        for len=1:length(control.observe.targets)
            bilinear.label{len} = tise.lab{control.observe.targets(len)};
            switch tise.obs
                case 'amo'
                    log.disp (['Observable ' int2str(len) ': Additional multiplicative operators: ' bilinear.label{len}])
                case 'prj'
                    log.disp (['Observable ' int2str(len) ': Populations as projectors onto eigenstates: ' bilinear.label{len}])
                otherwise
                    log.error ('Invalid choice of observable for LvNE')
            end
            % Transpose, vectorize, transpose such that C*x gives Tr(O*rho)
            op_mat = tise.mat{control.observe.targets(len)};
            op_mat = op_mat.';
            op_vec = op_mat(:);
            bilinear.C {len} = op_vec';
            bilinear.Q {len} = false; % use Re<c|x> as observable

        end
        % bilinear.C = bilinear.C'; % make C a row cell vector
        log.disp (' ')
        
        % if desired: transform full matrix problem
        % columnwise -> diagonal first (cw2df)
        if isfield (control.lvne,'order') && strcmpi(control.lvne.order,'df')
            U=math.cw2df(dim);
            
            bilinear.A = U*bilinear.A*U';
            
            for len=1:length(bilinear.C)
                bilinear.C{len} = bilinear.C{len}*U';
            end
            
            for len=1:length(bilinear.N)
                bilinear.N{len} = U*bilinear.N{len}*U';
                bilinear.B{len} = U*bilinear.B{len};
            end
            
            bilinear.x.initial = U*bilinear.x.initial;
            bilinear.x.equilib = U*bilinear.x.equilib;
        end
        
        % Values of observables for initial and equilibrium state
        for len=1:length(control.observe.targets)
            switch tise.obs
                case {'amo','prj'}
                    bilinear.y.equilib(len) = real(bilinear.C{len}*bilinear.x.equilib);
                    bilinear.y.initial(len) = real(bilinear.C{len}*bilinear.x.initial);
            end

        end
        
        % Shift initial state vector w.r.t. its equilibrium state
        bilinear.x.initial = bilinear.x.initial - bilinear.x.equilib;
        
    case 'tdse'
        
        log.disp (' ')
        log.disp (' coming from the time-dependent Schroedinger equation  ')
        log.disp (' using atomic units throughout, i.e. hbar = m_e = e = 1')
        log.disp (' for quantum systems interacting with electrical fields')
        log.disp ('                                                       ')
        log.disp ('    d                                                  ')
        log.disp (' i -- |psi(t)> = ( H  - F(t) mu ) |psi(t)>             ')
        log.disp ('   dt               0                                  ')
        log.disp ('                                                       ')
        log.disp (' H_0 is a diagonal matrix for the unperturbed system,  ')
        log.disp (' mu is the Hermitian matrix with dipole moments        ')
        log.disp (' and F(t) is the electric field.                       ')
        log.disp ('                                                       ')
        log.disp ('  <D>(t)= <psi(t)| D |psi(t)>                          ')
        log.disp ('                                                       ')
        
        % Initial and equilibrium state vectors/densities for propagation
        oct.tdse (real(tise.ham))
        
        % making title string for plots
        bilinear.title = 'TDSE: ' ;
        
        % Set up A; shift the eigenenergies by the ground state energy.
        bilinear.A = -1i*(diag(tise.ham)-tise.ham(1)*eye(dim));
        
        % Coupling(s) to control field(s)
        if isfield (tise, 'dip')
            for p = 1:length(tise.dip)
                if ~isempty (tise.dip{p})
                    bilinear.N{p}=tise.dip{p};
                    bilinear.B{p}=bilinear.N{p}*bilinear.x.equilib;
                end
            end
        end
        
        % Set up C vectors or D matrices
        % Choice of observables is given in control.observe.targets
        for len=1:length(control.observe.targets)
            bilinear.label{len} = tise.lab{control.observe.targets(len)};
            switch tise.obs
                case 'amo'
                    log.disp (['Observable ' int2str(len) ': Additional multiplicative operators: ' bilinear.label{len}])
                    bilinear.D{len} = tise.mat{control.observe.targets(len)};
                case 'prj'
                    log.disp (['Observable ' int2str(len) ': Populations as projectors onto eigenstates: ' bilinear.label{len}])
                    bilinear.D{len} = tise.mat{control.observe.targets(len)};
                case 'ovl'
                    log.disp (['Observable ' int2str(len) ': Populations from overlaps with eigenstates: ' bilinear.label{len}])
                    bilinear.C{len} = tise.vec{control.observe.targets(len)}';
                    bilinear.Q{len} = true; % use |<c|x>|^2 as observable
                otherwise
                    log.error ('Invalid choice of observable for TDSE')
            end
        end
        log.disp (' ')
                
        % Values of observables for initial and equilibrium state
        for len=1:length(control.observe.targets)
            switch tise.obs
                case {'amo', 'prj'}
                    bilinear.y.equilib(len) = dot ( bilinear.x.equilib, bilinear.D{len} * bilinear.x.equilib );
                    bilinear.y.initial(len) = dot ( bilinear.x.initial, bilinear.D{len} * bilinear.x.initial );
                case 'ovl'
                    bilinear.y.equilib(len) = abs ( bilinear.C{len}*bilinear.x.equilib )^2;
                    bilinear.y.initial(len) = abs ( bilinear.C{len}*bilinear.x.initial )^2;
            end
        end
        
        % Shift initial state vector w.r.t. its equilibrium state
        bilinear.x.initial = bilinear.x.initial - bilinear.x.equilib;
        
    otherwise
        log.error (['Wrong choice of equation-of-motion type : ' eom])
        
end

% Sparsity of (complex) matrix A, if density below threshold
density = nnz(bilinear.A)/numel(bilinear.A);
log.disp (['Density of matrix A: ' num2str(density)])
if density<0.1
    bilinear.A = sparse(bilinear.A);
    log.disp ('  ==> using sparse storage scheme')
else
    log.disp ('  ==> using full storage scheme')
end

% Sparsity of (real) matrices N, if density below threshold
for d=1:length(bilinear.N)
    density = nnz(bilinear.N{d})/numel(bilinear.N{d});
    log.disp (['Density of matrix N_'  int2str(d) ': ' num2str(density)])
    if density<0.1
        bilinear.N{d} = sparse(bilinear.N{d});
        log.disp ('  ==> using sparse storage scheme')
    else
        log.disp ('  ==> using full storage scheme')
    end
end

% Plot spectrum of A 
oct.spec_A (mfilename)

% Save ABNCD etc matrices to data file (index 0 stands for unbalanced, untruncated)
log.disp(['Saving matrices A, B, N, C|D and densities to file : ' eom '.mat'])
save (eom, 'bilinear');

% Output clock/date/time
log.clock;

end

