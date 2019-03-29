% Like "normal" Chebychev imaginary time propagator, but apply the operator
% (1-P) H (1-P) instead of H, where P is the projection on the ground state.
% It can be calculated by expanding a wavefunction
% |Psi> = c_0 |GS> + |Rest>,      (1-P) |Psi> = |Rest>
% with c_0 = <GS|Psi>
% => |Rest> = |Psi> - c_0 |GS>

% Copyright (C) 2017 Burkhard Schmidt
%               2008 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

classdef my_relax_2 < handle
    
    properties (Access = public)
        order           % order of polynomial expansion
        precision       % used for truncating the polynomial series
        automatic       % toggle
    end
    
    properties (Access = private)
        gs1
        gs2
         
        main_coeffs
        main_expon
        
   end
    
    methods (Access = public)
        
        % Construct propagator: Set default values
        function obj = my_relax_2
            obj.order = 0;
            obj.precision = eps;
            obj.automatic = false;
        end
        
        % Display propagator, overloading default disp method
        function disp(obj)
            global time
            log.disp ('Chebychev polynomial expansion in imaginary time')
            log.disp ('***************************************************************')
            log.disp (' ')
            log.disp (['Kosloff number dE*dt / (2*hbar) : ' num2str(time.steps.m_alpha)])
            if obj.automatic
                log.disp (['Automatic truncation if coefficients drop below : ' num2str(obj.precision)])
                log.disp (['Number of polynomials required : ' int2str(obj.order)])
            else
                log.disp (['Truncating polynomial series after : ' int2str(obj.order)])
            end
        end
        
        % Initialize propagator
        function wave_init (obj,psi)
            global gs hamilt space time
            
            % Set/check parameters
            if isfield(time,'pulse')
                log.error ('This propagator is only for time independent Hamiltonians (no electric fields)')
            end
            
            if time.steps.m_alpha <= 10
                log.error ('Kosloff number should not be below 10, then it becomes inefficient.')
            end
            
            % Automatically determine expansion order
            if obj.order == 0
                obj.automatic = true;
                
                % Search backward
                for ii = round(10*sqrt(time.steps.m_alpha)) : -1 : round(sqrt(time.steps.m_alpha))
                    if abs(besseli(ii,time.steps.m_alpha))>obj.precision
                        break
                    end
                end
                obj.order = ii;
            end
            
            % Use MODIFIED Bessel functions to get expansion coefficients
            obj.main_coeffs = besseli(0:obj.order,time.steps.m_alpha);
            obj.main_coeffs(2:obj.order+1) = 2 * obj.main_coeffs(2:obj.order+1);
            
            % Phase factors to compensate for normalization of Hamiltonian
            obj.main_expon = exp( -(hamilt.range.delta/2+hamilt.range.e_min) * time.steps.m_delta );
            
            % Ground and first excited state (from previous simulation)
            obj.gs1 = gs{1};
            obj.gs2 = gs{2};
            
            coeff = sum( conj(obj.gs1(:)) .* psi.dvr{1}(:) .* space.weight(:) );
            psi.dvr{1} = psi.dvr{1} - coeff * obj.gs1;
            
            coeff = sum( conj(obj.gs2(:)) .* psi.dvr{1}(:) .* space.weight(:) );
            psi.dvr{1} = psi.dvr{1} - coeff * obj.gs2;
            
            norm2 = sum(abs(psi.dvr{1}(:)).^2 .* space.weight(:));
            psi.dvr{1} = psi.dvr{1} / sqrt(norm2);
            
        end
        
        % Perform propagation
        function wave_propa (obj,psi)
            
            global hamilt space
            
            % Pre-allocate
            cheby0 = cell(hamilt.coupling.n_eqs,1);
            cheby1 = cell(hamilt.coupling.n_eqs,1);
            cheby2 = cell(hamilt.coupling.n_eqs,1);

            %-----------------------------------------------------------
            %  Zero-th Chebyshev polynomial : phi_0 = 1
            %-----------------------------------------------------------
            coeff = sum( conj(obj.gs1(:)) .* psi.dvr{1}(:) .* space.weight(:) );
            psi.dvr{1} = psi.dvr{1} - coeff * obj.gs1;
            coeff = sum( conj(obj.gs2(:)) .* psi.dvr{1}(:) .* space.weight(:) );
            psi.dvr{1} = psi.dvr{1} - coeff * obj.gs2;
            for m = 1:hamilt.coupling.n_eqs
                cheby0{m} = psi.dvr{m};
                psi.sum{m} = obj.main_coeffs(1) * cheby0{m};
            end
    
            %-----------------------------------------------------------
            %  First Chebychev polynomial phi_1 = - Hnorm
            %-----------------------------------------------------------
            
            % Projection before and after
            coeff = sum( conj(obj.gs1(:)) .* psi.dvr{1}(:) .* space.weight(:) );
            psi.dvr{1} = psi.dvr{1} - coeff * obj.gs1;
            coeff = sum( conj(obj.gs2(:)) .* psi.dvr{1}(:) .* space.weight(:) );
            psi.dvr{1} = psi.dvr{1} - coeff * obj.gs2;
            apply_ham(psi,[0 0],1);
            coeff = sum( conj(obj.gs1(:)) .* psi.new{1}(:) .* space.weight(:) );
            psi.new{1} = psi.new{1} - coeff * obj.gs1;
            coeff = sum( conj(obj.gs2(:)) .* psi.new{1}(:) .* space.weight(:) );
            psi.new{1} = psi.new{1} - coeff * obj.gs2;
            
            for m = 1:hamilt.coupling.n_eqs
                cheby1{m} = - psi.new{m};
                psi.sum{m} = psi.sum{m} + obj.main_coeffs(2) * cheby1{m};
            end
    
            
            %-----------------------------------------------
            %  Higher Chebychev polynomials (n>1) by recursion:
            %  phi_n = -2*Hnorm phi_{n-1} - phi_{n-2}
            %-----------------------------------------------
            for k=2:obj.order
                
                for m = 1:hamilt.coupling.n_eqs
                    psi.dvr{m} = cheby1{m};
                end
                coeff = sum( conj(obj.gs1(:)) .* psi.dvr{1}(:) .* space.weight(:) );
                psi.dvr{1} = psi.dvr{1} - coeff * obj.gs1;
                coeff = sum( conj(obj.gs2(:)) .* psi.dvr{1}(:) .* space.weight(:) );
                psi.dvr{1} = psi.dvr{1} - coeff * obj.gs2;
                apply_ham(psi,[0 0],1);
                coeff = sum( conj(obj.gs1(:)) .* psi.new{1}(:) .* space.weight(:) );
                psi.new{1} = psi.new{1} - coeff * obj.gs1;
                coeff = sum( conj(obj.gs2(:)) .* psi.new{1}(:) .* space.weight(:) );
                psi.new{1} = psi.new{1} - coeff * obj.gs2;
                for m = 1:hamilt.coupling.n_eqs
                    cheby2{m} = - 2 * psi.new{m} - cheby0{m};
                    psi.sum{m} = psi.sum{m} + obj.main_coeffs(k+1) * cheby2{m};
                    
                    cheby0{m} = cheby1{m};
                    cheby1{m} = cheby2{m};
                end
                
            end
            
            %  Multiply wave function with exponential factor
            for m = 1:hamilt.coupling.n_eqs
                psi.dvr{m} = psi.sum{m} * obj.main_expon;
            end
            
            %  Re-normalize wave function
            fac = 0;
            for m = 1:hamilt.coupling.n_eqs
                fac = fac + sum ( abs(psi.dvr{m}(:)).^2 .* space.weight(:) );
            end
            
            for m = 1:hamilt.coupling.n_eqs
                psi.dvr{m}(:) = psi.dvr{m}(:) / sqrt(fac);
            end
                        
        end
    end
end