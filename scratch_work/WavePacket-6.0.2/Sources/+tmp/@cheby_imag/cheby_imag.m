%----------------------------------------------------------------------
%
% Imaginary time Chebychev propagator for a quantum state vector
% ==============================================================
%
% By introducing an imaginary time variable, tau = i*t, the
% time-dependent Schrödinger equation can be cast into the 
% form of a diffusion equation. Then the propagation causes 
% an initially guessed wave function to relax to the lowest 
% eigenfunction, i.e., the quantum mechanical ground state.
%
% R. Kosloff and H. Tal-Ezer, Chem. Phys. Lett. 127(3), 223, (1986)
%
%                   (  de          )        (            )
% exp(-H*tau) = exp (-(--+e   )*tau)  * exp (-alpha*H    )
%                   (  2   min     )        (        norm)
%
% where
%         2    (       de       )
% H     = -- * (H - I*(--+e   ) )
%  norm = de   (       2   min  ) 
%
%          de * tau
% alpha  = --------        
%             2
%
% TAU     imaginary time step
% EMIN    minimum of the Hamiltonian
% DE      range of the Hamiltonian
% ALPHA   dimensionless parameter 
% Hnorm   Hamiltonian normalized to [-1,1]
% H       original Hamiltonian 
% I       unity operator
%
% Then the exponential of the time evolution operator 
% is expanded in a series of real Chebychev polynomials.
%
%     (             )     N
% exp ( -alpha*H    )  = Sum c (alpha) * phi (-H    )
%     (         norm)    n=0  n             n   norm 
%  
% 
% where the coefficients c_n are modified Bessel functions 
% and where the psi_n are the real Chebychev polynomials
% which are calculated using the recursion 
%
%    phi  = I   (unity)
%       0
%
%    phi  = -H
%       1     norm
%
%    phi  = -2*H     phi    - phi    , n>1
%       n       norm    n-1      n-2
%
%------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%
% see the README file for license details.

classdef cheby_imag < handle
    
    properties (Access = public)
        order           % order of polynomial expansion
        precision       % used for truncating the polynomial series
        automatic       % toggle
    end
    
    properties (Access = private)
        main_coeffs
        main_expon
    end
    
    methods (Access = public)
        
        % Construct propagator: Set default values
        function obj = cheby_imag
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
        function wave_init (obj, ~)
            global hamilt time
            
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
            
            % Factors to compensate for normalization of Hamiltonian
            obj.main_expon = exp( -(hamilt.range.delta/2+hamilt.range.e_min) * time.steps.m_delta );
            
            % Typical values for the wave function are of order 1 (give or take a few
            % orders of magnitude). The typical values of the Chebychev polynomials that
            % we are juggling around range from about 1 to 1/main_expon. Since the
            % largest absolute value that can be represented by doubles is around
            % 10^320, we have to check that we do not play around with too large
            % numbers. I use a very generous offset here, because numerical errors seem
            % to appear already earlier.
            if obj.main_expon < 1e-200
                log.error('Time step is too large. Chebychev expansion might not converge.');
            end
            
        end
        
        % Perform propagation
        function wave_propa (obj, psi)
            
            global hamilt space
            
            % Pre-allocate
            cheby0 = cell(hamilt.coupling.n_eqs,1);
            cheby1 = cell(hamilt.coupling.n_eqs,1);
            cheby2 = cell(hamilt.coupling.n_eqs,1);
            
            %-----------------------------------------------------------
            %  Zero-th Chebyshev polynomial : phi_0 = 1
            %-----------------------------------------------------------
            for m = 1:hamilt.coupling.n_eqs
                cheby0{m} = psi.dvr{m};
                psi.sum{m} = obj.main_coeffs(1) * cheby0{m};
            end
            
            
            %-----------------------------------------------------------
            %  First Chebychev polynomial phi_1 = - Hnorm
            %-----------------------------------------------------------
            apply_ham(psi,[0 0],1);
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
                apply_ham(psi,[0 0],1);
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
