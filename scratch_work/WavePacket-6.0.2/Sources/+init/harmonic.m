%------------------------------------------------------------------------------
%
% This class creates the initial state 
% as an eigenstate of a Harmonic oscillator.
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2007-2009 Ulf Lorenz
%               2008-2009 Burkhard Schmidt
%
% see the README file for license details.

classdef harmonic < init.generic & handle
    
    properties (Access = public)

        m_r         % reduced mass
        r_e         % equilibtrium distance
        n_q         % quantum number

        omega       % harmonic frequency
        v_2         % force constant
                
    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = harmonic
            obj.r_e = 0;
            obj.n_q = 0;
        end
        
        % Initialize pendular functions: Set/check parameters
        function init_dof (obj)
            global space
            
            % Default: Mass from the grid setting along the respective degree of freedom
            if ~isfield(obj, 'm_r')
                obj.m_r = space.dof{obj.dof}.mass;
            end
            
            % Check quantum number
            if round(obj.n_q)~=obj.n_q
               log.error ('Quantum number should be integer')
            end
            if (obj.n_q<0)
                log.error ('quantum number must not be negative')
            end
            
            % Use either v_2 or omega from the input file; calculate missing quantity
            if isfield(obj, 'v_2') && ~isfield(obj, 'omega')
                obj.omega = sqrt(obj.v_2 / obj.m_r);
            end
            if isfield(obj, 'omega') && ~isfield(obj, 'v_2')
                obj.v_2 = obj.omega^2 * obj.m_r;
            end
            
        end

       % Display Gaussian, overloading default disp method
       function disp(obj)
           log.disp ('Harmonic oscillator eigenstat')
           log.disp ('***************************************************************')
           log.disp ('   ' )
           log.disp ('Harmonic oscillator eigenstate')
           log.disp ( ['Angular Frequency            : ' num2str(obj.omega)] )
           log.disp ( ['corresp. Force constant      : ' num2str(obj.v_2)] )
           log.disp ( ['(Reduced) mass               : ' num2str(obj.m_r)] )
           log.disp ( ['Equilibrium position         : ' num2str(obj.r_e)] )
           log.disp ( ['Quantum number (eigenstate)  : ' int2str(obj.n_q)] )
       end
           
        % Evaluate  wave function on a grid ==> Q/M propagation
        function wave_dof (obj)
            global space

            % HO eigenstate \Psi(x) = exp(-m\omega/2 * x^2) * H_n(\sqrt(m\omega) * x)
            factor    = obj.m_r * obj.omega;
            position  = space.dvr{dir} - obj.r_e;
            obj.dvr = exp(- factor/2 * position.^2) ...
                .* math.hermite(sqrt(factor) * position, obj.n_q);
            
        end
        
        % Sample phase space density ==> Q/C propagation
        function traj_dof (obj, traj)
            log.error ('Code for phase space sampling still missing')
        end
        
    end
    
end
