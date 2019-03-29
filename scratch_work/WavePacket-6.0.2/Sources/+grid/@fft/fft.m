% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%               2017 Burkhard Schmidt
%
% see the README file for license details.

classdef fft < handle
    
    properties (Access = public)
        
        label       % labelling the degree of freedom
        dof         % which degree of freedom
        mass        % mass that enters the kinetic energy
        periodic    % use periodic boundary conditions or not
        nokin       % enable/disable the kinetic energy operator
        
        n_pts       % number of grid points (equal in position and momentum grid)
        x_min       % lower bound of the position grid
        x_max       % upper bound of the position grid
        
        weight      % weights in DVR ("position space")
        x_grid      % grid points in DVR ("position space")
        p_grid      % grid points in FBR ("momentum space")

        x_dlt       % (constant!) grid spacing in position space
        p_dlt       % (constant!) grid spacing in momentum space
        
    end
    
    properties (Access = private)
        
        p_min       % lower bound of the momentum grid
        p_max       % upper bound of the momentum grid
        
        kin         % DVR grid representation of the kinetic energy operator
        kin_shift   % an fftshifted version of kin
        kin_expo    % same as grid, but exponentiated, for split operator
        kin_factor  % normalize and removes a spurious complex phase
        
    end
    
    methods (Access = public) % in separate files within *this* directory
        
        dvr = fbr2dvr(obj, fbr)
        fbr = dvr2fbr(obj, dvr)
        init_grid  (obj)
        init_kin   (obj, fraction, output)
        kinetic    (obj, psi, new)
        kinetic_exp(obj, psi)
        dvrkin = kinetic2dvr(obj)
        dvr = matrix2dvr(obj, fbr)
        retval = momentum(obj, psi)
        obj = subsasgn(obj, index, val)
        retval = subsref(obj, s)
        
        % Constructor method: Initialize default property values
        function obj = fft
            obj.dof = 1;
            obj.mass = 1;
            obj.periodic = true;
            obj.nokin = false;
        end
        
        function disp (obj) % Diplay method, overloading default disp method
            % Works only if init_grid has been called previously
            log.disp ( '***************************************************************')
            log.disp ( [ 'DVR for the degree of freedom: ' obj.label] )
            log.disp ( 'Discretization scheme: Equally spaced grids (FFT) ')
            log.disp ( '***************************************************************')
            log.disp ( ' ' )
            log.disp ( [ 'Number of grid points  : ' num2str(obj.n_pts) ] )
            log.disp ( ' ' )
            log.disp ( 'Position space' )
            log.disp ( [ 'Minimum of grid values : ' num2str(obj.x_min) ] )
            log.disp ( [ 'Maximum of grid values : ' num2str(obj.x_max) ] )
            log.disp ( [ 'Spacing of grid values : ' num2str(obj.x_dlt) ] )
            log.disp ( ' ' )
            log.disp ( 'Momentum (wavenumber) space (Fourier transform)' )
            log.disp ( [ 'Maximum of grid values : ' num2str(obj.p_max) ] )
            log.disp ( [ 'Spacing of grid values : ' num2str(obj.p_dlt) ] )
            log.disp ( ' ' )
            log.disp ( 'Momentum (wavenumber) space (Wigner transform)' )
            log.disp ( [ 'Maximum of grid values : ' num2str(obj.p_max/2) ] )
            log.disp ( [ 'Spacing of grid values : ' num2str(obj.p_dlt/2) ] )
            log.disp ( ' ' )
            
        end
        
        % Evaluate kinetic energy function (for trajectories!)
        function eval_kin (obj, state)
            state.kin = state.kin + state.mom{obj.dof}.^2 / (2*obj.mass);
        end
        
    end
    
end
