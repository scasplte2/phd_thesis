% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2007-2008 Ulf Lorenz
%
% see the README file for license details.

function retval = subsref(obj, s)

% For now, allow access to all internal data
if s.type == '.'
    switch s.subs
        % first the things that are required from all grid classes
        case 'label'
            retval = obj.label;
        case 'dof'
            retval = obj.dof;
        case 'n_pts'
            retval = obj.n_pts;
        case 'dvr_min'
            retval = obj.x_min;
        case 'dvr_max'
            retval = obj.x_max;
        case 'fbr_min'
            retval = obj.p_min;
        case 'fbr_max'
            retval = obj.p_max;
        case 'kin_max'
            kmax = pi * obj.n_pts / (obj.x_max - obj.x_min);    % maximum k vector
            retval = kmax^2 / (2*obj.mass);                     % kinetic energy
	    if obj.nokin
	        retval = 0;
	    end
        case 'nokin'
            retval = obj.nokin;

        % Now the things that are public
        case 'mass'
            retval = obj.mass;
        case 'x_min'
            retval = obj.x_min;
        case 'x_max'
            retval = obj.x_max;
        case 'periodic'
            retval = obj.periodic;

        % Some derived stuff
        case 'kin'
            retval = obj.kin;
        otherwise
            util.error ( ['Tried to get invalid field ' s.subs 'in FFT grid object'] );
    end
else
    util.error( 'Tried unsupported access form of class grid_fft' );
end
