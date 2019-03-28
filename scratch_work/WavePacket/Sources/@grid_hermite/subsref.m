% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2007-2008 Ulf Lorenz
%
% see the README file for license details.

function retval = subsref(obj, index)

global space

if index.type == '.'
	switch index.subs
        % first those attributes that are required from all grid classes
		case 'label'
			retval = obj.label;
		case 'dof'
			retval = obj.dof;
        case 'n_pts'
            retval = obj.n_pts;
        case 'dvr_min'
            retval = min(space.dvr.grid_1D{obj.dof}(:));
        case 'dvr_max'
            retval = max(space.dvr.grid_1D{obj.dof}(:));
        case 'fbr_min'
            retval = 0;
        case 'fbr_max'
            retval = obj.n_pts - 1;
        case 'kin_max'
            if obj.nokin
                retval = 0;
            elseif obj.kin_max > 0
                retval = obj.kin_max;
            else
                % energy is omega * (n+0.5), which is not quite correct, since
                % we only use the_kinetic_ energy operator internally, but, well...
                % in lack of better numbers
                retval = obj.omega * (obj.n_pts - 0.5);
            end

        % some of the public data
        case 'mass'
            retval = obj.mass;
        case 'omega'
            retval = obj.omega;
        case 'v_2'
            retval = obj.v_2;
        case 'r_e'
            retval = obj.r_e;

        % some internal data that probably noone has a use for,
        % but let us hand it out anyway.
        case 'momentum'
            retval = obj.momentum;
		case 'kin'
			retval = obj.kin;
		case 'kinexpo'
			retval = obj.kinexpo;
        case 'trafo_expand'
            retval = obj.trafo_expand;
        case 'trafo_reconstruct'
            retval = obj.trafo_reconstruct;
        case 'nokin'
            retval = obj.nokin;
		otherwise
			util.error( ['unknown field ' index.subs ' of grid_hermite requested']);
	end
else
	util.error('Tried unsupported access form for class grid_hermite');
end
