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
            retval = obj.l_max - abs(obj.m_0) + 1;
        case 'dvr_min'
            retval = -1;
        case 'dvr_max'
            retval = +1;
        case 'fbr_min'
            retval = abs(obj.m_0);
        case 'fbr_max'
            retval = obj.l_max;
        case 'kin_max'
            % energy is l(l+1)/2mR^2
            if ~isempty(obj.R_0)
                Rmin = obj.R_0;
            else
                Rmin = min(space.dvr.grid_1D{obj.R_dof}(:));
            end
            retval = obj.l_max * (obj.l_max+1) / (2*obj.mass*Rmin^2);
	    if obj.nokin
	        retval = 0;
	    end

        % Now our public data
		case 'mass'
			retval = obj.mass;
		case 'R_dof'
			retval = obj.R_dof;
		case 'R_0'
			retval = obj.R_0;
		case 'l_max'
			retval = obj.l_max;
		case 'm_0'
			retval = obj.m_0;

        % and some internal data that probably noone has a use for,
        % but let us hand it out anyway.
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
			util.error( ['unknown field ' index.subs ' of grid_legendre requested']);
	end
else
	util.error('Tried unsupported access form for class grid_legendre');
end
