% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2007-2008 Ulf Lorenz
%
% see the README file for license details.

function obj = subsasgn(obj, index, val)
    
% We allow only setting of the parameters that are needed,
% otherwise always return an error.
if index.type == '.'
    switch index.subs
        case 'label'
            obj.label = val;
        case 'dof'
            obj.dof = val;
            if isempty(obj.label)
                obj.label = int2str(obj.dof);
            end
        case 'n_pts'
            % Check for an even number of grid points!
            if val ~= 2*round(val/2)
                util.error ('Number of grid points must be even');
            end
            obj.n_pts = val;
        case 'mass'
            if val <= 0
                util.error ('Non-positive mass assigned');
            end
            obj.mass = val;
        case 'periodic'
            obj.periodic = val;
        case 'x_min'
            if ~isempty(obj.x_max) && obj.x_max <= val
                util.error ('maximum grid point was smaller than minimum one');
            end
            obj.x_min = val;
        case 'x_max'
            if ~isempty(obj.x_min) && obj.x_min >= val
                util.error ('maximum grid point was smaller than minimum one');
            end
            obj.x_max = val;
        case 'nokin'
            if ~islogical(val) || length(val) ~= 1
                util.error('Need to supply a boolean value')
            end
            obj.nokin = val;
        otherwise
            util.error ( ['Tried to set invalid field ' index.subs 'in FFT DVR object'] );
    end
end
