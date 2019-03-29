% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2007-2008 Ulf Lorenz
%
% see the README file for license details.

function obj = subsasgn(obj, index, val)

global space

if index.type == '.'
	switch index.subs
		case 'label'
			obj.label = val;
		case 'dof'
            if val > space.n_dim
                log.error ('too large value for index "dof" assigned');
            end
			obj.dof = val;
            if isempty(obj.label)
                if space.n_dim > 1
                    obj.label = int2str(obj.dof);
                else
                    obj.label = '';
                end
            end
		case 'mass'
			if val <= 0
				log.error('non-positive mass assigned');
			end
			obj.mass = val;
		case 'R_dof'
			if ~isempty(obj.R_0)
				log.error('Legendre DVR must be either fixed rotator or given an additional degree of freedom');
			end
			obj.R_dof = val;
		case 'R_0'
			if ~isempty(obj.R_dof)
				log.error('Legendre DVR must be either fixed rotator or given an additional degree of freedom');
			end
			obj.R_0 = val;
		case 'l_max'
			if val < 0
				log.error('Cannot handle negative angular momenta in assignment');
			end
			if ~isempty(obj.m_0) && abs(obj.m_0) > val
				log.error('Error: Legendre polynomials must have l >= m');
			end
			obj.l_max = val;
		case 'm_0'
			if ~isempty(obj.l_max) && obj.l_max < abs(val)
				log.error('Error: Legendre polynomials must have l >= m');
			end
			obj.m_0 = val;
        case 'nokin'
            if ~islogical(val) || length(val) ~= 1
                log.err('Value has to be set to a single boolean.')
            end
            obj.nokin = val;
		otherwise
			log.error( ['Tried to set invalid field ' index.subs ' of Legendre DVR'] );
	end
end
