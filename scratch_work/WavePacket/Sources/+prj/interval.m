%--------------------------------------------------------------------------
%
% Rectangular interval projection function (PRJ)
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2007 Burkhard Schmidt's group
%
% see the README file for license details.

function interval

global space

util.disp (' ')
util.disp ( '****************************************' )
util.disp ( 'Rectangular interval projection function' )
util.disp ( '****************************************' )
util.disp ( [ 'Beginning of projection interval : ' num2str(space.prj.params.min) ] )
util.disp ( [ 'End of projection interval       : ' num2str(space.prj.params.max) ] )


% Check input parameters
if length(space.prj.params.min)~=space.size.n_dim
    util.error ('Incompatible dimensionality for beginning of projection')
end

if length(space.prj.params.max)~=space.size.n_dim
    util.error ('Incompatible dimensionality for end of projection')
end

if any (space.prj.params.max <= space.prj.params.min)
    util.error ( 'Wrong ordering of projection interval min/max parameters' )
end

space.prj.grid_ND = ones(size(space.dvr.grid_ND{1}));

% Multiplying contributions from each component of position vector
for k = 1:space.size.n_dim

    % Find grid points outside interval and set PRJ function to zero
    space.prj.grid_ND(space.dvr.grid_ND{k} < space.prj.params.min(k) ...
        | space.dvr.grid_ND{k} > space.prj.params.max(k)) = 0;

end
