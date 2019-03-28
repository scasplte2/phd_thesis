%--------------------------------------------------------------------------
%
% Gaussian bell-shaped projection function (PRJ)
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2007 Burkhard Schmidt's group
%
% see the README file for license details.

function gauss

global space

util.disp (' ')
util.disp ( '*************************************************' )
util.disp ( 'Gaussian bell-shaped projection function         ' )
util.disp ( '                                                 ' )
util.disp ( '           N      [   ( Ri-R0i )^2 ]             ' )
util.disp ( 'prj(R) = Prod exp [ - (--------)   ]             ' )
util.disp ( '          i=1     [   (  2*Wi  )   ]             ' )
util.disp ( '                                                 ' )
util.disp ( 'where the product extends over all dimensions    ' )
util.disp ( '*************************************************' )
util.disp ( [ 'Mean value position       R0 : ' num2str(     space.prj.params.pos_0) ] )
util.disp ( [ 'Position uncertainty      W  : ' num2str(     space.prj.params.width) ] )

% Allocate and initialize
space.prj.grid = ones ( size(space.dvr.grid_ND{1}) );

% Tensor product of one-dimensional Gaussians
for k = 1:space.size.n_dim
    space.prj.grid_ND = space.prj.grid_ND .* exp (  ...
        -((space.dvr.grid_ND{k}-space.prj.params.pos_0(k)) / (space.prj.params.width(k)*2)).^2  );
end


