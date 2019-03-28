% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2007 Burkhard Schmidt's group
%
% see the README file for license details.

function tully2
global hamilt space

util.disp (' ')
util.disp ('**************************************************************')
util.disp ('Dual crossing example: J.C.Tully, J.Chem.Phys. 93, 1061 (1990)')
util.disp ('**************************************************************')
util.disp ( [ 'Tully parameter A : ' num2str(hamilt.pot.params.A) ] )
util.disp ( [ 'Tully parameter B : ' num2str(hamilt.pot.params.B) ] )
util.disp ( [ 'Tully parameter C : ' num2str(hamilt.pot.params.C) ] )
util.disp ( [ 'Tully parameter D : ' num2str(hamilt.pot.params.D) ] )
util.disp ( [ 'Tully parameter E : ' num2str(hamilt.pot.params.E) ] )

% Check validity
if space.size.n_dim ~= 1
    util.error ('This potential is only for 1 dimension')
end

if hamilt.coupling.n_eqs ~= 2
    util.error ('This potential is only for 2 states')
end
    
hamilt.pot.grid_ND{1,1} = + zeros(size(space.dvr.grid_ND{1}));
hamilt.pot.grid_ND{2,2} = - hamilt.pot.params.A * exp ( - hamilt.pot.params.B * space.dvr.grid_ND{1}.^2 ) + hamilt.pot.params.E;
hamilt.pot.grid_ND{1,2} = + hamilt.pot.params.C * exp ( - hamilt.pot.params.D * space.dvr.grid_ND{1}.^2 );



