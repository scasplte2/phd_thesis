% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2007 Burkhard Schmidt's group
%
% see the README file for license details.

function single 
global hamilt space

util.disp (' ')
util.disp ('*****************************************************')
util.disp ('Single crossing example (Nettesheim, Schuette / 1999)')
util.disp (' ')
util.disp ('             ( A/R   C   ) ')
util.disp (' V_dia (R) = (           ) ')
util.disp ('             ( C    BR^2 ) ')
util.disp (' ')
util.disp ('*****************************************************')
util.disp ( [ 'Parameter of repulsion  A : ' num2str(hamilt.pot.params.repuls) ] )
util.disp ( [ 'Parameter of attraction B : ' num2str(hamilt.pot.params.attrac) ] )
util.disp ( [ 'Constant coupling       C : ' num2str(hamilt.pot.params.couple) ] )

% Check validity
if space.size.n_dim ~= 1
    util.error ('This potential is only for 1 dimension')
end

if hamilt.coupling.n_eqs ~= 2
    util.error ('This potential is only for 2 states')
end
    
hamilt.pot.grid_ND{1,1} = hamilt.pot.params.repuls ./ space.dvr.grid_ND{1};
hamilt.pot.grid_ND{2,2} = hamilt.pot.params.attrac  * space.dvr.grid_ND{1} .^ 2;
hamilt.pot.grid_ND{1,2} = hamilt.pot.params.couple  * ones(size(space.dvr.grid_ND{1}));



