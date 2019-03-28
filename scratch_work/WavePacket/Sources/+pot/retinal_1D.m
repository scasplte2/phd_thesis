% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2007 Burkhard Schmidt's group
%
% see the README file for license details.

function retinal_1D
global hamilt space

util.disp (' ')
util.disp ('*********************************************************************')
util.disp ('Retinal isomerization: S.Hahn, G.Stock, JPC B 104(6),1149 (2000)')
util.disp ('*********************************************************************')
util.disp ( ' ' )
util.disp ( 'Ground state ' )
util.disp ( [ 'Vertical shift    : ' num2str(hamilt.pot.params.shift(1)) ] )
util.disp ( [ 'Barrier height    : ' num2str(hamilt.pot.params.barrier(1)) ] )
util.disp ( [ 'Frequency         : ' num2str(hamilt.pot.params.omega(1)) ] )
util.disp ( [ 'Linear term       : ' num2str(hamilt.pot.params.kappa(1)) ] )
util.disp ( ' ' )
util.disp ( 'Excited state ' )
util.disp ( [ 'Vertical shift    : ' num2str(hamilt.pot.params.shift(2)) ] )
util.disp ( [ 'Barrier height    : ' num2str(hamilt.pot.params.barrier(2)) ] )
util.disp ( [ 'Frequency         : ' num2str(hamilt.pot.params.omega(2)) ] )
util.disp ( [ 'Linear term       : ' num2str(hamilt.pot.params.kappa(2)) ] )
util.disp ( ' ' )
util.disp ( [ 'Coupling mode     : ' num2str(hamilt.pot.params.xc) ] )
util.disp ( [ 'Coupling constant : ' num2str(hamilt.pot.params.lambda) ] )

% Check validity
if space.size.n_dim ~= 1
    util.error ('This potential is only for 1 dimension')
end

if hamilt.coupling.n_eqs ~= 2
    util.error ('This potential is only for 2 states')
end


hamilt.pot.grid_ND{1,1} = hamilt.pot.params.shift(1) ...
           + hamilt.pot.params.barrier(1) * (1 - cos(space.dvr.grid_ND{1})) / 2 ...
           + hamilt.pot.params.omega(1)*hamilt.pot.params.xc^2/2 ...
           + hamilt.pot.params.kappa(1)*hamilt.pot.params.xc;
       
hamilt.pot.grid_ND{2,2} = hamilt.pot.params.shift(2) ...
           + hamilt.pot.params.barrier(2) * (1 - cos(space.dvr.grid_ND{1})) / 2 ...
           + hamilt.pot.params.omega(2)*hamilt.pot.params.xc^2/2 ...
           + hamilt.pot.params.kappa(2)*hamilt.pot.params.xc;
       
hamilt.pot.grid_ND{1,2} = hamilt.pot.params.lambda * hamilt.pot.params.xc * ones(size(space.dvr.grid_ND{1}));



