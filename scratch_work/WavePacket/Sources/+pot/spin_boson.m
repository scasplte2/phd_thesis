% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008 Burkhard Schmidt
%
% see the README file for license details.

function spin_boson
global hamilt space

util.disp (' ')
util.disp ('*******************************************************')
util.disp ('Potential energy: Spin boson model                     ')
util.disp ('                                                       ')
util.disp ('           ( f(R)+kappa*R+delta/2]  gamma    )         ')
util.disp (' V   (R) = (                                 )         ')
util.disp ('  dia      ( gamma      f(R)-kappa*R-delta/2 )         ')
util.disp ('                                                       ')
util.disp ('                                                       ')
util.disp (' E   (R)=f(R) +/- sqrt[ gamma^2 + (delta/2+kappa*R)^2 ]')
util.disp ('  adi                                                  ')
util.disp ('                                                       ')
util.disp (' f(R) = omega^2/2 R^2                                  ')
util.disp ('                                                       ')
util.disp ('*******************************************************')
util.disp ([ 'Harmonic frequency omega  : ' num2str(hamilt.pot.params.omega)])
util.disp ([ 'Interstate coupling gamma : ' num2str(hamilt.pot.params.gamma)])
util.disp ([ 'Asymmetry parameter kappa : ' num2str(hamilt.pot.params.kappa)])
util.disp ([ 'Energy gap delta          : ' num2str(hamilt.pot.params.delta)])

% Check validity

if space.size.n_dim ~= 1
   util.error ('This potential only for one dimension')
end

if hamilt.coupling.n_eqs ~= 2
   util.error ('This potential only for two coupled states')
end


hamilt.pot.grid_ND{1,1} = space.dvr.grid_ND{1}.^2 * hamilt.pot.params.omega^2/2;
hamilt.pot.grid_ND{1,1} = hamilt.pot.grid_ND{1,1} + space.dvr.grid_ND{1}*hamilt.pot.params.kappa;
hamilt.pot.grid_ND{1,1} = hamilt.pot.grid_ND{1,1} + hamilt.pot.params.delta/2;

hamilt.pot.grid_ND{2,2} = space.dvr.grid_ND{1}.^2 * hamilt.pot.params.omega^2/2;
hamilt.pot.grid_ND{2,2} =hamilt.pot.grid_ND{2,2}  - space.dvr.grid_ND{1}*hamilt.pot.params.kappa;
hamilt.pot.grid_ND{2,2} =hamilt.pot.grid_ND{2,2}  - hamilt.pot.params.delta/2;

hamilt.pot.grid_ND{1,2} = hamilt.pot.params.gamma * ones(size(space.dvr.grid_ND{1}));
hamilt.pot.grid_ND{2,1} =hamilt.pot.grid_ND{1,2};







