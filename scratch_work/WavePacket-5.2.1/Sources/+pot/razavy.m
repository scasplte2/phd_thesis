% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2015 Burkhard Schmidt
%
% see the README file for license details.

function razavy

global hamilt space

util.disp (' ')
util.disp ('*******************************************************')
util.disp ('Potential energy for a Razayv single/double well       ')
util.disp ('see American Journal of Physics, 48(4), 285-288 (1980) ')
util.disp ('DOI:10.1119/1.12141, Eq. (2.2) with beta=1             ')
util.disp ('                                                       ')
util.disp (' V (R) = - s*(n+1)*cosh(2*x) + s^2/8*cosh(4*x) - s^2/8 ')
util.disp ('                                                       ')
util.disp ('*******************************************************')
util.disp ( [ 'Strength parameter  s : ' num2str(hamilt.pot.params.s) ] )
util.disp ( [ 'Magic number n        : ' num2str(hamilt.pot.params.n) ] )

% Check validity
if space.size.n_dim ~= 1
    util.error ('This potential only for one dimension')
end

if hamilt.coupling.n_eqs ~= 1
    util.error ('This potential only for single Schrödinger equation')
end

s = hamilt.pot.params.s;
n = hamilt.pot.params.n;
x = space.dvr.grid_ND{1};

hamilt.pot.grid_ND{1,1} = ...
    - s*(n+1)*cosh(2*x) ...
    + s^2/8*cosh(4*x) ...
    - s^2/8;








