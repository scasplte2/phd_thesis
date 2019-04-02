% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2007 Burkhard Schmidt's group
%
% see the README file for license details.

function metiu2

global hamilt space

util.disp(' ')
util.disp('**********************************************************')
util.disp('Model potential: Dateo&Metiu, J.Chem.Phys. 95, 7392 (1991)')
util.disp(' ')
util.disp('                                                      2   ')
util.disp(' V(r,\Theta) = V     (r) + 0.1*D*(1-cos \Theta)*(r /r)    ')
util.disp('                Morse                             e       ')
util.disp('**********************************************************')
util.disp([ 'Dissociation energy d_e:' num2str(hamilt.pot.params.d_e)])
util.disp([ 'range parameter a      :' num2str(hamilt.pot.params.alf)])
util.disp([ 'Equilibrium dist. r_e  :' num2str(hamilt.pot.params.r_e)])

if space.size.n_dim ~= 2 || ~strcmp(class(space.dof{1}), 'grid_fft') ...
        || ~strcmp(class(space.dof{2}), 'grid_legendre')
    util.err('This potential works only for one FFT and one Legendre DOF')
end

if hamilt.coupling.n_eqs > 1
    util.err('This potential works only for one Schroedinger equation')
end

% first apply the morse potential
hamilt.pot.grid_ND{1,1} = hamilt.pot.params.d_e * (1 - exp(-hamilt.pot.params.alf ...
            * (space.dvr.grid_ND{1} - hamilt.pot.params.r_e))).^2;

% then add the cosine \theta term
hamilt.pot.grid_ND{1,1} = hamilt.pot.grid_ND{1,1} + hamilt.pot.params.d_e * 0.1 ...
            * (1 - space.dvr.grid_ND{2}) .* (hamilt.pot.params.r_e ./ space.dvr.grid_ND{1}).^2;