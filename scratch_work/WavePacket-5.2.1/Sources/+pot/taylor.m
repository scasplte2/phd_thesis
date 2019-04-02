% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2016 Burkhard Schmidt
%
% see the README file for license details.

function taylor 
global hamilt space

util.disp (' ')
util.disp ('*******************************************************')
util.disp ('Potential energy: Taylor series (diag. in N dimensions)')
util.disp ('                                                       ')
util.disp ('          inf   N   v_mn,jk              j             ')
util.disp (' V  (R) = Sum  Sum  ------- ( R  - s    )   + c        ')
util.disp ('  mn      j=1  k=1     j!      k    mn,k        mn     ')
util.disp ('                                                       ')
util.disp ('*******************************************************')
if ~isfield(hamilt.pot,'params')
    hamilt.pot.params = [];
end
if ~isfield(hamilt.pot.params,'s')
    hamilt.pot.params.s = cell (hamilt.coupling.n_eqs);
end
if ~isfield(hamilt.pot.params,'c')
    hamilt.pot.params.c = cell (hamilt.coupling.n_eqs);
end
if ~isfield(hamilt.pot.params,'v')
    hamilt.pot.params.v = cell (hamilt.coupling.n_eqs);
end

%Diagonal and offdiagonal (coupling) elements of Hamiltonian ( m <= n )
for m = 1:hamilt.coupling.n_eqs
    for n = m:hamilt.coupling.n_eqs
        util.disp (['Taylor series for potential(' int2str(m) ',' int2str(n) ')'])
        
        % Reference position: default is zero
        if isempty(hamilt.pot.params.s{m,n})
            hamilt.pot.params.s{m,n} = zeros(1,space.size.n_dim); % row vector
        end
        
        % Constant energy offset; default is zero
        if isempty(hamilt.pot.params.c{m,n})
            hamilt.pot.params.c{m,n} = 0; % scalar
        end
        
        % Evaluate Taylor series only for desired {m,n} combinations
        if isempty (hamilt.pot.params.v{m,n})
           hamilt.pot.grid_ND{m,n} = ...
               hamilt.pot.params.c{m,n} * ...
               ones(size(space.dvr.grid_ND{1}));
        else
            hamilt.pot.grid_ND{m,n} = util.taylor (...
                space.dvr.grid_ND, ...
                hamilt.pot.params.s{m,n}, ...
                hamilt.pot.params.c{m,n}, ...
                hamilt.pot.params.v{m,n} ); 
        end
        
       util.disp ('   ')
    end
end