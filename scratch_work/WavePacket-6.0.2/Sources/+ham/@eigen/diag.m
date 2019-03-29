%-------------------------------------------------------------------------
%
% Solve the time-independent Schroedinger equation to get eigenstates
% and energies in pseudospectral representation using DVR/FBR techniques 
% 
% Part 3/3: Diagonalizing the Hamiltonian
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2011 Ulf Lorenz
%               2011 Ulf Lorenz
%
% see the README file for license details.

function diag (obj)

global space

%% Find eigenvalues and eigenvectors of Hamiltonian matrix; Sort them!
tim = cputime;

if obj.storage == 'f'
    log.disp('Start diagonalizing full matrix ...')
    [ obj.eig_vecs, obj.eig_vals ] = eig ( obj.matrix );
else
    log.disp(['Density of matrix :' num2str(nnz(obj.matrix) / space.n_tot^2)])
    log.disp('Start diagonalizing sparse matrix ...')
    opts.disp = 0;
    [ obj.eig_vecs, obj.eig_vals ] = ...
        eigs ( obj.matrix, obj.stop+1, 'sm', opts );
end

log.disp('Sorting eigenvalues...')

%[sorted, order] = sort(diag(obj.eig_vals),'ComparisonMethod','real');
[sorted, order] = sort(real(diag(obj.eig_vals)));
obj.eig_vals = sorted;
sortvecs = zeros(size(obj.eig_vecs));
for k= 1:length(order)
    sortvecs(:, k) = obj.eig_vecs(:, order(k));
end
obj.eig_vecs = sortvecs;


log.disp (['Finished after [CPU seconds] : ' num2str(cputime-tim)])


%% List eigenvalues
log.disp(' ')
log.disp('Table of eigenvalues')
for ii=obj.start : obj.stop
    energy = obj.eig_vals(ii+1);
    log.disp([ int2str(ii) ' : ' num2str(math.real(energy))])
end
log.disp(' ')
