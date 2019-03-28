%--------------------------------------------------------------------------
%
% Negative imaginary potential (NIP) to be used as 
% (smoothly) absorbing boundary conditions
%
% Power function: (x-x0)^n
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2011 Ulf Lorenz
%
% see the README file for license details.

function power
global hamilt space

util.disp (' ')
util.disp ( '**************************************************' )
util.disp ( 'Absorbing boundary conditions' )
util.disp ( '**************************************************' )
util.disp ( 'Negative imaginary potential  : Power function' )
util.disp ( ' ' )
util.disp ( [ 'Exponent of power function : ' int2str(hamilt.nip.params.exp) ] )
util.disp ( ' ' )
util.disp ( [ 'Beginning of inner grid region : ' num2str(hamilt.nip.params.min) ] )
util.disp ( [ 'End of inner grid region       : ' num2str(hamilt.nip.params.max) ] )
util.disp ( '   ' )

% Check input parameters
if hamilt.nip.params.exp<1
    util.error ( 'Exponent of negative imaginary potential must be at least 1' )
end

if length(hamilt.nip.params.exp)~=space.size.n_dim
    util.error ('Incompatible dimensionality for NIP exponenents')
end

if length(hamilt.nip.params.min)~=space.size.n_dim
    util.error ('Incompatible dimensionality for beginning of inner grid')
end

if length(hamilt.nip.params.max)~=space.size.n_dim
    util.error ('Incompatible dimensionality for end of inner grid')
end

if any (hamilt.nip.params.max <= hamilt.nip.params.min)
    util.error ( 'Wrong ordering of "NIP" min/max parameters' )
end

hamilt.nip.grid = zeros(size(space.dvr.grid_ND{1}));

% Small parameter: exp(-NIP)=epsilon at grid boundaries
epsilon = eps('single');

% Multiplying contributions from each component of position vector
for k = 1:space.size.n_dim

    % Lower "gobbler" region
    mask = find ( space.dvr.grid_ND{k} < hamilt.nip.params.min(k) );
    if ~isempty ( mask )
        factor = -log(epsilon)/(hamilt.nip.params.min(k)-space.dof{k}.dvr_min)^hamilt.nip.params.exp(k);
        hamilt.nip.grid(mask) = hamilt.nip.grid(mask) + factor * ...
            (hamilt.nip.params.min(k)-space.dvr.grid_ND{k}(mask)).^hamilt.nip.params.exp(k);
    end

    % Upper "gobbler" region
    mask = find ( space.dvr.grid_ND{k} > hamilt.nip.params.max(k) );
    if  ~isempty ( mask )
        factor = -log(epsilon)/(space.dof{k}.dvr_max-hamilt.nip.params.max(k))^hamilt.nip.params.exp(k);
        hamilt.nip.grid(mask) = hamilt.nip.grid(mask) + factor * ...
            (space.dvr.grid_ND{k}(mask)-hamilt.nip.params.max(k)).^hamilt.nip.params.exp(k);
    end

end
