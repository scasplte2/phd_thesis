% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2008,2011 Ulf Lorenz
%
% see the README file for license details.

function init_kin (obj, fraction, output)

% initialises the kinetic energy

global hamilt space time

if nargin < 3
    output = true;
end

%% Checks

if isempty( obj.mass_R ) || isempty( obj.mass_r )
     log.error ('masses not set')
end
if isempty( obj.dof_c ) 
     log.error ('c degree of freedom not set')
end
if isempty( obj.dof_R ) 
     log.error ('R degree of freedom not set')
end
if isempty( obj.dof_r ) && isempty( obj.r_0 )
     log.error ('r degree of freedom not set')
end

% Check that we have a reasonable grid (i.e. Legendre polynomials)
if ~isa ( space.dof{obj.dof_c}, 'grid.legendre' )
     log.error ('Jacobi kinetic energy only works for Legendre grids')
end


%% Informational output

if output
    log.disp ('***************************************************************')
    log.disp ('Kinetic energy: Triatomic molecule ABC            ')
    log.disp ('***************************************************************')
    log.disp ('                                                  ')
    log.disp ('           [     1             1    ] ^ 2         ')
    log.disp (' T (c) = - [ --------  +   -------- ] L           ')
    log.disp ('           [ 2 M R^2       2 m r^2  ]             ')
    log.disp ('                                                  ')
    log.disp ('where R, M is the distance and reduced mass of BC,')
    log.disp ('r distance of A to CMS of BC, m the reduced mass  ')
    log.disp ('of A and BC. See eg J. Chem. Phys. 116,4403 (2002)')
    log.disp ('                                                  ')
    log.disp ( [ 'M     : ' num2str(obj.mass_R) ] )
    log.disp ( [ 'm     : ' num2str(obj.mass_r) ] )
    log.disp ( [ 'DOF c : ' num2str(obj.dof_c)  ] )
    log.disp ( [ 'DOF R : ' num2str(obj.dof_R)  ] )
    if ~isempty(obj.dof_r)
        log.disp ( [ 'DOF r : ' num2str(obj.dof_r) ] )
    else
        log.disp ( [ 'constant r : ' num2str(obj.r_0) ] )
    end
    log.disp (' ')
end


%% Create all the grids

% Prefactor
if isempty( obj.r_0 )
    obj.grid = 1 ./ (2 * space.dvr{obj.dof_R}.^2 * obj.mass_R) ...
               + 1 ./ (2 * space.dvr{obj.dof_r}.^2 * obj.mass_r);
else
    obj.grid = 1 ./ (2 * space.dvr{obj.dof_R}.^2 * obj.mass_R) ...
               + 1 ./ (2 * obj.r_0.^2 * obj.mass_r);
end

% L^2
obj.grid = space.fbr{obj.dof_c}.^2 .* obj.grid;

% Truncation
obj.grid(obj.grid > hamilt.truncate.delta) = hamilt.truncate.delta;

% Short-time propagator
obj.grid_exp = exp(-1i * obj.grid * time.steps.s_delta * fraction);
