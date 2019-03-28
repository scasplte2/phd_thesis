% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%
% see the README file for license details.

function [ obj, weight, x_grid, p_grid ] = init(obj)

% formerly known as space.pos.delta, this is the (constant) weight
% function
weight = ones(obj.n_pts, 1) * (obj.x_max - obj.x_min) / obj.n_pts;

% the grid points in coordinate space as a column (!) vector
x_grid = linspace(obj.x_min, obj.x_max - weight(1), obj.n_pts)';

% The mapping from grid points in the DVR pseudospectral basis (coordinate space)
% to the FBR spectral basis (momentum space)

obj.p_min = - pi / weight(1);
obj.p_max = -obj.p_min;

p_grid = linspace(obj.p_min, obj.p_max - 2*obj.p_max/obj.n_pts, obj.n_pts)';


%% Informational output
util.disp ( ' ' )
util.disp ( '**************************************************' )
util.disp ( [ 'DVR for the degree of freedom: ' obj.label] )
util.disp ( 'Discretization scheme: Equally spaced grids (FFT)' )
util.disp ( '**************************************************' )
util.disp ( [ 'Number of grid  points : ' num2str(obj.n_pts)    ] )
util.disp ( ' ' )
util.disp ( 'Position space' )
util.disp ( [ 'Minimum of grid  : ' num2str(obj.x_min)    ] )
util.disp ( [ 'Maximum of grid  : ' num2str(obj.x_max)    ] )
util.disp ( [ 'Grid spacing     : ' num2str(weight(1))  ] )
util.disp ( ' ' )
util.disp ( 'Momentum (wavenumber) space (Fourier transform)' )
util.disp ( [ 'Maximum of grid  : ' num2str(obj.p_max)    ] )
util.disp ( [ 'Grid spacing     : ' num2str(2*obj.p_max/obj.n_pts)  ] )
util.disp ( ' ' )
util.disp ( 'Momentum (wavenumber) space (Wigner transform)' )
util.disp ( [ 'Maximum of grid  : ' num2str(obj.p_max/2)  ] )
util.disp ( [ 'Grid spacing     : ' num2str(obj.p_max/obj.n_pts)] )
