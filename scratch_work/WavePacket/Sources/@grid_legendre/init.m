% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2007 Martin Winter
%               2007-2008 Ulf Lorenz
%
% see the README file for license details.

function [ obj, weights, x_grid, p_grid ] = init(obj)

% The stuff itself is calculated by Gaussian quadrature in an internal function.
[weights, x_grid, obj.trafo_expand] = quadrature(obj.l_max - abs(obj.m_0) + 1, abs(obj.m_0));

% The spectral grid is sqrt( l(l+1) )
p_grid = (abs(obj.m_0):obj.l_max)';
p_grid = sqrt(p_grid.*(p_grid+1));

% The inverse trafo is just the transpose of the simple trafo, since it is
% unitary (even orthogonal)
obj.trafo_reconstruct = obj.trafo_expand';

% Ah, yes. Don't forget to add the weights to the transformation matrix, since it
% involves an integration/summation over the DVR grid.
obj.trafo_expand = obj.trafo_expand .* repmat(weights', [ length(weights) 1] );

%% Informational output
util.disp ( ' ' )
util.disp ( '**************************************************' )
util.disp ( [ 'DVR for the degree of freedom: ' obj.label] )
util.disp ( 'Discretization scheme: Gauss-Legendre DVR' )
util.disp ( '**************************************************' )
util.disp ( ' ' )
util.disp ( 'Parameters of the (angular) grid' )
util.disp ( [ 'Number of grid points       : ' int2str(obj.l_max - abs(obj.m_0) + 1) ] )
util.disp ( [ 'Lower bound for L           : ' int2str(abs(obj.m_0)) ] )
util.disp ( [ 'Upper bound for L           : ' int2str(obj.l_max) ] )
util.disp ( [ 'Azimuthal quantum number m  : ' int2str(obj.m_0) ] )
