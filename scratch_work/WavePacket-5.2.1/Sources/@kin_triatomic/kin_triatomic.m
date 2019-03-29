% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2007-2008 Ulf Lorenz
%
% see the README file for license details.

function obj = kin_triatomic

% Set the publicly available data
obj.dof = [];
obj.mass = [];
obj.theta = [];


% Set the private data
obj.grid = [];
obj.grid_exp = [];

obj = class(obj, 'kin_triatomic');
