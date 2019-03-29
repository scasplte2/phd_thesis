% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2008 Ulf Lorenz
%
% see the README file for license details.

function obj = kin_jacobi

% Kinetic energy in Jacobi coordinates for the angle.
% See for example J. Chem. Phys 116:4403

% Set the publicly available data
obj.dof_R = [];
obj.dof_r = [];
obj.r_0   = [];
obj.dof_c = [];

obj.mass_R = [];
obj.mass_r = [];

% internal data
obj.grid         = [];
obj.grid_exp     = [];


obj = class(obj, 'kin_jacobi');
