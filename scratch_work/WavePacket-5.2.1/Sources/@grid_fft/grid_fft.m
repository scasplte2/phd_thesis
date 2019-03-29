% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%
% see the README file for license details.

function obj = grid_fft

%% fill some standard values and just return the object

obj.label = [];
obj.dof = 1;
obj.mass = 1;
obj.periodic = true;

% Grid in possition space
obj.n_pts = [];
obj.x_min = [];
obj.x_max = [];

% Grid in momentum space
obj.p_min = [];
obj.p_max = [];

% Grid representation of kinetic operator
obj.kin = [];
obj.nokin = false;

% internal grid representation of the kinetic energy operator/propagator
obj.intern_kin = [];
obj.intern_kinexpo = [];
obj.intern_factor = [];

% Create object
obj = class(obj, 'grid_fft');
