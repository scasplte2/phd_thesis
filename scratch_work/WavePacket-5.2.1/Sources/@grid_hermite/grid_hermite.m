% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2007-2008 Ulf Lorenz
%
% see the README file for license details.

function obj = grid_hermite

% set some data that all grid have to supply
obj.label = [];
obj.dof = 1;

% Some values for the kinetic energy operator and the scaling
obj.mass  = 1;
obj.omega = 1;
obj.v_2   = [];

% grid options
obj.n_pts = [];
obj.r_e   = 0;

% data related to the internal energy calculations
obj.momentum = [];
obj.kin = [];
obj.kinexpo  = [];
obj.kin_max = 0;
obj.nokin = false;

% arrays for transforming the wave function
obj.trafo_expand = [];
obj.trafo_reconstruct = [];

% finally, create the object
obj = class(obj, 'grid_hermite');
