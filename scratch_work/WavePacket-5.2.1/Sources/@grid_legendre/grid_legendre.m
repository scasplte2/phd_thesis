% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2007 Martin Winter
%               2007-2008 Ulf Lorenz
%
% see the README file for license details.

function obj = grid_legendre

% Just fill some standard values and create the object.

% Required by all grids
obj.label = [];
obj.dof = 1;

% Some values for the kinetic energy operator. Leave the R empty, might lead to
% clashes otherwise.
obj.mass = [];
obj.R_dof = [];
obj.R_0 = [];

% grid options
obj.l_max = [];
obj.m_0 = [];

% internal stuff
obj.kin = [];
obj.kinexpo = [];
obj.nokin = false;

% Arrays for expanding the wave function
obj.trafo_expand = [];
obj.trafo_reconstruct = [];

% And create the object.
obj = class(obj, 'grid_legendre');
