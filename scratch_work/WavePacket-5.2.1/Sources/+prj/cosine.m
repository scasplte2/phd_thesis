%--------------------------------------------------------------------------
%
% Degree of alignment (cos^n) as projection function (PRJ)
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2007-2008 Ulf Lorenz
%               2008 Burkhard Schmidt
%
% see the README file for license details.

function cosine

global space

%% set default values
if ~isfield(space.prj, 'params')
    space.prj.params.dof = [];
end

if ~isfield(space.prj.params, 'dof') || isempty(space.prj.params.dof)
    space.prj.params.dof = 1;
end
if ~isfield(space.prj.params, 'exp') || isempty(space.prj.params.exp)
    space.prj.params.exp = 1;
end

%% Output

util.disp (' ')
util.disp ( '****************************************************' )
util.disp ( 'Directional property (cos^n) as projection function ' )
util.disp ( '                                                    ' )
util.disp ( '           N      n                                 ' )
util.disp ( 'prj(R) =  Sum  cos  Theta_i                         ' )
util.disp ( '          i=1                                       ' )
util.disp ( '                                                    ' )
util.disp ( 'where the sum extends over all grid points in the   ' )
util.disp ( 'selected coordinate.                                ' )
util.disp ( '****************************************************' )

util.disp( ['Exponent n: ' num2str(space.prj.params.exp)])
if space.size.n_dim > 1
    util.disp ( ['Degree of freedom: ' space.prj.params.dof] )
end
if strcmp(class(space.dof{space.prj.params.dof}), 'grid_legendre')
    util.disp( 'Treating position variable x as cos \theta' )
end

%% Create the projection grid
if strcmp(class(space.dof{space.prj.params.dof}), 'grid_legendre')
    space.prj.grid_ND = space.dvr.grid_ND{space.prj.params.dof}.^space.prj.params.exp;
else
    space.prj.grid_ND = cos(space.dvr.grid_ND{space.prj.params.dof}).^space.prj.params.exp;
end
