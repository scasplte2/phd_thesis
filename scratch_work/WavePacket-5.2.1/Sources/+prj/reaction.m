% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2008 Ulf Lorenz
%               2008 Burkhard Schmidt
%
% see the README file for license details.


function reaction

global space

%% Set default values and do some checks
if space.size.n_dim < 2
    util.error('prj.react requires at least two degrees of freedom!');
end


if ~isfield(space.prj, 'params') || ~isfield(space.prj.params, 'reac')
    space.prj.params.reac = 1;
end

if ~isfield(space.prj.params, 'prod')
    space.prj.params.prod = 2;
end

if ~isfield(space.prj.params, 'side')
    space.prj.params.side = 'p';
end

%% Output and more checks
util.disp(' ')
util.disp('*************************************')
util.disp('Projection on educt/product region   ')
util.disp('for a chemical exchange reaction     ')
util.disp('          A + BC -> AB + C           ')
util.disp('*************************************')
util.disp(' ')
util.disp(['Index of educt distance coordinate AB :' num2str(space.prj.params.reac)])
util.disp(['Index of product distance coordinate BC :' num2str(space.prj.params.prod)])
if space.prj.params.side == 'r'
    util.disp('Projecting on reactant part')
elseif space.prj.params.side == 'p'
    util.disp('Projecting on product part')
else
    util.error('space.prj.params.side has to be "r" or "p"')
end

%% Set the grid
if space.prj.params.side == 'r'
    space.prj.grid_ND = ...
        (space.dvr.grid_ND{space.prj.params.reac} > space.dvr.grid_ND{space.prj.params.prod});
else
    space.prj.grid_ND = ...
        (space.dvr.grid_ND{space.prj.params.reac} < space.dvr.grid_ND{space.prj.params.prod});
end
