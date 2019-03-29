%--------------------------------------------------------------------------
%
% Plot evolution of input/state/output of bilinear control problem 
% in one figure with 2 or 3 subplots 
% (1) input=control field(s) u(t) - not drawn for the field-free case - 
% (2) state vector x(t)
% (3) output=observable(s) y(t) 
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%
% see the README file for license details.

function  initial ( obj )

global bilinear control

% Get the number of subplots: 2 (without control) or 3 (with control)
obj.numPlots = 2 + isfield(control,'u');

% Initialize figure and resize
h = figure(7);
clf;
set(h,'units','pixels', ...
    'position',[ ...
    obj.w_left ...
    obj.w_lower ...
    obj.w_width ...
    obj.w_height] );

% Optionally open a movie file
if control.plot.mov
    log.disp ('Opening animated uxy plot file : uxy.mp4')
    obj.writerObj = VideoWriter ('uxy', 'MPEG-4');
    open(obj.writerObj);
end

% Prepare legends for u
if isfield(control,'u')
    control.u.legends = cell(control.u.dim,1);
    for i=1:control.u.dim
        control.u.legends{i}=['u_{' int2str(i) '}(t)'];
    end
end

% Prepare legends for x
control.x.max = 7;
control.x.legends = cell(min(control.x.dim,control.x.max),1);
for i=1:min(control.x.dim,control.x.max)
    control.x.legends{i}=['x_{' int2str(i) '}(t)'];
end

% Prepare legends for y
control.y.legends = bilinear.label;


end

