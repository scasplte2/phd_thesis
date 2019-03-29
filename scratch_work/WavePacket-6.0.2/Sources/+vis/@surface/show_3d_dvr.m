%--------------------------------------------------------------
%
% Visualize wavepacket in 3 dimensions in position (DVR) space 
% using Matlab's builtin isosurface plots
%
%--------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%
% see the README file for license details.

function show_3d_dvr ( obj, state, step )
global hamilt info space

if hamilt.coupling.n_eqs~= 1
    log.error ('Can not yet draw a 3D isosurface for coupled TDSEs')
end

% Get position density
rho = abs ( state.dvr{1} ) .^2;

% Deletes all visible graphics objects from previous time steps
if step>1
    cla 
end

% Surface plot of densities
p = patch(isosurface ( space.dvr{1},space.dvr{2},space.dvr{3}, ...
    rho, obj.scale_dvr/10));
p.FaceColor = 'green';
p.EdgeColor = 'none';

% Specify view point in terms of azimuth and elevation
view (obj.srf_view(1),obj.srf_view(2))

% Lighting of surfaces
camlight;  
lighting PHONG;
if obj.srf_look(2)
    lightangle(obj.srf_light(1),obj.srf_light(2));
end

% Axes, labels, etc
if ~obj.range
    axis ([ ...
        space.dof{1}.dvr_min space.dof{1}.dvr_max ...
        space.dof{2}.dvr_min space.dof{2}.dvr_max ...
        space.dof{3}.dvr_min space.dof{3}.dvr_max ]);
else
    axis ([ ...
        obj.x_min obj.x_max ...
        obj.y_min obj.y_max ...
        obj.z_min obj.z_max ]);
end

h = gca;
h.LineWidth  = obj.l_thick;
h.FontName   = obj.f_name;
h.FontSize   = obj.f_large;
h.FontWeight = obj.f_heavy;

title  ( {info.header1;info.header2} )
xlabel ( [ 'R_{', space.dof{1}.label, '}' ] )
ylabel ( [ 'R_{', space.dof{2}.label, '}' ] )
zlabel ( [ 'R_{', space.dof{3}.label, '}' ] )

end