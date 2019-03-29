%-------------------------------------------------
% Plot position density
%-------------------------------------------------
function show_1d_dvr ( obj, state, step )
global expect hamilt info space 


% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2009 Ulf Lorenz
%
% see the README file for license details.

%% Loop over individual wavefunctions
for  m=1:hamilt.coupling.n_eqs

    % If population not too small
    if ( expect.pop.cha{m}(step) > expect.min_pop )
        
        % Get density from wavefunction
        rho = abs   ( state.dvr{m} ) .^2 / obj.scale_dvr;

        % Get phase of wavefunction: Map interval [-pi,pi] into [0,1]
        phi = angle ( state.dvr{m} ) / (2*(pi+0.001)) + 1/2;

        % Get the angle. If we have a grid.legendre, we assume the
        % coordinate is the cosine of the angle, otherwise it is
        % the angle
        if isa(space.dof{1}, 'grid.legendre')
            theta = acos(space.dvr{1});
            % Cover the entire range [0,2pi] instead of [0,pi]
            theta = cat(1, theta, 2*pi - flipud(theta));
            rho   = cat(1, rho, flipud(rho));
            phi   = cat(1, phi, flipud(phi));
        else
            theta = space.dvr{1};
            % Add the first point to the end, so that we have a periodic drawing
            theta = cat(1, theta, theta(1));
            rho   = cat(1, rho, rho(1));
            phi   = cat(1, phi, phi(1));
        end

        % Convert position to Cartesian coordinates; note theta=0 is vertical
        x = rho .* sin(theta);
        y = rho .* cos(theta);

        % Draw colored curve; extra-thick!
        vis.styles.plot_color(x, y, phi, obj.colors(m,:), obj.l_extra, 0, 2)
    end
end

% Draw a few radial lines to guide the eye in steps of 30 degrees
for theta = linspace(0, pi, 7)
    h = line([-sin(theta) sin(theta)], [-cos(theta) cos(theta)] );
    h.LineStyle = obj.patterns{3};
    h.LineWidth = obj.l_thin;
end

% Draw a few circles to guide the eye in steps of 0.25 radius
for r = (0.25:0.25:1)
    angles = linspace(0,2*pi,100);
    h = plot(r*cos(angles), r*sin(angles) );
    h.LineStyle = obj.patterns{3};
    h.LineWidth = obj.l_thin;
end

%% Axes, labels, etc
axis ( [ -1 1 -1 1] )

h = gca; 
h.XAxisLocation = 'top';
h.XTickLabel    = {};
h.YTickLabel    = {};
h.LineWidth     = obj.l_thick;
h.FontName      = obj.f_name;
h.FontSize      = obj.f_large;
h.FontWeight    = obj.f_heavy;

xlabel ( {info.header1;info.header2} )

end
