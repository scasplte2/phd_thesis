%--------------------------------------------------------------------------
%
% Visualize only the density either in DVR (coordinate space) or FBR
% ("momentum" space). This plot accepts essentially all kinds of grids, but
% does not look good with all of them.
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2009 Ulf Lorenz
%               2009 Burkhard Schmidt
%
% see the README file for license details.

function curve ( step, wide )
global plots space

if wide % Wide format: 16:9
    w=16; h=09;
else % Square format: 9:9
    w=09; h=09;
end

if space.size.n_dim ~= 1
    util.error('Cannot draw a curve plot for >1 dimensions. Use the contour plot instead')
end

subplot ( 'Position', [1/w 3/h 7/w 5/h] );
hold off; plot( 1234567890, 1234567890 ); hold on;

if strcmp(plots.density.representation, 'dvr')
    density_dvr ( step );
else
    density_fbr ( step );
end

subplot ( 'Position', [1/w 1/h 7/w 1/h] )
my_colors

%-------------------------------------------------
% Plot position density and potential energy curve
%-------------------------------------------------
function density_dvr ( step )
global expect hamilt info plots psi space 


%% Loop over individual wavefunctions
for  m=1:hamilt.coupling.n_eqs;

    % If population not too small
    if ( expect.ind.pop{m}(step) > expect.min_pop )
        
        % Get density from wavefunction
        rho = abs   ( psi.dvr.grid_ND{m} ) .^2;

        % Get phase of wavefunction: Map interval [-pi,pi] into [0,1]
        phi = angle ( psi.dvr.grid_ND{m} ) / (2*(pi+0.001)) + 1/2;

        % Plot density and phase (with horizontal offset and base line)
        if plots.density.energy.on
            offset = expect.ind.pot{m}(step);
        else
            offset = 0;
        end
        
        plot.color ( space.dvr.grid_1D{1}, ...
                     rho*plots.pot.delta/plots.dvr.rho_max, ...
                     phi, ...
                     plots.style.colors(m,:), ...
                     plots.style.line.extrathick, ...
                     offset, ...
                     0 )
       if plots.density.energy.on
            line ( [space.dof{1}.dvr_min space.dof{1}.dvr_max], ...
                   [offset offset], ...
                   'LineStyle', '-', ...
                   'Color',     plots.style.colors(m,:), ...
                   'LineWidth', plots.style.line.thin) 
        end
        
    end

    % Plot potential energy curve
    if plots.density.energy.on
        plot ( space.dvr.grid_1D{1}, ...
               hamilt.pot.grid_ND{m,m}, ...
               'LineStyle', '-', ...
               'Color', plots.style.colors(m,:), ...
               'LineWidth', plots.style.line.thin )
    end
    
end

%% Axes, labels, etc
if plots.density.range.on == false
    axis ( [ space.dof{1}.dvr_min space.dof{1}.dvr_max ...
            plots.pot.min-plots.pot.delta/10 plots.pot.max+plots.pot.delta/10 ] )
else
    axis ( [ plots.density.range.x_min plots.density.range.x_max ...
            plots.pot.min-plots.pot.delta/10 plots.pot.max+plots.pot.delta/10 ] )
end

set ( gca, 'XAxisLocation', 'bottom', ...
           'LineWidth',     plots.style.line.thick, ...
           'FontName',      plots.style.font.name,  ...
           'FontSize',      plots.style.font.large, ...
           'FontWeight',    plots.style.font.heavy)
xlabel('R [a_0]')
title ( {info.header1;info.header2} )
if plots.density.energy.on

    if hamilt.coupling.n_eqs==1
        ylabel ( 'V(R) [E_h]' )
    else
        if strcmpi ( hamilt.coupling.representation,'adi')
            ylabel ( 'V_{adi}(R) [E_h]' )
        elseif strcmpi ( hamilt.coupling.representation,'dia')
            ylabel ( 'V_{dia}(R) [E_h]' )
        end
    end
    
else
    
    if hamilt.coupling.n_eqs==1
        ylabel ( '\rho(R)' )
    else
        if strcmpi ( hamilt.coupling.representation,'adi')
            ylabel ( '\rho_{adi}(R)' )
        elseif strcmpi ( hamilt.coupling.representation,'dia')
            ylabel ( '\rho_{dia}(R)' )
        end
    end

end

% Negative imaginary potential (as absorbing boundary condition)
if ~isempty(hamilt.nip.grid)
    if hamilt.nip.params.min(1) > space.dof{1}.dvr_min
        line ( [hamilt.nip.params.min(1) hamilt.nip.params.min(1)], ...
               [plots.pot.min-plots.pot.delta/10 plots.pot.max+plots.pot.delta/10], ...
               'LineStyle', '--', ...
               'Color',     'k', ...
               'LineWidth', plots.style.line.thin)
    end
    
    if hamilt.nip.params.max(1) < space.dof{1}.dvr_max
        line ( [hamilt.nip.params.max(1) hamilt.nip.params.max(1)], ...
               [plots.pot.min-plots.pot.delta/10 plots.pot.max+plots.pot.delta/10], ...
               'LineStyle', '--', ...
               'Color',     'k', ...
               'LineWidth', plots.style.line.thin)
    end
end

%----------------------------------------------------------
% Plot momentum density and kinetic energy curve
%----------------------------------------------------------
function density_fbr ( step ) 
global expect hamilt plots psi space 

%% Loop over individual wavefunctions
for  m=1:hamilt.coupling.n_eqs

    % If population not too small
    if ( expect.ind.pop{m}(step) > expect.min_pop )

        psi.fbr.grid_ND{m} = dvr2fbr(space.dof{1}, psi.dvr.grid_ND{m});
        
        % Get density from wavefunction
        rho = abs ( psi.fbr.grid_ND{m} ) .^2;

        % Get phase of wavefunction: Map interval [-pi,pi] into [0,1]
        phi = angle ( psi.fbr.grid_ND{m} ) / (2*(pi+0.001)) + 1/2;

        % Plot density and phase (with horizontal offset and base line)
        if plots.density.energy.on
            offset = expect.ind.kin{m}(step);
        else
            offset = 0;
        end
            
        plot.color ( space.fbr.grid_1D{1}, ...
                      rho*plots.kin.delta/plots.fbr.rho_max, ...
                      phi, ...
                      plots.style.colors(m,:), ...
                      plots.style.line.extrathick, ...
                      offset, ...
                      0 )
                  
        if plots.density.energy.on
               line ( [space.dof{1}.fbr_min space.dof{1}.fbr_max], ...
                      [ offset        offset       ], ...
                      'LineStyle', '-', ...
                      'Color',     plots.style.colors(m,:), ...
                      'LineWidth', plots.style.line.thin)
        end
                      
    end

end

%% Plot kinetic energy curve
if plots.density.energy.on
    plot ( -space.fbr.grid_1D{1}, ...
            space.dof{1}.kin, ...
            'Color',    'k', ...
            'LineWidth', plots.style.line.thin )
end
    
%% Axes, labels, etc
if plots.density.range.on == false
    axis ( [ space.dof{1}.fbr_min space.dof{1}.fbr_max ...            
             plots.kin.min-plots.kin.delta/10 plots.kin.max+plots.kin.delta/10 ] )
else
    axis ( [ plots.density.range.y_min plots.density.range.y_max ...
            plots.kin.min-plots.kin.delta/10 plots.kin.max+plots.kin.delta/10 ] )
end

set ( gca, 'XAxisLocation', 'top', ...
           'YAxisLocation', 'left', ...
           'LineWidth',     plots.style.line.thick, ...
           'FontName',      plots.style.font.name,  ...
           'FontSize',      plots.style.font.large, ...
           'FontWeight',    plots.style.font.heavy )

if plots.density.energy.on
        ylabel ( 'T(P) [E_h]' )
else
    
    if hamilt.coupling.n_eqs==1
        ylabel ( '\rho(P)' )
    else
        if strcmpi ( hamilt.coupling.representation,'adi' )
            ylabel ( '\rho_{adi}(P)' )
        elseif strcmpi ( hamilt.coupling.representation,'dia' )
            ylabel ( '\rho_{dia}(P)' )
        end
    end

end
xlabel ('P [hbar / a_0]')
title ( {info.header1;info.header2} )

%----------------------------------------------------------
% Plot color mapping
%----------------------------------------------------------
function my_colors 
global plots

x = linspace ( -1, 1, 50 );
y = ones (size(x));
c = linspace ( 0, 1, 50 );
plot.color ( x, y, c, [1 0 0], plots.style.line.extrathick, 0, 0) 
axis ( [-1 1 0.9 1.1] )
set ( gca, 'XTick',          -0.5:0.5:1, ...
           'YTick',         [-123456789 123456789], ... % Fake: suppress tick labels
           'YAxisLocation', 'left', ...
           'LineWidth',     plots.style.line.thick, ...
           'FontName',      plots.style.font.name,  ...
           'FontSize',      plots.style.font.large, ...
           'FontWeight',    plots.style.font.heavy )
xlabel ('phase [\pi]')
ylabel (['color';' map '])

