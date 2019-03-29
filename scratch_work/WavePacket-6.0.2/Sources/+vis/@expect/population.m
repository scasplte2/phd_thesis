%--------------------------------------------------------------------
% Plot populations and optional AMOs vs. time
%--------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%
% see the README file for license details.

function population (obj)
global expect hamilt info time

% Plot total population (black curves, solid)
h = plot ( time.steps.m_grid (obj.mask_tot), expect.pop.tot (obj.mask_tot));
h.LineStyle = obj.patterns{1};
h.LineWidth = obj.l_thick;
h.Color     = 'black';
h.DisplayName = 'Norm'; % 'Population';
hold on

% Plot additional multiplicative operators (black curves, different linestyles)
if isfield(hamilt, 'amo')
    for p = 1:length(hamilt.amo)
        if ~isempty (hamilt.amo{p})
            h = plot ( time.steps.m_grid (obj.mask_tot), expect.amo{p}.tot (obj.mask_tot));
            h.LineStyle = obj.patterns{p+1};
            h.LineWidth = obj.l_thick;
            h.Color     = 'black';
            h.DisplayName = hamilt.amo{p}.label;
        end
    end
end


% Plot populations for each channel (individual colors)
if hamilt.coupling.n_eqs > 1
    for m=1:hamilt.coupling.n_eqs
        
        % If populations exceed threshold, at least for some time steps
        if ~isempty (obj.mask_cha{m})
            
            % Plot populations  (colored curves, solid)
            h = plot ( time.steps.m_grid (obj.mask_cha{m}), expect.pop.cha{m} (obj.mask_cha{m}));
            h.LineStyle = obj.patterns{1};
            h.LineWidth = obj.l_thick;
            h.Color     = obj.colors(m,:);
            
            % Legend explaining the colors
            if strcmpi (hamilt.coupling.represent,'dia')
                h.DisplayName = hamilt.coupling.labels{m};
            elseif strcmpi (hamilt.coupling.represent,'adi')
                h.DisplayName = ['adi-' int2str(m)];
            end
            
            % Plot additional multiplicative operators
            if isfield(hamilt, 'amo')
                for p = 1:length(hamilt.amo)
                    if ~isempty (hamilt.amo{p})
                        h = plot ( time.steps.m_grid (obj.mask_cha{m}), expect.amo{p}.cha{m} (obj.mask_cha{m}));
                        h.LineStyle = obj.patterns{p+1};
                        h.LineWidth = obj.l_thick;
                        h.Color     = obj.colors(m,:);
                        
                        % Legend explaining the linestyles
                        h.DisplayName = hamilt.amo{p}.label;
                    end
                end
            end
        end
    end
end

% Draw legend
if obj.legends
    legend('Location','SouthWest')
end

% Axes and labels
axis ( [ 0 time.steps.t_total obj.p_min obj.p_max ] );
h = gca;
h.LineWidth   = obj.l_thick;
h.FontName    = obj.f_name;
h.FontSize    = obj.f_large; 
h.FontWeight  = obj.f_heavy;

h.XAxisLocation = 'top';
switch lower (info.program)
    case {'qm_propa'} % Time dependent simulations
        xlabel ('t')
    case{'qm_bound'} % Time independent: Loop over eigenstates
        xlabel ('n')
    otherwise
        log.error ('Wrong choice of program')     
end

if hamilt.coupling.n_eqs==1
    ylabel ( 'N(t)' )
else
    if strcmpi (hamilt.coupling.represent,'adi')
        ylabel ( 'N_{adi}(t)' )
    elseif strcmpi (hamilt.coupling.represent,'dia')
        ylabel ( 'N_{dia}(t)' )
    end
end

hold off

end
