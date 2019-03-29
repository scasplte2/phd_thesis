%------------------------------------------------------------------------------------
% Plot expectation values: sum of energies versus time
%------------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%
% see the README file for license details.

function energies (obj)
global expect hamilt info time

%% Plot energies for total wavefunction (black curves)

% Total energy
h = plot ( time.steps.m_grid(obj.mask_tot), expect.total(obj.mask_tot));
h.LineStyle = obj.patterns{1};
h.LineWidth = obj.l_thick;
h.Color     = 'black';
h.DisplayName = '<E>';
hold on

% Potential energy
if hamilt.coupling.n_eqs==1
    if obj.errorbar
        h = errorbar ( time.steps.m_grid(obj.mask_tot), ...
            expect.pot.tot   (obj.mask_tot), ...
            expect.pot.unc{1}(obj.mask_tot) );
        h.LineStyle = obj.patterns{2};
        h.LineWidth = obj.l_thick;
        h.Color     = 'black';
        h.DisplayName = '<V>';
    else
        h = plot ( time.steps.m_grid(obj.mask_tot), ...
            expect.pot.tot   (obj.mask_tot) );
        h.LineStyle = obj.patterns{2};
        h.LineWidth = obj.l_thick;
        h.Color     = 'black';
        h.DisplayName = '<V>';
    end
else
    h = plot ( time.steps.m_grid(obj.mask_tot), ...
        expect.pot.tot(obj.mask_tot) );
    h.LineStyle = obj.patterns{2};
    h.LineWidth = obj.l_thick;
    h.Color     = 'black';
    h.DisplayName = '<V>';
end

% Kinetic energy
if hamilt.coupling.n_eqs==1
    if obj.errorbar
        h = errorbar ( time.steps.m_grid(obj.mask_tot), ...
            expect.kin.tot   (obj.mask_tot), ...
            expect.kin.unc{1}(obj.mask_tot) );
        h.LineStyle = obj.patterns{3};
        h.LineWidth = obj.l_thick;
        h.Color     = 'black';
        h.DisplayName = '<T>';
    else
        h = plot ( time.steps.m_grid(obj.mask_tot), ...
            expect.kin.tot   (obj.mask_tot) );
        h.LineStyle = obj.patterns{3};
        h.LineWidth = obj.l_thick;
        h.Color     = 'black';
        h.DisplayName = '<T>';
    end
else
    h = plot ( time.steps.m_grid(obj.mask_tot), ...
        expect.kin.tot(obj.mask_tot) );
    h.LineStyle = obj.patterns{3};
    h.LineWidth = obj.l_thick;
    h.Color     = 'black';
    h.DisplayName = '<T>';
end

%% Plot all/potential/kinetic energies for individual wavefunctions (colored curves)
if hamilt.coupling.n_eqs>1
    for m=1:hamilt.coupling.n_eqs
        
        if strcmpi (hamilt.coupling.represent,'dia')
            my_label = hamilt.coupling.labels{m};
        elseif strcmpi (hamilt.coupling.represent,'adi')
            my_label = ['adi-' int2str(m)];
        end
        
        
        % If populations exceed threshold, at least for some time steps
        if ~isempty (obj.mask_cha{m})
            
            % Total energy
            h = plot ( time.steps.m_grid   (obj.mask_cha{m}), ...
                expect.pot.cha{m}(obj.mask_cha{m})+ ...
                expect.kin.cha{m}(obj.mask_cha{m}) );
            h.LineStyle = obj.patterns{1};
            h.LineWidth = obj.l_thick;
            h.Color     = obj.colors(m,:);
            h.DisplayName = my_label;
            
            % Potential, kinetic energy
            if obj.errorbar
                h = errorbar ( time.steps.m_grid   (obj.mask_cha{m}), ...
                    expect.pot.cha{m}(obj.mask_cha{m}), ...
                    expect.pot.unc{m}(obj.mask_cha{m}) );
                h.LineStyle = obj.patterns{2};
                h.LineWidth = obj.l_thick;
                h.Color     = obj.colors(m,:);
                h.DisplayName = my_label;
                h = errorbar ( time.steps.m_grid   (obj.mask_cha{m}), ...
                    expect.kin.cha{m}(obj.mask_cha{m}), ...
                    expect.kin.unc{m}(obj.mask_cha{m}) );
                h.LineStyle = obj.patterns{3};
                h.LineWidth = obj.l_thick;
                h.Color     = obj.colors(m,:);
                h.DisplayName = my_label;
            else
                h = plot ( time.steps.m_grid   (obj.mask_cha{m}), ...
                    expect.pot.cha{m}(obj.mask_cha{m}) );
                h.LineStyle = obj.patterns{2};
                h.LineWidth = obj.l_thick;
                h.Color     = obj.colors(m,:);
                h.DisplayName = my_label;
                h = plot ( time.steps.m_grid   (obj.mask_cha{m}), ...
                    expect.kin.cha{m}(obj.mask_cha{m}) );
                h.LineStyle = obj.patterns{3};
                h.LineWidth = obj.l_thick;
                h.Color     = obj.colors(m,:);
                h.DisplayName = my_label;
            end
        end
    end
end

%% Axes, labels, etc
if ~isempty(obj.e_min)
    e_min = obj.e_min;
elseif isfield (hamilt,'pot_min')
    e_min = hamilt.pot_min;
else
    e_min = -inf;
end

if ~isempty(obj.e_max)
    e_max = obj.e_max;
elseif isfield (hamilt,'pot_max') && isfield (hamilt,'kin_max')
    e_max = hamilt.pot_max + hamilt.kin_max;
else
    e_max = +inf;
end

axis ( [ 0 time.steps.t_total e_min e_max ] )

% Legend explaining the line styles / colors
if obj.legends
    legend('Location','West')  % '<E>','<V>','<T>'
end

% Fonts
h = gca;
h.LineWidth  = obj.l_thick;
h.FontName   = obj.f_name;
h.FontSize   = obj.f_large; 
h.FontWeight = obj.f_heavy;

if isfield(time,'pulse') ||  ~isempty(time.steps.acf)
    h.XTickLabels = []; % Suppress tick labels
else
    switch lower (info.program)
        case {'qm_propa'} % Time dependent simulations
            xlabel ('t')
        case{'qm_bound'} % Time independent: Loop over eigenstates
            xlabel ('n')
        otherwise
            log.error ('Wrong choice of program')
    end
end

if hamilt.coupling.n_eqs==1
    ylabel ( '<V,T,E>' )
else
    if strcmpi (hamilt.coupling.represent,'adi')
        ylabel ( '<V>, <T>, <E>_{adi}' )
    elseif strcmpi (hamilt.coupling.represent,'dia')
        ylabel ( '<V>, <T>, <E>_{dia}' )
    end
end
hold off

end
