%--------------------------------------------------------------------
% Plot autocorrelation vs. time
%--------------------------------------------------------------------
function correlation (obj)
global time

% Not too many time steps: plot colored curve
if time.steps.s_number * time.steps.m_number<=1000
    
    % Get modulus of autocorrelation function
    rho = abs   ( time.steps.acf ) .^2;
    
    % Get phase of autocorrelation function: Map interval [-pi,pi] into [0,1]
    phi = angle ( time.steps.acf ) / (2*(pi+0.001)) + 1/2;
    
    % Plot autocorrelation function (colored curve)
    vis.styles.plot_color ( time.steps.s_grid, rho, phi, [1 1 0]/2, obj.l_thick, 0, 0 )
    
    % Labelling the ordinate
    my_label = 'C(t)';
    
    % Too many time steps: plot black curve showing abs(C) instead
else

    h = plot ( time.steps.s_grid, abs(time.steps.acf) );
    h.LineStyle = obj.patterns{1};
    h.LineWidth = obj.l_thick;
    h.Color     = 'black';
    
    % Labelling the ordinate
    my_label = '|C(t)|';
    
end

% Axes and labels
axis ( [ 0 time.steps.t_total -0.1 1.1 ] )
h = gca;
h.LineWidth  = obj.l_thick;
h.FontName   = obj.f_name;
h.FontSize   = obj.f_large; 
h.FontWeight = obj.f_heavy;

xlabel ('t')
ylabel (my_label)

end