%--------------------------------------------------------------------------
% Plot power spectral density (PSD), frequency-resolved optical gating (FROG)
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2015 Burkhard Schmidt
%
% see the README file for license details.

function plot_psd ( iter )
global plots control time

% Initialize figure; draw logos, and resize
if iter==1
    h = figure(9);
    clf;
    show_logo (plots.control);
    set(h,'units','pixels', ...
        'position',[ ...
        plots.control.w_left+plots.control.w_width*2.2 ...
        plots.control.w_lower ...
        plots.control.w_width ...
        plots.control.w_height] );
else
    figure(9)
end

% Length of time series; enforce even number of time steps
N = control.t.n;
if mod(N,2) 
    N = N-1;
end

% Maximal frequency
fmax = control.f.forward(round((N/2)/time.frog.zoom));

%% Curve plot of PSD (left panel)
if ~strcmpi(time.frog.choice,'none')
    subplot( 'Position', [1 1 3 6]/8 )
end
plot(...
    control.p.forward(1:round((N/2)/time.frog.zoom)), ...
    control.f.forward(1:round((N/2)/time.frog.zoom)), ...
    'LineStyle',    '-', ...
    'LineWidth',    plots.control.l_thin)

axis ([0 inf 0 fmax])
title (control.title);
xlabel('power')
ylabel('frequency')

set ( gca, ...
    'LineWidth',    plots.control.l_thick, ...
    'FontName',     plots.control.f_name,  ...
    'FontSize',     plots.control.f_large, ...
    'FontWeight',   plots.control.f_heavy )

%% Contour plot of FROG (right panel)
if ~strcmpi(time.frog.choice,'none')
    
    subplot( 'Position', [4 1 3 6]/8 )
    contour(...
        control.t.steps  (1:N), ...
        control.f.forward(1:round((N/2+1)/time.frog.zoom)),...
        control.g.forward(1:round((N/2+1)/time.frog.zoom),1:N),...
        50)
    
    axis ([0 time.steps.t_total 0 fmax])
    if strcmpi(time.frog.choice,'gauss')
        title ({'FROG',time.frog.choice,['width = ' num2str(time.frog.width)]});
    else
        title ({'FROG',time.frog.choice});
    end
    xlabel('delay time')
    ylabel('frequency')
    
    set ( gca, ...
        'LineWidth',     plots.control.l_thick, ...
        'FontName',      plots.control.f_name,  ...
        'FontSize',      plots.control.f_large, ...
        'FontWeight',    plots.control.f_heavy, ...
        'YAxisLocation', 'right')
end

drawnow

saveas (9,'psd.fig')
saveas (9,'psd.jpg')
