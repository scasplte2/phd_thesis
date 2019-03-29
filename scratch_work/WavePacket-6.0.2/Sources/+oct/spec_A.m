%--------------------------------------------------------------------------
% Calculate spectrum of A-matrix and
% plot it in the complex number plane
% -------------------------------------------------------------------------

function spec_A (choice)

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2012 Boris Schaefer-Bung, Burkhard Schmidt, Ulf Lorenz
%
% see the README file for license details.

global balanced bilinear 

%% Distinguish different choices
switch lower(choice)
    case 'qm_abncd'
        A = bilinear.A;
        T = {'Matrix A: spectrum';[bilinear.title int2str(size(A,1)) ' coupled ODEs']};
        figure(91);
    case 'qm_balance'
        A = balanced.A;
        T = {'Matrix A: spectrum after balancing';[balanced.title int2str(size(A,1)) ' coupled ODEs']};
        figure(92);
    case 'qm_truncate'
        A = bilinear.A;
        T = {'Matrix A: spectrum after reduction';[bilinear.title int2str(size(A,1)) ' coupled ODEs']};
        figure(93);
    case 'qm_h2model'
        A = bilinear.A;
        T = {'Matrix A: spectrum after reduction';[bilinear.title int2str(size(A,1)) ' coupled ODEs']};
        figure(94);
    otherwise
        log.error (['Invalid choice for plotting spectrum of A matrix: ' choice])
end

clf;
thisplot = vis.styles; % Construct object
show_logo (thisplot)

%% Get spectrum and plot it in the complex number plane
spec = eig(full(A));

plot(real(spec),imag(spec),'o', ...
    'MarkerSize'     ,thisplot.m_large, ...
    'MarkerEdgeColor',thisplot.colors(2,:),...
    'MarkerFaceColor',thisplot.colors(2,:));

set ( gca, 'LineWidth',     thisplot.l_thick, ...
    'FontName',      thisplot.f_name,  ...
    'FontSize',      thisplot.f_large, ...
    'FontWeight',    thisplot.f_heavy )
xlabel('Real part')
ylabel('Imaginary part')
title(T)
