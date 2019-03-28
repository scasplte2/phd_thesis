%--------------------------------------------------------------------------
%
% Visualize reduced densities from multidimensional wavepacket propagation
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2007-2010 Ulf Lorenz
%
% see the README file for license details.

function reduced ( step, wide )
global space psi hamilt

if wide % Wide format: 16:9
    w=16; h=09;
else % Square format: 9:9
    w=09; h=09;
end

% Since the data is displayed using a Wigner transform, this for now only works
% for the simple case of FFT grids until the next rewrite
for k = 1:space.size.n_dim
    if ~strcmp(class(space.dof{k}), 'grid_fft')
        util.error('Reduced density plot only available for FFT grids')
    end
end

%% Calculate reduced densities

% Pre-allocate cell arrays. psi.dvr.redu{m,k} is the reduced density on
% surface m reduced to the degree of freedom k.
psi.dvr.redu = cell (hamilt.coupling.n_eqs, space.size.n_dim);

% Loop over (coupled) wavefunctions
for m = 1:hamilt.coupling.n_eqs
    
    for k = 1:space.size.n_dim
        psi.dvr.redu{m,k} = zeros(space.dof{k}.n_pts);
    end
    
    switch space.size.n_dim
        case 1
            rho_redu_1d(m);
        case 2
            rho_redu_2d(m);
        case 3
            rho_redu_3d(m);
        case 4
            rho_redu_4d(m);
        case 5
            rho_redu_5d(m);
        otherwise
            util.error ('Reduced density matrix only up to 5 dimensions')
    end
   
    % Weights are very simple for FFT
	% Formally not really correct: As Max B. keeps telling us,
	% the reduced densities should be normalized with the 
	% weights in N-1 dimensions. However, this error is comp-
	% ensated in the calculation of the purities!
    for k = 1:space.size.n_dim
        psi.dvr.redu{m,k} = psi.dvr.redu{m,k} * space.dvr.weight_ND(1);
    end

end

%% Visualize reduced densities
switch space.size.n_dim
    
    case 1 % note that tr(rho^2)=1 in one dimension
        
        % Horizontal = Vertical = 2:5:2
        subplot ( 'Position', [1/w 1/h 7/w 7/h] ); redu_dens ( step, 1 );
        
        
    case 2 % note that tr(rho_1^2)=tr(rho_2^2) in two dimensions

        % Horizontal = vertical = 2:5:2:5:2
        subplot ( 'Position', [1/w 5/h 3/w 3/h] ); redu_dens ( step, 1 );
        subplot ( 'Position', [5/w 5/h 3/w 3/h] ); redu_dens ( step, 2 );
        
        subplot ( 'Position', [1/w 1/h 7/w 3/h] ); purity ( step, 1 );
        
    otherwise

        % Horizontal = 2:5:2:5:2:5:2..., Vertical = 2:5:2:5:2
        for k=1:space.size.n_dim
            
            n = space.size.n_dim;
            
            subplot ( 'Position', [ (2+7*(k-1))/(2+7*n)*9/w 5/h 5/(2+7*n)*9/w 3/h] );
            redu_dens ( step, k );
            
            subplot ( 'Position', [ (2+7*(k-1))/(2+7*n)*9/w 1/h 5/(2+7*n)*9/w 3/h] );
            purity ( step, k );
        end
end
        
%------------------------------------------------------------------------
% Obtain density matrix from 1-dimensional wavepacket propagation
%-------------------------------------------------------------------------
function rho_redu_1d (m)
global psi space
    
for ii=1:space.dof{1}.n_pts
    for jj=1:space.dof{1}.n_pts
        psi.dvr.redu{m,1}(ii,jj) = conj(psi.dvr.grid_ND{m}(ii)) .* psi.dvr.grid_ND{m}(jj);
    end
end

%------------------------------------------------------------------------
% Reduce total density matrix (from 2-dimensional wavepacket propagation)
% to k-th dimension (k=1,2) by tracing over all the remaining dimensions
%-------------------------------------------------------------------------
function rho_redu_2d (m)   
global psi space
    
% First coordinate
for ii=1:space.dof{1}.n_pts
    for jj=1:space.dof{1}.n_pts
        psi.dvr.redu{m,1}(ii,jj) = sum ( conj(psi.dvr.grid_ND{m}(ii,:)) .* psi.dvr.grid_ND{m}(jj,:) );
    end
end

% Second coordinate
for ii=1:space.dof{2}.n_pts
    for jj=1:space.dof{2}.n_pts
        psi.dvr.redu{m,2}(ii,jj) = sum ( conj(psi.dvr.grid_ND{m}(:,ii)) .* psi.dvr.grid_ND{m}(:,jj) );
    end
end

%------------------------------------------------------------------------
% Reduce total density matrix (from 3-dimensional wavepacket propagation)
% to k-th dimension (k=1,2,3) by tracing over all the remaining dimensions
%-------------------------------------------------------------------------
function rho_redu_3d (m)
global psi space

% First coordinate
for ii=1:space.dof{1}.n_pts
    for jj=1:space.dof{1}.n_pts
        psi.dvr.redu{m,1}(ii,jj) = sum ( sum ( conj(psi.dvr.grid_ND{m}(ii,:,:)) .* psi.dvr.grid_ND{m}(jj,:,:) ) );
    end
end

% Second coordinate
for ii=1:space.dof{2}.n_pts
    for jj=1:space.dof{2}.n_pts
        psi.dvr.redu{m,2}(ii,jj) = sum ( sum ( conj(psi.dvr.grid_ND{m}(:,ii,:)) .* psi.dvr.grid_ND{m}(:,jj,:) ) );
    end
end

% Third coordinate
for ii=1:space.dof{3}.n_pts
    for jj=1:space.dof{3}.n_pts
        psi.dvr.redu{m,3}(ii,jj) = sum ( sum ( conj(psi.dvr.grid_ND{m}(:,:,ii)) .* psi.dvr.grid_ND{m}(:,:,jj) ) );
    end
end

%------------------------------------------------------------------------
% Reduce total density matrix (from 4-dimensional wavepacket propagation)
% to k-th dimension (k=1,2,3,4) by tracing over all the remaining dimensions
%-------------------------------------------------------------------------
function rho_redu_4d (m)
global psi space

% First coordinate
for ii=1:space.dof{1}.n_pts
    for jj=1:space.dof{1}.n_pts
        psi.dvr.redu{m,1}(ii,jj) = sum ( sum ( sum ( conj(psi.dvr.grid_ND{m}(ii,:,:,:)) .* psi.dvr.grid_ND{m}(jj,:,:,:) ) ) );
    end
end

% Second coordinate
for ii=1:space.dof{2}.n_pts
    for jj=1:space.dof{2}.n_pts
        psi.dvr.redu{m,2}(ii,jj) = sum ( sum ( sum ( conj(psi.dvr.grid_ND{m}(:,ii,:,:)) .* psi.dvr.grid_ND{m}(:,jj,:,:) ) ) );
    end
end

% Third coordinate
for ii=1:space.dof{3}.n_pts
    for jj=1:space.dof{3}.n_pts
        psi.dvr.redu{m,3}(ii,jj) = sum ( sum ( sum ( conj(psi.dvr.grid_ND{m}(:,:,ii,:)) .* psi.dvr.grid_ND{m}(:,:,jj,:) ) ) );
    end
end

% Fourth coordinate
for ii=1:space.dof{4}.n_pts
    for jj=1:space.dof{4}.n_pts
        psi.dvr.redu{m,4}(ii,jj) = sum ( sum ( sum ( conj(psi.dvr.grid_ND{m}(:,:,:,ii)) .* psi.dvr.grid_ND{m}(:,:,:,jj) ) ) );
    end
end

%------------------------------------------------------------------------
% Reduce total density matrix (from 5-dimensional wavepacket propagation)
% to k-th dimension (k=1,2,3,4) by tracing over all the remaining dimensions
%-------------------------------------------------------------------------
function rho_redu_5d (m)
global psi space

% First coordinate
for ii=1:space.dof{1}.n_pts
    for jj=1:space.dof{1}.n_pts
        psi.dvr.redu{m,1}(ii,jj) = sum ( sum ( sum ( sum ( conj(psi.dvr.grid_ND{m}(ii,:,:,:,:)) .* psi.dvr.grid_ND{m}(jj,:,:,:,:) ) ) ) );
    end
end

% Second coordinate
for ii=1:space.dof{2}.n_pts
    for jj=1:space.dof{2}.n_pts
        psi.dvr.redu{m,2}(ii,jj) = sum ( sum ( sum ( sum ( conj(psi.dvr.grid_ND{m}(:,ii,:,:,:)) .* psi.dvr.grid_ND{m}(:,jj,:,:,:) ) ) ) );
    end
end

% Third coordinate
for ii=1:space.dof{3}.n_pts
    for jj=1:space.dof{3}.n_pts
        psi.dvr.redu{m,3}(ii,jj) = sum ( sum ( sum ( sum ( conj(psi.dvr.grid_ND{m}(:,:,ii,:,:)) .* psi.dvr.grid_ND{m}(:,:,jj,:,:) ) ) ) );
    end
end

% Fourth coordinate
for ii=1:space.dof{4}.n_pts
    for jj=1:space.dof{4}.n_pts
        psi.dvr.redu{m,4}(ii,jj) = sum ( sum ( sum ( sum ( conj(psi.dvr.grid_ND{m}(:,:,:,ii,:)) .* psi.dvr.grid_ND{m}(:,:,:,jj,:) ) ) ) );
    end
end

% Fifth coordinate
for ii=1:space.dof{5}.n_pts
    for jj=1:space.dof{5}.n_pts
        psi.dvr.redu{m,5}(ii,jj) = sum ( sum ( sum ( sum ( conj(psi.dvr.grid_ND{m}(:,:,:,:,ii)) .* psi.dvr.grid_ND{m}(:,:,:,:,jj) ) ) ) );
    end
end

%------------------------------------------------------------
% Plot reduced density matrices for k-th dimension 
%------------------------------------------------------------
function redu_dens ( step, k )       
global info plots psi space expect hamilt
persistent rho_max

% Initialize calculation of max density
if step==1
    rho_max(k) = 0;
end

% Main loop over (coupled) densities
for m=1:hamilt.coupling.n_eqs

    % Find maximal density at first step
    if step==1
        rho_max(k) = max ( rho_max(k), max(max(abs(psi.dvr.redu{m,k}))) );
    end

    % Contour plots

    if expect.ind.pop{m}(step)>expect.min_pop
        contour ( ...
            space.dvr.grid_1D{k}, ...
            space.dvr.grid_1D{k}, ...
            abs(psi.dvr.redu{m,k}), ...
            linspace(0, rho_max(k), plots.density.contour.nlev(1)), ...   % use even number of contours to avoid zero!
            'LineStyle', '-', ...
            'LineWidth', plots.style.line.thin, ...
            'Color',     plots.style.colors(m,:)   );
    end

    if m==1
        hold on
    end

end
hold off

% Dotted line along diagonal
line ( [ space.dof{k}.dvr_min space.dof{k}.dvr_max  ], ...
       [ space.dof{k}.dvr_min space.dof{k}.dvr_max ], ...
       'LineStyle', ':', ...
       'Color', 'k', ...
       'LineWidth', plots.style.line.thin)

% Axes, labels, etc
axis ( [ space.dof{k}.dvr_min space.dof{k}.dvr_max space.dof{k}.dvr_min space.dof{k}.dvr_max ] )
set ( gca, 'LineWidth',     plots.style.line.thick, ...
           'FontName',      plots.style.font.name,  ...
           'FontSize',      plots.style.font.large, ...
           'FontWeight',    plots.style.font.heavy )

% place the header approximately in the middle
if k==ceil(space.size.n_dim/2)
    title ( {info.header1;info.header2} )
end

if space.size.n_dim>1
    xlabel ( ['R_{', space.dof{k}.label, '} [a_0]'] )
    ylabel ( ['R_{', space.dof{k}.label, '} [a_0]'] )
else
    xlabel ( 'R [a_0]' )
    ylabel ( 'R [a_0]' )
end
    
% Negative imaginary potential (as absorbing boundary conditions)
if ~isempty(hamilt.nip.grid)
    if hamilt.nip.params.min(k) > space.dof{k}.dvr_min % Left/lower border
        line ( [ hamilt.nip.params.min(k) hamilt.nip.params.min(k) ], ...
               [ space.dof{k}.dvr_min space.dof{k}.dvr_max ], ...
               'LineStyle', '--', ...
               'Color', 'k', ...
               'LineWidth', plots.style.line.thin)
        line ( [ space.dof{k}.dvr_min space.dof{k}.dvr_max ], ...
               [ hamilt.nip.params.min(k) hamilt.nip.params.min(k) ], ...
               'LineStyle', '--', ...
               'Color', 'k', ...
               'LineWidth', plots.style.line.thin)
    end
    
    if hamilt.nip.params.max(k) < space.dof{k}.dvr_max % Right/upper border
        line ( [ hamilt.nip.params.max(k) hamilt.nip.params.max(k) ], ...
               [ space.dof{k}.dvr_min space.dof{k}.dvr_max ], ...
               'LineStyle', '--', ...
               'Color', 'k', ...
               'LineWidth', plots.style.line.thin)
        line ( [ space.dof{k}.dvr_min space.dof{k}.dvr_max ], ...
               [ hamilt.nip.params.max(k) hamilt.nip.params.max(k) ], ...
               'LineStyle', '--', ...
               'Color', 'k', ...
               'LineWidth', plots.style.line.thin)
    end
    
end

%--------------------------------------------------------------------
% Plot purity measure for k-th dimension vs. time
%--------------------------------------------------------------------
function purity ( step, k )
global expect info plots time psi space hamilt

% Calculate tr(rho^2)
expect.tot.pur(step,k) = 0;
for m = 1:hamilt.coupling.n_eqs
    if expect.ind.pop{m}(step) > expect.min_pop
        expect.ind.pur{m}(step,k) = abs(trace(psi.dvr.redu{m,k}^2));
        expect.tot.pur   (step,k) = expect.tot.pur(step,k) + expect.ind.pur{m}(step,k) / expect.ind.pop{m}(step)^2;
    end
end

% Plot evolving curve
mask = find ( expect.tot.pur(:,k) > eps );
plot ( time.main.grid (mask), ...
       expect.tot.pur(mask,k), ...
       'LineStyle', '-', ...
       'LineWidth', plots.style.line.thick, ...
       'Color',     'k' )

% Line styles and fonts
set ( gca, 'LineWidth',  plots.style.line.thick, ...
           'FontName',      plots.style.font.name,  ...
           'FontSize',   plots.style.font.large, ...
           'FontWeight', plots.style.font.heavy )   

% Axes range and labels
axis ( [ 0 time.main.total -0.1 1.1 ] )
if strcmpi(info.program,'qm_propa')
    xlabel ('t [hbar / E_h]')
elseif strcmpi(info.program,'qm_bound')
    xlabel ('n')
end
    
switch space.size.n_dim
    case 1
        ylabel ( 'tr(\rho^2)' )
    case 2
        ylabel ( 'tr(\rho_1^2) = tr(\rho_2^2)' )
    otherwise
        ylabel ( ['tr(\rho_' int2str(k) '^2)'] )
end   



