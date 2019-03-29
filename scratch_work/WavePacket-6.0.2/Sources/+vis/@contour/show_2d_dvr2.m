%------------------------------------------------------------
% Densities projected on second coordinate (horizontal plot)
%------------------------------------------------------------
function show_2d_dvr2 ( obj, state, step )       
global expect hamilt space
persistent rho_max

% Preallocate
rho  = cell (hamilt.coupling.n_eqs,1);

% Get reduced densities where populations not too small
for  m=1:hamilt.coupling.n_eqs
    if expect.pop.cha{m}(step)>expect.min_pop
        switch(lower(class(state)))
            case ('traj')
                [c,~] = histcounts (state.pos{2}(state.cha==m), space.dof{2}.x_grid);
                rho{m} = c';
            case ('wave')
                rho{m} = sum( abs(state.dvr{m}).^2 .* space.weight, 1 );
        end
    end
end
        
% Get maximum density
% if step==1
rho_max = 0;
for  m=1:hamilt.coupling.n_eqs
    if expect.pop.cha{m}(step)>expect.min_pop
        rho_max = max ( rho_max, max(rho{m}) );
    end
end
% end

% Plot projected densities where populations not too small
for  m=1:hamilt.coupling.n_eqs
    if expect.pop.cha{m}(step)>expect.min_pop
        switch(lower(class(state)))
            case ('traj')
                x  = space.dof{2}.x_grid;
                n  = space.dof{2}.n_pts;
                dx = space.dof{2}.x_dlt;
                h = plot ( rho{m}/rho_max, x(1:n-1)+dx/2 );
            case ('wave')
                h = plot ( rho{m}/rho_max, space.dof{2}.x_grid );
        end
        h.LineStyle = obj.patterns{1};
        h.LineWidth = obj.l_thick;
        h.Color     = obj.colors(m,:);
    end
end

% Axes and labels
h = gca;
h.XAxisLocation = 'top';
h.YAxisLocation = 'right';
h.LineWidth     = obj.l_thick;
h.FontName      = obj.f_name;
h.FontSize      = obj.f_large;
h.FontWeight    = obj.f_heavy;

xlabel ( [ '\rho (R_{', space.dof{2}.label, '})' ] )
ylabel ( [ 'R_{', space.dof{2}.label, '}' ] )
if ~obj.range
    axis ( [ 0 1 space.dof{2}.dvr_min space.dof{2}.dvr_max ] )
else
    axis ( [ 0 1 obj.y_min obj.y_max ] )
end

% Absorbing boundary conditions
if isa(state,'wave') && isfield (hamilt,'nip')
    for  m=1:hamilt.coupling.n_eqs
        if ~isempty(hamilt.nip{m}.dvr)
            if hamilt.nip{m}.min(2) > space.dof{2}.dvr_min % Lower
                h = line ( [ 0 1 ], [ hamilt.nip{m}.min(2) hamilt.nip{m}.min(2) ] );
                 h.LineStyle= obj.patterns{2};
                h.Color= obj.colors(m,:);
                h.LineWidth= obj.l_thin;
            end
            
            if hamilt.nip{m}.max(2) < space.dof{2}.dvr_max % Upper
                h = line ( [ 0 1 ], [ hamilt.nip{m}.max(2) hamilt.nip{m}.max(2) ] );
                h.LineStyle= obj.patterns{2};
                h.Color= obj.colors(m,:);
                h.LineWidth= obj.l_thin;
            end
        end
    end
end

