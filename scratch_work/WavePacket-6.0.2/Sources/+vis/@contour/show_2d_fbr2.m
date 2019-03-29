%------------------------------------------------------------
% Densities projected on second coordinate (horizontal plot)
% FBR plot
%------------------------------------------------------------
function show_2d_fbr2 ( obj, state, step )          
global expect hamilt space
persistent rho_max

% Preallocation
rho  = cell (hamilt.coupling.n_eqs,1);

% Get reduced densities where populations not too small
for  m=1:hamilt.coupling.n_eqs
    if expect.pop.cha{m}(step)>expect.min_pop
        switch(lower(class(state)))
            case ('traj')
                [c,~] = histcounts (state.mom{2}(state.cha==m), space.dof{2}.p_grid);
                rho{m} = c';
            case ('wave')
                fbr = dvr2fbr(space.dof{1}, state.dvr{m});
                fbr = dvr2fbr(space.dof{2}, fbr);
                rho{m} = sum( abs(fbr).^2, 1 );
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
                y  = space.dof{2}.p_grid;
                n  = space.dof{2}.n_pts;
                dy = space.dof{2}.p_dlt;
                h = plot ( rho{m}/rho_max, y(1:n-1)+dy/2 );
            case ('wave')
                h = plot ( rho{m}/rho_max, space.dof{2}.p_grid );
        end
        h.LineStyle = obj.patterns{1};
        h.LineWidth = obj.l_thick;
        h.Color     = obj.colors(m,:);
    end
end

if ~obj.range
    axis ( [ 0 1 space.dof{2}.fbr_min space.dof{2}.fbr_max ] )
else
    axis ( [ 0 1 obj.y_min obj.y_max ] )
end

% Axes and labels
h = gca;
h.XAxisLocation = 'top';
h.YAxisLocation = 'right';
h.LineWidth     = obj.l_thick;
h.FontName      = obj.f_name;
h.FontSize      = obj.f_large;
h.FontWeight    = obj.f_heavy;

xlabel ( [ '\rho (P_{', space.dof{2}.label, '})' ] )
ylabel ( [ 'P_{', space.dof{2}.label, '}' ] )
