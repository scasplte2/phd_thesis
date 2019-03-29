%------------------------------------------------------------
% Densities projected on first coordinate (vertical plot)
%------------------------------------------------------------
function show_2d_dvr1 ( obj, state, step )       
global expect hamilt space
persistent rho_max

% Preallocate
rho = cell (hamilt.coupling.n_eqs,1);

% Get reduced densities where populations not too small
for  m=1:hamilt.coupling.n_eqs
    if expect.pop.cha{m}(step)>expect.min_pop
        switch(lower(class(state)))
            case ('traj')
                [c,~] = histcounts (state.pos{1}(state.cha==m), space.dof{1}.x_grid);
                rho{m} = c';
            case ('wave')
                rho{m} = sum( abs(state.dvr{m}).^2 .* space.weight, 2 );
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
                x  = space.dof{1}.x_grid;
                n  = space.dof{1}.n_pts;
                dx = space.dof{1}.x_dlt;
                h = plot ( x(1:n-1)+dx/2, rho{m}/rho_max ); 
            case ('wave')
                h = plot ( space.dof{1}.x_grid, rho{m}/rho_max );
        end
        h.LineStyle = obj.patterns{1};
        h.LineWidth = obj.l_thick;
        h.Color     = obj.colors(m,:);
    end
end

% Axes and labels
h = gca;
h.LineWidth     = obj.l_thick;
h.FontName      = obj.f_name;
h.FontSize      = obj.f_large;
h.FontWeight    = obj.f_heavy;

xlabel ( [ 'R_{', space.dof{1}.label, '}' ] )
ylabel ( [ '\rho (R_{', space.dof{1}.label, '})' ] )

if ~obj.range
    axis ( [ space.dof{1}.dvr_min space.dof{1}.dvr_max 0 1 ] )
else
    axis ( [ obj.x_min obj.x_max 0 1 ] )
end

% Absorbing boundary conditions
if isa(state,'wave') && isfield (hamilt,'nip')
    for  m=1:hamilt.coupling.n_eqs
        if ~isempty(hamilt.nip{m}.dvr)
            if hamilt.nip{m}.min(1) > space.dof{1}.dvr_min % Left
                h = line ( [ hamilt.nip{m}.min(1) hamilt.nip{m}.min(1) ], [ 0 1 ] );
                h.LineStyle= obj.patterns{2};
                h.Color= obj.colors(m,:);
                h.LineWidth= obj.l_thin;
            end
            
            if hamilt.nip{m}.max(1) < space.dof{1}.dvr_max % Right
                h = line ( [ hamilt.nip{m}.max(1) hamilt.nip{m}.max(1) ], [ 0 1 ] );
                h.LineStyle= obj.patterns{2};
                h.Color= obj.colors(m,:);
                h.LineWidth= obj.l_thin;
            end
        end
    end
end

