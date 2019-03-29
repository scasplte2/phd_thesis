%--------------------------------------------------------------------------
% Plot position density and potential energy curve (from wavefunctions)
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%
% see the README file for license details.

function show_1d_dvr ( obj, state, step )
global expect hamilt space


%% Loop over (coupled) states
for  m=1:hamilt.coupling.n_eqs
    
    % If population not too small
    if expect.pop.cha{m}(step) > expect.min_pop
        
        % Horizontal offset (= kinetic + potential energy)
        if obj.energy
            offset = expect.pot.cha{m}(step) + expect.kin.cha{m}(step);
            scaling = obj.pot_delta;
            y_min = obj.pot_min-obj.pot_delta/10;
            y_max = obj.pot_max+obj.pot_delta/10;

        else
            offset = 0;
            scaling = 1;
            y_min = -0.1;
            y_max = +1.1;
        end
      
        switch(lower(class(state)))
            
            case ('traj')
                
                % Sort trajectories into histogram
                x  = space.dof{1}.x_grid;
                n  = space.dof{1}.n_pts;
                dx = space.dof{1}.x_dlt;
                
                [rho,~] = histcounts (state.pos{1}(state.cha==m), x);
                
                h = plot ( x(1:n-1)+dx/2, rho*scaling/obj.scale_dvr+offset);
                h.LineStyle = obj.patterns{1};
                h.Color     = obj.colors(m,:);
                h.LineWidth = obj.l_thick;

            case ('wave')
                
                % Get density from wavefunction
                rho = abs   ( state.dvr{m} ) .^2;
                
                % Get phase of wavefunction: Map interval [-pi,pi] into [0,1]
                phi = angle ( state.dvr{m} ) / (2*(pi+0.001)) + 1/2;
                                
                vis.styles.plot_color ( space.dvr{1}, ...
                    rho*scaling/obj.scale_dvr, ...
                    phi, ...
                    obj.colors(m,:), ...
                    obj.l_thick, ...
                    offset, ...
                    0 )
                
        end
        
        % Plot horizontal base line
        if obj.energy
            h = line ( [space.dof{1}.dvr_min space.dof{1}.dvr_max], [1 1]*offset );
            h.LineStyle = obj.patterns{1};
            h.Color     = obj.colors(m,:);
            h.LineWidth = obj.l_thin;
        end
        
    end
    
    % Plot potential energy curve
    if obj.energy
        if ~isempty(hamilt.pot{m,m}.dvr)
            h = plot ( space.dvr{1}, hamilt.pot{m,m}.dvr );
        else
            h = plot ( space.dvr{1}, zeros(size(space.dvr{1})));
        end
        h.LineStyle = obj.patterns{1};
        h.Color     = obj.colors(m,:);
        h.LineWidth = obj.l_thin;
    end
    
end

%% Axes, labels, etc
if ~obj.range
    axis ( [ space.dof{1}.dvr_min space.dof{1}.dvr_max y_min y_max ] )
else
    axis ( [ obj.x_min obj.x_max y_min y_max ] )
end

h = gca;
h.LineWidth     = obj.l_thick;
h.FontName      = obj.f_name;
h.FontSize      = obj.f_large;
h.FontWeight    = obj.f_heavy;

xlabel ( [ 'R_{', space.dof{1}.label, '}' ] )

if obj.energy
    
    if hamilt.coupling.n_eqs==1
        ylabel ( 'V(R)' )
    else
        if strcmpi ( hamilt.coupling.represent,'adi')
            ylabel ( 'V_{adi}(R)' )
        elseif strcmpi ( hamilt.coupling.represent,'dia')
            ylabel ( 'V_{dia}(R)' )
        end
    end
    
else
    
    if hamilt.coupling.n_eqs==1
        ylabel ( '\rho(R)' )
    else
        if strcmpi ( hamilt.coupling.represent,'adi')
            ylabel ( '\rho_{adi}(R)' )
        elseif strcmpi ( hamilt.coupling.represent,'dia')
            ylabel ( '\rho_{dia}(R)' )
        end
    end
    
end

% Negative imaginary potential (as absorbing boundary condition)
if isa (state,'wave') && isfield (hamilt,'nip')
    for  m=1:hamilt.coupling.n_eqs
        if ~isempty(hamilt.nip{m}.dvr)
            if hamilt.nip{1}.min(1) > space.dof{1}.dvr_min
                h = line ( [hamilt.nip{m}.min(1) hamilt.nip{m}.min(1)], ...
                    [obj.pot_min-obj.pot_delta/10 obj.pot_max+obj.pot_delta/10] );
                h.LineStyle = obj.patterns{2};
                h.Color     = obj.colors(m,:);
                h.LineWidth = obj.l_thin;
            end
            
            if hamilt.nip{m}.max(1) < space.dof{1}.dvr_max
                h = line ( [hamilt.nip{m}.max(1) hamilt.nip{m}.max(1)], ...
                    [obj.pot_min-obj.pot_delta/10 obj.pot_max+obj.pot_delta/10] );
                h.LineStyle = obj.patterns{2};
                h.Color     = obj.colors(m,:);
                h.LineWidth = obj.l_thin;
            end
        end
    end
end

end
