%--------------------------------------------------------------------------
% Create animated graphics of time evolving densities
% Create animated graphics of time evolving expectation values
% Create graphics of absorption spectrum (FFT of autocorrelation)
%
% If "state" is an object of class "wave": densities from wavefunctions
%
% If "state" is an object of class "traj": densities from trajectories
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2018 Burkhard Schmidt's group
%               2007-2010 Ulf Lorenz
%
% see the README file for license details.

function show ( state, step )
global expect hamilt plots space time
persistent writerObj

% Fake time axis for the case of only 1 step
if time.steps.t_total==0
    time.steps.t_total=1;
end

% Various initializations of density plots
if step==1 && isfield(plots,'density')
    
    %% Q/C only: Find maxima of binned densities
    if isa(state,'traj')
        if ~isa(plots.density,'vis.reduced_1d') && ~isa(plots.density,'vis.reduced_2d')
            
            dvr_max = 0;
            fbr_max = 0;
            wig_max = 0;
            
            for m=1:hamilt.coupling.n_eqs
                if expect.pop.cha{m}(step)>expect.min_pop
                    
                    switch space.n_dim
                        case 1
                            [c_dvr,~] = histcounts  (state.pos{1},               space.dof{1}.x_grid                     );
                            [c_fbr,~] = histcounts  (              state.mom{1},                      space.dof{1}.p_grid);
                            [c_wig,~] = histcounts2 (state.pos{1}, state.mom{1}, space.dof{1}.x_grid, space.dof{1}.p_grid);
                        case 2
                            [c_dvr,~] = histcounts2 (state.pos{1}, state.pos{2}, space.dof{1}.x_grid, space.dof{2}.x_grid);
                            [c_fbr,~] = histcounts2 (state.mom{1}, state.mom{2}, space.dof{1}.p_grid, space.dof{2}.p_grid);
                            c_wig    = 0; % dummy
                        otherwise
                            log.error ('Binning of trajectories in more than 2 dimensions not yet implemented')
                    end
                    
                    dvr_max = max ( dvr_max, max(c_dvr(:)) );
                    fbr_max = max ( fbr_max, max(c_fbr(:)) );
                    wig_max = max ( wig_max, max(c_wig(:)) );
                    
                end
            end
            
        else
            
            dvr_max = zeros (space.n_dim,1);
            fbr_max = zeros (space.n_dim,1);
            wig_max = zeros (space.n_dim,1);
            
            for m=1:hamilt.coupling.n_eqs
                if expect.pop.cha{m}(step)>expect.min_pop
                    
                    for k=1:space.n_dim
                        [c_dvr,~] = histcounts  (state.pos{1},               space.dof{1}.x_grid                     );
                        [c_fbr,~] = histcounts  (              state.mom{1},                      space.dof{1}.p_grid);
                        [c_wig,~] = histcounts2 (state.pos{1}, state.mom{1}, space.dof{1}.x_grid, space.dof{1}.p_grid);
                        
                        dvr_max(k) = max ( dvr_max(k), max(c_dvr(:)) );
                        fbr_max(k) = max ( fbr_max(k), max(c_fbr(:)) );
                        wig_max(k) = max ( wig_max(k), max(c_wig(:)) );
                    end
                    
                end
            end
            
        end
        
        if isempty(plots.density.scale_dvr)
            plots.density.scale_dvr = dvr_max;
        end
        if isempty(plots.density.scale_fbr)
            plots.density.scale_fbr = fbr_max;
        end
        if isempty(plots.density.scale_wig)
            plots.density.scale_wig = wig_max;
        end
        
    end
    
    %% Q/M only: Determine ranges of kinetic/potential/total energy
    if isa(state,'wave') || ( isa(state,'traj') && space.n_dim < 3 )
        if isempty(plots.density.kin_min)
            plots.density.kin_min = 0;
        end
        if isempty(plots.density.kin_max)
            plots.density.kin_max = hamilt.kin_max;
        end
        plots.density.kin_delta = plots.density.kin_max - plots.density.kin_min;
        
        if isempty(plots.density.pot_min)
            plots.density.pot_min = hamilt.pot_min;
        end
        if isempty(plots.density.pot_max)
            plots.density.pot_max = hamilt.pot_max;
        end
        if plots.density.pot_min==plots.density.pot_max % Dirty trick for free particle
            plots.density.pot_max = plots.density.kin_max;
        end
        plots.density.pot_delta = plots.density.pot_max - plots.density.pot_min;
        
        plots.density.tef_min = min(plots.density.pot_min,0);
        plots.density.tef_max = plots.density.pot_max + plots.density.kin_max;
        plots.density.tef_delta = plots.density.tef_max - plots.density.tef_min;
    end
    
    % Determine maximal values of densities in pos/mom representation
    if isa (state,'wave') && space.n_dim <= 3
        
        dvr_max = 0;
        fbr_max = 0;
        
        for m=1:hamilt.coupling.n_eqs
            if expect.pop.cha{m}(step)>expect.min_pop
                dvr_max = max ( dvr_max, max(abs(state.dvr{m}(:)).^2) );
                
                fbr = state.dvr{m};
                for k = 1:space.n_dim
                    fbr = dvr2fbr(space.dof{k}, fbr);
                end
                fbr_max = max ( fbr_max, max(abs(fbr(:)).^2) );
                
            end
        end
        if isempty(plots.density.scale_dvr)
            plots.density.scale_dvr = dvr_max;
        end
        if isempty(plots.density.scale_fbr)
            plots.density.scale_fbr = fbr_max;
        end
        
    end
    
end

%% Animate densities in DVR/FBR/phase space

% Toggle plotting
if isfield (plots,'density')
    
    % First figure
    h1=figure(1);
    
    % First call only
    if step==1
        
        % Clear figure and set figure size
        clf;
        if isfield (plots,'expect')
            set(h1,'units','pixels', ...
                'position',[...
                plots.density.w_left ...
                plots.density.w_lower ...
                plots.density.w_width + plots.expect.w_width ...
                plots.density.w_height] );
            plots.density.wide = true;
        else
            set(h1,'units','pixels', ...
                'position',[...
                plots.density.w_left ...
                plots.density.w_lower ...
                plots.density.w_width ...
                plots.density.w_height] );
            plots.density.wide = false;
        end
        
        % Logos in all four corners of the plots
        if plots.density.logo
            show_logo (plots.density)
        end
        
        % Initialize movie export (first step only)
        if plots.density.export
            
            if isempty(plots.density.file)
                switch lower(class(state))
                    case 'wave'
                        plots.density.file = 'wave';
                    case 'traj'
                        if isfield (time,'hop')
                            hop_meth = class(time.hop);
                            plots.density.file = hop_meth(5:end);
                        else
                            plots.density.file = 'traj';
                        end
                end
                plot_type = class(plots.density);
                plot_type = plot_type (5:end);
                plots.density.file = strcat(plots.density.file,'_',plot_type);
            end
            
            log.disp ('***************************************************************')
            if ~plots.density.images
                if ispc||ismac
                    extension = '.mp4';
                    myprofile = 'MPEG-4';
                elseif isunix
                    extension = '.avi';
                    myprofile = 'Motion JPEG AVI';
                end
                log.disp ( ['Saving animated density plot : ' strcat(plots.density.file,extension)] )
                if exist (strcat(plots.density.file,extension), 'file')
                    delete( strcat(plots.density.file,extension) );
                end
                close(writerObj);
                writerObj = VideoWriter (plots.density.file, myprofile);
                open(writerObj);
            else
                log.disp( 'Saving animated density plot as sequence of jpeg files')
            end
            log.disp ('***************************************************************')
        end
    end
    
    % Various types of plots (with/without marginals)
    show_plot (plots.density, state, step)
    
    % Draw also expectation values (optionally)
    if isfield (plots,'expect')
        plots.expect.wide = true;
        show_expect (plots.expect);
        show_detail (plots.expect,state)
    end
    
    % Info about rendering, double buffering
    if step==1
        r = get (gcf, 'Renderer');
        d = get (gcf, 'DoubleBuffer');
        log.disp (' ')
        log.disp (['Type of density plot             : ' class(plots.density) ])
        log.disp (['Rendering method                 : ' r])
        log.disp (['Double buffering (painters only) : ' d])
        log.disp (' ')
    end
    
    % Save last snapshot
    if step==time.steps.m_number
        if plots.density.export
            full_file_name1 = strcat(plots.density.file, '.jpg');
            full_file_name2 = strcat(plots.density.file, '.fig');
            log.disp ( '***************************************************************' )
            log.disp ( ['Saving last snapshot to file : ' full_file_name1] )
            log.disp ( ['Saving last snapshot to file : ' full_file_name2] )
            log.disp ( '***************************************************************' )
            log.disp (' ')
            saveas(gcf,full_file_name1)
            saveas(gcf,full_file_name2)
        end
    end
    
    % Add current figure as movie frame. If we export single images, we output
    % a new file instead
    if plots.density.export
        if ~plots.density.images
            frame = getframe(gcf);
            writeVideo(writerObj,frame);
        else
            % pad the filename with as many zeros as required
            padding = floor(log10(time.steps.m_number)) + 1;
            pattern = strcat('%s%0', num2str(padding), 'i.jpg');
            full_file_name = sprintf(pattern, plots.density.file, step);
            saveas(gcf, full_file_name);
        end
    end
    
    % Close/clear movie object (last step only)
    if plots.density.export && step==time.steps.m_number
        if ~plots.density.images
            close(writerObj);
        end
    end
    
end

%% Animate expectation values (with uncertainties) if no densities are displayed

% Toggle plotting
if isfield(plots,'expect') && ~isfield(plots,'density')
    
    % Second figure
    h2=figure(2);
    
    %% Upon first call only ...
    if step == 1
        
        % Clear figure and set size
        clf
        set(h2,'units','pixels', ...
            'position',[plots.expect.w_left ...
            plots.expect.w_lower ...
            plots.expect.w_width ...
            plots.expect.w_height] );
        
        % Logos in the corners of the plots
        if plots.expect.logo
            show_logo (plots.expect)
        end
        
    end
    
    % Draw expectation values
    plots.expect.wide = false;
    show_expect (plots.expect);
    show_detail (plots.expect,state)
    
    % Last step only
    if step==time.steps.m_number
        
        
        % Export graphics to file (if desired)
        if plots.expect.export
            if isempty(plots.expect.file)
                switch lower(class(state))
                    case 'wave'
                        plots.expect.file = 'wave';
                    case 'traj'
                        if isfield (time,'hop')
                            hop_meth = class(time.hop);
                            plots.expect.file = hop_meth(5:end);
                        else
                            plots.expect.file = 'traj';
                        end
                end
                plot_type = class(plots.expect);
                plot_type = plot_type (5:end);
                plots.expect.file = strcat(plots.expect.file,'_',plot_type);
                full_file_name1 = strcat(plots.expect.file, '.jpg');
                full_file_name2 = strcat(plots.expect.file, '.fig');
            end
            log.disp ( '***************************************************************' )
            log.disp ( ['Saving plots of expectation values: ' full_file_name1] )
            log.disp ( ['Saving plots of expectation values: ' full_file_name2] )
            log.disp ( '***************************************************************' )
            log.disp ( ' ' )
            saveas(gcf,full_file_name1)
            saveas(gcf,full_file_name2)
        end
        
    end
    
    drawnow;
end


%% Q/M only: Absorption spectrum 
if isa (state,'wave')
    
    % Toggle plotting (last step only, TDSE only)
    if isfield (plots,'spectrum') && step==time.steps.m_number && step>1
        
        % Fourier transform of autocorrelation
        spectrum (time.steps);
        
        % Second figure
        h3=figure(3);
        
        % Clear figure and set size
        clf
        set(h3,'units','pixels', ...
            'position',[plots.spectrum.w_left ...
            plots.spectrum.w_lower ...
            plots.spectrum.w_width ...
            plots.spectrum.w_height] );
        
        % Logos in the corners of the plots
        if plots.spectrum.logo
            show_logo (plots.spectrum)
        end
        
        % Draw spectrum
        show_spec (plots.spectrum);
        
        % Export graphics to file (if desired)
        if plots.spectrum.export
            if isempty(plots.spectrum.file)
                switch lower(class(state))
                    case 'wave'
                        plots.spectrum.file = 'wave';
                    case 'traj'
                        plots.spectrum.file = state.q_c;
                end
                plot_type = class(plots.spectrum);
                plot_type = plot_type (5:end);
                plots.spectrum.file = strcat(plots.spectrum.file,'_',plot_type);
                full_file_name1 = strcat(plots.spectrum.file, '.jpg');
                full_file_name2 = strcat(plots.spectrum.file, '.fig');
            end
            log.disp ( '***************************************************************' )
            log.disp ( ['Saving plot of spectrum to file : ' full_file_name1] )
            log.disp ( ['Saving plot of spectrum to file : ' full_file_name2] )
            log.disp ( '***************************************************************' )
            log.disp ( ' ' )
            saveas(gcf,full_file_name1)
            saveas(gcf,full_file_name2)
        end
        
        drawnow;
    end
end
