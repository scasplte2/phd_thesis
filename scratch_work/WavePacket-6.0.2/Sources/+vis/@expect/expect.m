%--------------------------------------------------------------------------
%
% Compose figure with two or three subplots
% displaying various expectation values as time series
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%
% see the README file for license details.

classdef expect < vis.styles & handle
    
    properties (Access = public)
        
        w_left     % Left border of plot window
        w_lower    % Left border of plot window
        w_width    % Width of plot window
        w_height   % Height of plot window

        errorbar    % Toggle plotting errorbars
        legends     % Toggle curve legends
        logo        % Toggle logos in all four corners
        wide        % Toggle wide format: 16:9
        
        p_min       % Lower bound for population
        p_max       % Upper bound for populations
        
        e_min       % Lower bound for energies
        e_max       % Upper bound for energies
        
        export      % Toggle graphics export
        file        % Set the output to a custom filename (w/o suffix)
        
    end
    
    properties (Access = private)
        
        mask_tot    % Fake: Find steps already existing
        mask_cha    % Find steps where populations of individual channels are non-negligible

    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = expect
            
            obj = obj@vis.styles;
            
            % Position and size of plot window: Format 7:9
            obj.w_left   = round(11*obj.s_height/16);  % Left border
            obj.w_lower  = round(01*obj.s_height/16);  % Lower border
            obj.w_width  = round(07*obj.s_height/16);  % Width
            obj.w_height = round(09*obj.s_height/16);  % Height
            obj.w_width =  2*round(obj.w_width /2);    % Should be even
            obj.w_height = 2*round(obj.w_height/2);    % Should be even
            
            % Appearence of plot
            obj.errorbar = false;        % Toggle plotting errorbars
            obj.legends  = true;         % Toggle legends
            obj.logo     = true;         % Toggle logos in all four corners
            obj.wide = [];               % Toggle wide format: 16:9
            
            % Range settings for populations
            obj.p_min = -0.1;            % Minimum of population plot
            obj.p_max = +1.1;            % Maximum of population plot
            
            % Range settings for energies
            obj.e_min   = [];            % Lower bound of the energy plot
            obj.e_max   = [];            % Upper bound of the energy plot
            
            % Export to graphics file: jpg/tiff/epsc/pdf, etc ...
            obj.export  = false;         % Toggle graphics export
            obj.file    = [];            % Set custom filename (suffix determines image type)
            
            % Masks where individual/total populations do exist
            obj.mask_tot = [];           % Fake: Find steps already existing
            obj.mask_cha = [];           % Find steps where individual populations 
           
        end
        
        %--------------------------------------------
        % Show expectation values
        %--------------------------------------------
        function show_expect (obj)
            
            global expect hamilt info time
            
            if obj.wide % Wide format: 16:9
                w=16; h=09; o=9/16;
            else % Narrow format: 7:9
                w=07; h=09; o=0;
            end
            
            % Find steps populations of individual channels are non-negligible
            for m=1:hamilt.coupling.n_eqs
                obj.mask_cha{m} = find ( expect.pop.cha{m} > expect.min_pop );
            end
            
            % Fake: Find steps already existing
            obj.mask_tot = find ( expect.pop.tot > expect.min_pop );
            
            % Time dependent simulations
            switch lower (info.program)
                case {'qm_propa'}
                    
                    % Populations/projections vs. time
                    subplot ( 'Position', [o+1/w (1+14/3)/h 5/w 7/(3*h)] )
                    cla;
                    population (obj)
                    if obj.wide; set ( gca, 'YAxisLocation', 'right'); end
                    
                    % Expectation values: Energies versus time
                    subplot ( 'Position', [o+1/w (1+7/3)/h 5/w 7/(3*h)] )
                    cla;
                    energies (obj)
                    if obj.wide; set ( gca, 'YAxisLocation', 'right'); end
                    
                    if isfield(time,'pulse')
                        % External electric field versus time
                        subplot ( 'Position',[o+1/w 1/h 5/w 7/(3*h)] )
                        cla;
                        efield (obj)
                        if obj.wide; set ( gca, 'YAxisLocation', 'right'); end
                    elseif ~isempty(time.steps.acf)
                        % Autocorrelation function versus time
                        subplot ( 'Position', [o+1/w 1/h 5/w 7/(3*h)] )
                        cla;
                        correlation (obj)
                        if obj.wide; set ( gca, 'YAxisLocation', 'right'); end
                    end
                
                    % Time independent simulations
                case {'qm_bound'}
                    
                    % Populations/projections and energies
                    if isfield(hamilt, 'amo')
                        subplot ( 'Position', [o+1/w 5/h 5/w 3/h] )
                        population (obj)
                        if obj.wide; set ( gca, 'YAxisLocation', 'right'); end
                        
                        subplot ( 'Position', [o+1/w 1/h 5/w 3/h] )
                        energies (obj)
                        if obj.wide; set ( gca, 'YAxisLocation', 'right'); end
                    else
                        
                        % Energies only
                        subplot ( 'Position', [o+1/w 1/h 5/w 7/h] )
                        energies (obj)
                        if obj.wide; set ( gca, 'YAxisLocation', 'right'); end
                    end
                    
                otherwise
                    log.error ('Wrong choice of program')
                
            end
            
        end
 
 
        % Details of surface hopping trajectory methodology
        function show_detail (obj,state)
            
            global time
            
            if obj.wide % Wide format: 16:9
                w=16; h=09; o=9/16;
            else % Narrow format: 7:9
                w=07; h=09; o=0;
            end
            
            if strcmpi(class(state),'wave')
                return
            end
            
            pos_siz = [o+1/w 1.5/h 5/w 1/h];  % x y w h
            
            % First string: Q/C method
            if isfield (time,'hop')
                hop_meth = class(time.hop);
                hop_meth = upper(hop_meth(5:end));
                string1 = [ 'Quantum-Classical method      : ' hop_meth ];
            else
                string1 = 'Purely classical trajectories';
            end
            
            % Second string: Number of trajectories
            string2 = ['Number of trajectories used   : ' num2str(state.n_p)];
            
            % Third string: Seeding random numbers
            if ~isempty (state.rnd_seed)
                string3 = ['Seeding random numbers        : ' num2str(state.rnd_seed)];
            else
                string3 = 'Not seeding random numbers';
            end
            
            % Fourth string: Momentum rescaling
            if isfield (time,'hop')
                if time.hop.rescale
                    if time.hop.sca_nac
                        string4 = 'Momentum rescaling            : along NAC vector';
                    else
                        string4 = 'Momentum rescaling            : along momenta';
                    end
                else
                    string4 = 'Momentum rescaling    : none';
                end
            else
                string4 = ' ';
            end
            
            % Fifth string: Propagation scheme
            prop_meth = class(time.propa);
            prop_meth = prop_meth (5:end);
            string5 = ['Propagation scheme            : ' prop_meth ', order ' int2str(time.propa.order)];
            
            h = annotation('textbox');
            h.Position            = pos_siz;
            h.HorizontalAlignment = 'left';
            h.Interpreter         = 'none';
            h.LineStyle           = 'none';
            h.String              = { string1, string2, string3, string4, string5 };
            h.FitHeightToText     = 'on';
            h.FontName            = obj.f_name;
            h.FontSize            = obj.f_large;
            h.FontWeight          = obj.f_heavy;
            
        end
        
    end
    
    
    methods (Access = private)
        
        population  (obj)
        energies     (obj)
        correlation (obj)
        efield      (obj)
        
    end
    
end
