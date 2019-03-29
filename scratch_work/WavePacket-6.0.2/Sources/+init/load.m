%------------------------------------------------------------------------------
%
% Loads the initial wavefunction from a file generated
% in a previous run of qm_propa or qm_bound. The user has to
% specify the directory, filename, and (time step) index,
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2017 - .... Burkhard Schmidt
%               2008 - 2010 Ulf Lorenz
%               
%
% see the README file for license details.

classdef load < handle

    properties (Access = public)
        dir         % Name of directory  where the savefile resides
        file        % Name of the file where the wavefunction is stored
        index       % If set, gives the index of time step or bound state
        channel     % If set, gives the index/indices of the channel(s)
    end

    properties (Access = private)
    end
    
    methods (Access = public)
        
        % Constructor: Set default values
        function obj = load
            obj.dir = [];
            obj.file = [];
            obj.index = [];
            obj.channel = [];
        end
        
        % Initialize loading: Set/check parameters
        function init_corr (obj)
            global hamilt
            if isempty (obj.index)
                log.error ('Index of time step or bound state has to be specified');
            end
            if isempty (obj.channel)
                obj.channel = 1:hamilt.coupling.n_eqs;
            end
        end
        
        % Display Gaussian, overloading default disp method
        function disp(obj)
            log.disp ('Loading initial wavefunction from file            ')
            log.disp ('***************************************************************')
            log.disp (' ')
            log.disp ( ['Directory name           : '         obj.dir     ])
            log.disp ( ['File name                : '         obj.file    ])
            log.disp ( ['Time step or bound state : ' num2str(obj.index)  ])
            log.disp ( ['Channel(s) to load       : ' num2str(obj.channel)])
        end

        % Evaluate wave function on a grid
        function wave_corr (obj, state)
            
            global hamilt space
            
            % Load wavefunction from previous calculation
            if ~isempty(obj.dir)
                state.sav_dir  = obj.dir;
            end
            if ~isempty(obj.file)
                state.sav_file = obj.file;
            end
            i_o.load_0 (state, false);
            i_o.load_n (state, obj.index);
            
            % Check dimensionality
            if size(state.dvr{1}) ~= size(space.dvr{1})
                log.error('Loaded initial wavefunction has wrong dimensions')
            end
            
            % Check number of coupled channels
            if length(state.dvr) ~= hamilt.coupling.n_eqs
                log.error('Loaded initial wavefunction has different number of channels')
            end
            
            % Zero wavefunctions, except in channels specified
            for m = 1:hamilt.coupling.n_eqs
                if ~any(m==obj.channel)
                    state.dvr{m} = zeros(size(state.dvr{1}));
                end
            end
                                         
        end
        
    end
end
