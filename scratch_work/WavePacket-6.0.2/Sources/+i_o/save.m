% ------------------------------------------------------------------------
%
% Stores states (trajectories or wavefunctions) in an internal array and 
% finally saves it to MATLAB® formatted binary file (MAT-file, ....mat) 
%
% ------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2008-2010 Ulf Lorenz
%
% see the README file for license details.

function save ( state, step )

global expect hamilt info plots space time

persistent arr_dvr arr_pos arr_mom arr_cha index

%% If we do not want to save the states, exit immediately
if ~state.sav_export
    return;
end
    
%% Initialize everything
if step == 1
    
    % Output of file/directory
    log.disp('***************************************************************');
    switch class(state)
        case ('traj')
            log.disp('Saving trajectories to mat-files')
        case ('wave')
            log.disp('Saving wave functions to mat-files')
    end
    log.disp('***************************************************************');
    log.disp(' ');
    log.disp(['Directory name     : ' log.shorten(state.sav_dir)  ])
    log.disp(['File name template : ' log.shorten(state.sav_file) ])
    if ~isempty(state.sav_step)
        log.disp(['Maximum steps per file    : ' int2str(state.sav_step) ])
    end
    if state.sav_mem>2^30
        log.disp(['Maximum file size  : ' num2str(state.sav_mem/2^30) ' GB' ])
    elseif state.sav_mem>2^20
        log.disp(['Maximum file size  : ' num2str(state.sav_mem/2^20) ' MB' ])
    elseif state.sav_mem>2^10
        log.disp(['Maximum file size  : ' num2str(state.sav_mem/2^10) ' kB' ])
    else
        log.disp(['Maximum file size  : ' int2str(state.sav_mem     ) ' B' ])
    end
    
    % Memory per time step:
    switch class(state)
        case ('traj') % position (8 byte), momentum (8 byte), channel (4 byte)
            mem_per_step = state.n_p * ( space.n_dim * 16 + 4 );
        case ('wave') % double complex = 16 byte for each grid point, channel
            mem_per_step = numel(space.dvr{1}) * 16 * hamilt.coupling.n_eqs;
    end
    if mem_per_step>2^30
        log.disp(['Memory per step    : ' int2str(mem_per_step/2^30) ' GB' ])
    elseif mem_per_step>2^20
        log.disp(['Memory per step    : ' int2str(mem_per_step/2^20) ' MB' ])
    elseif mem_per_step>2^10
        log.disp(['Memory per step    : ' int2str(mem_per_step/2^10) ' kB' ])
    else
        log.disp(['Memory per step    : ' int2str(mem_per_step     ) ' B' ])
    end

    % Determine the stepsize;
   steps_per_file = min(int32(ceil(state.sav_mem / mem_per_step)), ...
        time.steps.m_number);
    
    % Optionally: User defined step size
    if isempty(state.sav_step) || state.sav_step > steps_per_file
        state.sav_step = steps_per_file;
    end
    log.disp(['Steps per mat-file : ' int2str(state.sav_step) ])
    log.disp(' ');
    
    % Create cell arrays for states to be saved 
    switch class(state)
        case ('traj')
            arr_cha = cell(1,                     state.sav_step);
            arr_pos = cell(space.n_dim,           state.sav_step);
            arr_mom = cell(space.n_dim,           state.sav_step);
        case ('wave')
            arr_dvr = cell(hamilt.coupling.n_eqs, state.sav_step);
    end
    
    % Initialize counting of states
    index = 1;
    
end


%% Save the states and increase index
switch class(state)
    case ('traj')
        arr_cha{1,index} = state.cha;
        for k=1:space.n_dim
            arr_pos{k,index} = state.pos{k};
            arr_mom{k,index} = state.mom{k};
        end
    case ('wave')
        for m=1:hamilt.coupling.n_eqs
            arr_dvr{m,index} = state.dvr{m};
        end
end
index = index + 1;

%% Write out the WF and possible additional data.
if index > state.sav_step || step == time.steps.m_number
	if ~isdir(state.sav_dir)
		mkdir(state.sav_dir);
    end

    % Get file name (containing index)
	file_indx = int32(floor( (step-1)/double(state.sav_step) ))+1; 
    file_name = strcat(state.sav_file, '_', int2str(file_indx), '.mat');
    log.disp('***************************************************************');
    switch class(state)
        case ('traj')
            log.disp('Saving trajectories to mat-file')
        case ('wave')
            log.disp('Saving wave functions to mat-file')
    end
    log.disp('***************************************************************');
    log.disp(' ');
    log.disp(['Directory : ' log.shorten(state.sav_dir) ])
    log.disp(['File name : ' log.shorten(file_name) ])
    log.disp(' ');
    
    % Save arrays to data file and reset them afterwards
    switch class(state)
        case ('traj')
            save(fullfile(state.sav_dir,file_name), 'arr_pos', 'arr_mom', 'arr_cha', '-mat');
            arr_cha = cell(1,                     state.sav_step);
            arr_pos = cell(space.n_dim,           state.sav_step);
            arr_mom = cell(space.n_dim,           state.sav_step);
        case ('wave')
            save(file_name, 'arr_dvr', '-mat');
            arr_dvr = cell(hamilt.coupling.n_eqs, state.sav_step);
    end
    index = 1;
    
end

% Last time step only
if step == time.steps.m_number
    file_name = strcat(state.sav_file, '_0.mat');
    log.disp('***************************************************************');
    log.disp('Saving general information to mat-file:')
    log.disp('Expectation values, Hamiltonian, discretizations etc')
    log.disp('***************************************************************');
    log.disp(' ');
    log.disp(['Directory : ' log.shorten(state.sav_dir) ])
    log.disp(['File name : ' log.shorten(file_name) ])
    log.disp(' ');
    
    % Apparently, Matlab cannot save class properties
    sav_step = state.sav_step;
    save(fullfile(state.sav_dir, file_name), 'sav_step', ...
        'expect', 'hamilt', ...
        'info', 'plots', ...
        'space', 'time', '-mat');
end

end
