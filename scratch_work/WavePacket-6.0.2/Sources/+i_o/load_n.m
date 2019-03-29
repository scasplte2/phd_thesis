% ------------------------------------------------------------------------
%
% Loads states (trajectories or wavefunctions) 
% previously stored with save from an external mat-file.
%
% To be done: Cache the above arrays, i.e. make them persistent.
% Upon every subsequent call check if we already cached this file, 
% otherwise load and cache it
%
% ------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2008-2010 Ulf Lorenz
%
% see the README file for license details.

function load_n ( state, step )

global hamilt space

% Calculate indices for the file suffix and the position in the file.
file_indx = int32(floor( (step-1)/double(state.sav_step) ))+1;
step_indx = step - (file_indx-1)*state.sav_step;

% Get file name (containing index)
file_name = strcat(state.sav_file, '_', int2str(file_indx), '.mat');
log.disp('***************************************************************');
switch class(state)
    case ('traj')
        log.disp('Loading trajectories from mat-file')
    case ('wave')
        log.disp('Loading wave functions from mat-file')
end
log.disp('***************************************************************');
log.disp(' ');
log.disp(['Directory/folder   : ' log.shorten(state.sav_dir)  ])
log.disp(['File name template : ' log.shorten(file_name) ])
log.disp(['Index of this step : ' int2str(step_indx) ])
log.disp(' ');

% Load arrays from data file
switch class(state)
    case ('traj')
        load(fullfile(state.sav_dir,file_name), 'arr_pos', 'arr_mom', 'arr_cha', '-mat');
        state.cha = arr_cha{1,step_indx};
        for k=1:space.n_dim
            state.pos{k} = arr_pos{k,step_indx};
            state.mom{k} = arr_mom{k,step_indx};
        end

    case ('wave')
        load(fullfile(state.sav_dir,file_name), 'arr_dvr', '-mat');
        for m=1:hamilt.coupling.n_eqs
            state.dvr{m} = arr_dvr{m,step_indx};
        end
end

end
