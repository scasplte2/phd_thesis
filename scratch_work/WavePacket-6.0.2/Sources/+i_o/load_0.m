% ------------------------------------------------------------------------
%
% Loads the setup previously stored with save from an external file. 
% The path and filename of saved calculation to be loaded are taken from
% properties of object "state", which can be of class "traj" (trajectories)
% or of class "wave" (wavefunctions) which is provided as an input
% argument. 
%
% The class property "sav_step" indicates how many time steps are stored 
% together in each of the mat-files. This information is needed, e.g. in
% load_n which actually retrieves the trajectories or wavefunctions 
% (class "traj" or "wave") from those files.
%
% If a second parameter "toggle" is supplied (boolean with true value),  
% this function will also overwrite the following global variables with   
% those from the saved calculation: expect, hamilt, info, space, time
% Global variables that are not set: plot. 
%
% ------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2017-20xy Burkhard Schmidt
%               2008-2010 Ulf Lorenz
%
% see the README file for license details.

function load_0 ( state, toggle )

global expect hamilt space time;

% Get variable "sav_step" from file "traj_0.mat" or "wave_0.mat"
file_name = strcat(state.sav_file, '_0.mat');
load(fullfile(state.sav_dir, file_name), 'sav_step', '-mat');
state.sav_step = sav_step;

% Console & logfile output
log.disp('***************************************************************');
log.disp('Loading general information from mat-file:')
log.disp('Expectation values, Hamiltonian, discretizations etc')
log.disp('***************************************************************');
log.disp(' ');
log.disp(['Directory mame     : ' log.shorten(state.sav_dir)  ])
log.disp(['File name template : ' log.shorten(file_name) ])
log.disp(['Steps per mat-file : ' int2str(state.sav_step) ])

% Optionally set global variables
if nargin>1 && toggle
    load(fullfile(state.sav_dir, file_name), ...
        'expect', 'hamilt', 'info', 'space', 'time', ...
        '-mat');
    log.disp('Setting global variables: expect, hamilt, info, space, time');
end
log.disp(' ');


end