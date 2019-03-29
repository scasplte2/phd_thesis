%--------------------------------------------------------------------------
%
% General information: Names of program, user, host, etc
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-20xy Burkhard Schmidt
%               2007-2014 Ulf Lorenz
%
% see the README file for license details.

function init (filename)
global info

% Path and file name
[pathstr, filestr, ~] = fileparts (filename);

% Program name and version number
info.package = 'WavePacket';
info.program = filestr;

% Get SHA-1 hash value: WavePacket installation
workdir = pwd;
cd (pathstr)
sha = log.getSHA1 ('./');
if ~isempty(sha)
    info.version1 = sha;
else
    info.version1 = 'V6.0.1';
end
cd (workdir)

% Get SHA-1 hash value: Current working directory
sha = log.getSHA1 ('./');
if ~isempty(sha)
    info.version2 = sha;
else
    info.version2 = '(unversioned)';
end

% Set/unset these two lines before/after release of a new version
info.version1 = 'V6.0.2';
info.version2 = 'V6.0.2';

% Get MATLAB release
info.release = version ('-release');

% Get user name 
result = license('inuse', 'matlab');
info.user_name = result.user;
% info.user_name = result.feature;

% Get host name
if isunix 
    [status,result] = unix('hostname');
    if ~status
        info.host_name = strcat(result,' (Unix)');
    else
        info.host_name = 'Running on Unix';
    end
elseif ispc
    [status,result] = dos('hostname');
    if ~status
        info.host_name = strcat(result,' (Windows)');
    else
        info.host_name = 'Running on Windows';
    end
end

    
% Get path and file name 
info.path_name = pwd;
addpath(info.path_name);

% Open log file (or create a new one) for writing in text mode
% Discard existing contents, if any.
filename = fullfile(info.path_name, strcat(info.program, '.log'));
fclose('all'); % needed when qm_xyz is called inside a loop
info.stdout = fopen(filename, 'wt');

if info.stdout == -1
    error(strcat('Could not open log file ',filename));
end

% Output program name/version etc
log.disp ( ' ' )
log.disp ( '***************************************************************')
log.disp ( 'About this WavePacket simulation' )
log.disp ( '***************************************************************')
log.disp ( ' ' )
log.disp (['Program package : ', info.package] )
log.disp (['Program name    : ', info.program] )
log.disp (['SHA1 hash value : ', num2str(info.version1),' (MATLAB)'])
log.disp (['SHA1 hash value : ', num2str(info.version2),' (work dir.)'])
log.disp ( ' ' )
log.disp (['Path name       : ', info.path_name] )
log.disp (['User name       : ', info.user_name] )
log.disp (['Host name       : ', info.host_name] )
log.disp (['MATLAB release  : ', info.release  ] )
log.disp ( ' ' )

% Initialize stopwatch timer; Output clock/date/time
info.start_time = cputime;
log.clock;

% Output GPL and copyright info
log.disp ( '***************************************************************')
log.disp ( 'Licensing information ' )
log.disp ( '***************************************************************')
log.disp ( ' ' )
log.disp ( 'This program is subject to the GNU General ')
log.disp ( 'Public License v2, see http://www.gnu.org or the ')
log.disp ( 'root path for the license text')
log.disp ( ' ' )
log.disp ( ' (C) 2004-2018 Burkhard Schmidt, Ulf Lorenz, FU Berlin' )
log.disp ( '     http://sourceforge.net/projects/matlab.wavepacket.p' )
log.disp ( ' ' )
