%--------------------------------------------------------------------------
%
% Handling of logfile output
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-20xy Burkhard Schmidt
%               2007-2014 Ulf Lorenz
%
% see the README file for license details.

classdef logfile < handle
    
    properties (Access = public)
        
        package     % Name of this package is "WavePacket"
        program     % Name of program, e.g. qm_bound, qm_propa, ...
        
        release     % Matlab release
        version1    % SHA-1 hash value: WavePacket installation
        version2    % SHA-1 hash value: Current working directory
        
        user_name   % Name of user
        host_name   % Name of host (computer)
        path_name   % Name of path (working directory)
        
        fileID      % ID of logfile: HERE IS THE PROBLEM
        
        start_time  % Starting the clock 
        
    end
        
   methods (Access = public)
        
        % Constructor: Names of program, user, host, etc
        function obj = logfile (filename)

            % Path and file name
            [pathstr, filestr, ~] = fileparts (filename);
            
            % Program name and version number
            obj.package = 'WavePacket';
            obj.program = filestr;

            % Get SHA-1 hash value: WavePacket installation
            workdir = pwd;
            cd (pathstr)
            sha = log.getSHA1 ('./');
            if ~isempty(sha)
                obj.version1 = sha;
            else
                obj.version1 = 'V5.3.0';
            end
            cd (workdir)

            % Get SHA-1 hash value: Current working directory
            sha = log.getSHA1 ('./');
            if ~isempty(sha)
                obj.version2 = sha;
            else
                obj.version2 = '(unversioned)';
            end

            % Get MATLAB release
            obj.release = version ('-release');
            
            % Get user name
            result = license('inuse', 'matlab');
            obj.user_name = result.user;

            % Get host name
            if isunix
                [status,result] = unix('hostname');
                if ~status
                    obj.host_name = strcat(result,' (Unix)');
                else
                    obj.host_name = 'Running on Unix';
                end
            elseif ispc
                [status,result] = dos('hostname');
                if ~status
                    obj.host_name = strcat(result,' (Windows)');
                else
                    obj.host_name = 'Running on Windows';
                end
            end
    
            % Get path and file name
            obj.path_name = pwd;
            addpath(obj.path_name);
            
            % Open log file (or create a new one) for writing in text mode
            % Discard existing contents, if any.
            filename = fullfile(obj.path_name, strcat(obj.program, '.log'));
            fclose('all'); % needed when qm_xyz is called inside a loop
            obj.fileID = fopen(filename, 'wt');
            
            if obj.fileID == -1
                error(strcat('Could not open log file ',filename));
            end

            % Output program name/version etc
            i_o.disp ( ' ' )
            i_o.disp ( '***************************************************************')
            i_o.disp ( 'About this WavePacket simulation' )
            i_o.disp ( '***************************************************************')
            i_o.disp ( ' ' )
            i_o.disp (['Program package : ', obj.package] )
            i_o.disp (['Program name    : ', obj.program] )
            i_o.disp (['SHA1 hash value : ', num2str(obj.version1),' (MATLAB)'])
            i_o.disp (['SHA1 hash value : ', num2str(obj.version2),' (work dir.)'])
            i_o.disp ( ' ' )
            i_o.disp (['Path name       : ', obj.path_name] )
            i_o.disp (['User name       : ', obj.user_name] )
            i_o.disp (['Host name       : ', obj.host_name] )
            i_o.disp (['MATLAB release  : ', obj.release  ] )
            i_o.disp ( ' ' )
            
            % Initialize stopwatch timer; Output clock/date/time
            obj.start_time = cputime;
            clock (obj);

            % Output GPL and copyright info
            i_o.disp ( '***************************************************************')
            i_o.disp ( 'Licensing information ' )
            i_o.disp ( '***************************************************************')
            i_o.disp ( ' ' )
            i_o.disp ( 'This program is subject to the GNU General ')
            i_o.disp ( 'Public License v2, see http://www.gnu.org or the ')
            i_o.disp ( 'root path for the license text')
            i_o.disp ( ' ' )
            i_o.disp ( ' (C) 2004-2018 Burkhard Schmidt, Ulf Lorenz, FU Berlin' )
            i_o.disp ( '     http://sourceforge.net/projects/matlab.wavepacket.p' )
            i_o.disp ( ' ' )
            
        end
        
        function clock(obj)
            
            % CPU time elapsed
            tt = cputime-obj.start_time;
            
            % Date and time (nicely formatted)
            dd = date;
            c = clock;
            hh = num2str(c(4),'%2.2i');
            mm = num2str(c(5),'%2.2i');
            ss = num2str(c(6),'%5.2f');
            
            % Output
            i_o.disp ('***************************************************************');
            if tt<0.1
                i_o.disp ( 'Starting the clock ...' )
            else
                i_o.disp (['Elapsed seconds : ' num2str(tt,'%11.5e') ' seconds'] );
            end
            i_o.disp ('***************************************************************');
            i_o.disp ( ' ' )
            i_o.disp (['Date            : ', dd] )
            i_o.disp (['Time            : ', hh, ':', mm, ':', ss] );
            i_o.disp ( ' ' )
            
            
        end
      
   end
   
   methods (Static)
       
       % Shorten long text string: First and last 21 characters only
       function str_output = shorten( str_input )
           
           L = length(str_input);
           
           if L>45
               str1 = str_input(1:21);
               str2 = str_input(L-20:L);
               str_output = [str1,' ... ',str2];
           else
               str_output = str_input;
           end
           
       end
       
       % Write text string to log file and console
       % How can this method know the fileID if it is a static method?
       function disp (string)
           
           fprintf (obj.fileID,'%s\n',string);
           disp (string);
           
       end
       
       % Write error message to log file, console and terminate
       % How can this method know the fileID if it is a static method?
       function error (string)
           
           fprintf (obj.fileID,'%s\n',string);
           error (string);
           
       end
       
       retval = getSHA1 (pathstr)        % see extra file
       retval = getSVN  (pathstr)        % see extra file
       
   end
   
end