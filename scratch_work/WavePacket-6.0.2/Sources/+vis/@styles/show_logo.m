%--------------------------------------------------------------------------
%
% Include logos in the four corners of plot specified by figure handle H
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%
% see the README file for license details.

function show_logo (obj)
global info

% Upper left: WavePacket name, program, and GIT hash value
string1 = info.package;
string2 = info.program;
string3 = info.version1;
h = annotation('textbox');
h.Position            = [0.00 0.95 0.05 0.05];  % x y w h
h.HorizontalAlignment = 'left';
h.Interpreter         = 'none';
h.LineStyle           = 'none';
h.String              = { string1, string2, string3 };
h.FitHeightToText     = 'on';
h.FontName            = obj.f_name;
h.FontSize            = obj.f_small;
h.FontWeight          = obj.f_heavy;

% Upper right: Host name, date, and time
string1 = info.host_name;
string2 = date;
c = clock;
string3 = strcat ( num2str(c(4),'%2.2i'), ':', num2str(c(5),'%2.2i'));
h = annotation('textbox');
h.Position            = [0.70 0.95 0.30 0.05];  % x y w h
h.HorizontalAlignment =  'right';
h.Interpreter         = 'none';
h.LineStyle           = 'none';
h.String              = { string1, string2, string3 };
h.FitHeightToText     = 'on';
h.FontName            = obj.f_name;
h.FontSize            = obj.f_small;
h.FontWeight          = obj.f_heavy;

% Lower left: User name and MATLAB release
string1 = ['User: ' info.user_name];
string2 = ['MATLAB Release R' info.release];
h = annotation('textbox');
h.Position            = [0.00 0.0 0.15 0.05];  % x y w h
h.HorizontalAlignment = 'left';
h.Interpreter         = 'none';
h.LineStyle           = 'none';
h.String              = { string1, string2};
h.FitHeightToText     = 'on';
h.FontName            = obj.f_name;
h.FontSize            = obj.f_small;
h.FontWeight          = obj.f_heavy;

% Lower right: Full path name (if too long, truncate it) and GIT hash value
string1 = info.path_name;
string2 = info.version2;
if numel(string1) > 48
    string1 = log.shorten(string1);
end
h = annotation('textbox');
h.Position            = [0.70 0.0 0.30 0.05];  % x y w h
h.HorizontalAlignment = 'right';
h.Interpreter         = 'none';
h.LineStyle           = 'none';
h.String              = { string1, string2};
h.FitHeightToText     = 'on';
h.FontName            = obj.f_name;
h.FontSize            = obj.f_small;
h.FontWeight          = obj.f_heavy;

end
